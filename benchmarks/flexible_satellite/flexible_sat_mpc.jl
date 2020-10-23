import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()


using SparseArrays, Random
using TrajectoryOptimization
using Altro
using Statistics
using LinearAlgebra
using MATLAB
using OSQP
using RobotDynamics
using RobotZoo
using StaticArrays
const RD = RobotDynamics

function eye(n::Int)
    return diagm(ones(n))
end
function speye(n::Int)
    sparse(eye(n))
end
function mat_from_vec(a)
    "Turn a vector of vectors into a matrix"


    rows = length(a[1])
    columns = length(a)
    A = zeros(rows,columns)

    for i = 1:columns
        A[:,i] = a[i]
    end

    return A
end
function OSQP_postprocess_mpc(results,N,nx,nu)
    X = zeros(nx,N+1)
    U = zeros(nu,N)
    x_vec = results.x[1:(N+1)*nx]
    u_vec = results.x[(1+(N+1)*nx):end]
    for i = 1:N+1
        index_first = 1 + (i-1)*nx
        index_last = index_first + (nx-1)
        X[:,i] = x_vec[index_first:index_last]
    end
    for i = 1:N
        index_first = 1 + (i-1)*nu
        index_last = index_first + (nu-1)
        U[:,i] = u_vec[index_first:index_last]
    end
    return X, U
end

function c2d(A,B,dt)
    n = size(A,1)
    p = size(B,2)

    expAB = exp([A*dt B*dt; zeros(p,n+p)])

    A_d = expAB[1:n,1:n]
    B_d = expAB[1:n, (n+1):end ]

    return A_d, B_d
end

function generate_AB()
    # inertia matrix
    J = diagm([1;2;3])

    # reaction wheel jacobian
    B_sc = diagm(ones(3))


    # linear momentum coupling matrix
    phi = [0 1 0;
           1 0 0;
           0 .2 -.8];

    # angular momentum coupling matrix
    delta = [0 0 1;
             0 1 0;
            -.7 .1 .1]

    # store this matrix for faster computations
    T = inv(J-delta'*delta)

    j = 3; # 3 modes

    # damping and stiffness
    zeta = [.001;.001;.001]
    Delta = [.05; .2; .125] * (2*pi)

    # damping and stiffness matrices
    C = zeros(j,j)
    K = zeros(j,j)
    for i =1:j
        C[i,i] = 2*zeta[i]*Delta[i];
        K[i,i] = Delta[i]^2;
    end


               #   mrp        w                  n                       ndot
    pdot_row = [zeros(3,3) .25*eye(3)       zeros(3,j)                 zeros(3,j)];
    wdot_row = [zeros(3,3) zeros(3,3)     T*delta'*K                  T*delta'*C];
    ndot_row = [zeros(j,3) zeros(j,3)     zeros(j,j)                  eye(j)];
    nddot_row = [zeros(j,3) zeros(j,3) (-K - delta*T*delta'*K)    (-C - delta*T*delta'*C)];

    # analytical A
    A_analytical = [pdot_row;wdot_row;ndot_row;nddot_row];

    # analytical B
    B_analytical = [zeros(3,3);
              -T*B_sc;
              zeros(j,3);
              delta*T*B_sc];

    # sample time
    dt = .5;
    Ad, Bd = c2d(A_analytical,B_analytical,dt)

    return Ad, Bd
end


function test_mpc()
# Random.seed!(1234)
Ad, Bd = generate_AB()
A = copy(Ad)
B = copy(Bd)
N = 80
n,m = size(B)


Q = Diagonal(10*ones(n))
R = Diagonal(0.1*ones(m))

Qf = copy(Q)


ū = .01*ones(m)

xf = zeros(n)
x0 = [.1;.1;.1;zeros(n-3)]
dt = 0.1 # doesn't matter, just needs to be non-zero
model = RD.LinearModel(A, B;dt = .1)
objective = LQRObjective(Q, R, Qf, xf, N)

constraints = ConstraintList(n, m, N)
bound = BoundConstraint(n, m, u_min=-ū, u_max=ū)
add_constraint!(constraints, bound, 1:N)

tf = (N-1)*dt

problem = Problem(model, objective, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)
solver = ALTROSolver(problem)
set_options!(solver, projected_newton=false)
solve!(solver)

b = benchmark_solve!(solver)

@show maximum(b)
@show max_violation(solver)
@show cost(solver)

X_altro = states(solver)
U_altro = controls(solver)



println("OSQP")

#
N = N-1
Ad = copy(A)
Bd = copy(B)
nx, nu = size(Bd)
umin = -ū
umax = ū

xmin = -Inf*ones(nx)
xmax = Inf*ones(nx)
# Objective function
QN = Qf
xr = xf

# cast mpc problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
Q_part = kron(speye(N),sparse(Q))
R_part = kron(speye(N),sparse(R))
# quadratic objective
P = blockdiag(Q_part,sparse(QN),R_part)
# linear objective
q = ([kron(ones(N),-Q*xr); -QN*xr; zeros(N*nu)])
# linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diagm(-1 => ones(N))), Ad)
B_first_part = diagm(-1 => ones(N))
Bu = kron(B_first_part[:,1:end-1], Bd)
Aeq = [Ax Bu]
leq = [-x0; zeros(N*nx)]
ueq = copy(leq)
# input ad state constraints
Aineq = speye((N+1)*nx + N*nu)
lineq = ([kron(ones(N+1),xmin) ; kron(ones(N),umin)])
uineq = ([kron(ones(N+1),xmax) ; kron(ones(N),umax)])
# OSQP constraints
A = [Aeq; Aineq]
l = [leq; lineq]
u = [ueq; uineq]
# create OSQP object
m = OSQP.Model()
# setup problem
OSQP.setup!(m; P=P, q=q, A=A, l=l, u=u, eps_abs = 1e-6, eps_rel = 1e-6, eps_prim_inf = 1e-6, eps_dual_inf = 1e-6)
# OSQP.setup!(m; P=P, q=q, A=A, l=l, u=u, eps_abs = 1e-7)
# solve
results = OSQP.solve!(m)

X_osqp,U_osqp = OSQP_postprocess_mpc(results,N,nx,nu)

# mat"
# figure
# hold on
# plot($X_osqp')
# hold off"
#
# X_altro = mat_from_vec(X_altro)
# mat"
# figure
# hold on
# plot($X_altro')
# hold off"
#
# @infiltrate
# error()


mpc_iterations = 45

osqp_times = zeros(mpc_iterations)
altro_times = zeros(mpc_iterations)
x0_new = copy(x0)
# Random.seed!(2)
for i = 1:mpc_iterations

      # now we compare an MPC start
      u_new = controls(solver);u_new = u_new[1]
      x0_new = Ad*x0_new + Bd*u_new + .0002*randn(length(x0))

      # ALTRO MPC start
      problem.x0 .= x0_new
      set_options!(solver, projected_newton=false,penalty_initial = 10.0)
      solve!(solver)
      altro_times[i] = solver.stats.tsolve


      # OSQP mpc start
      leq = [-x0_new; zeros(N*nx)]
      ueq = copy(leq)
      l_new = [leq; lineq]
      u_new = [ueq; uineq]

      OSQP.update!(m; l = l_new,u = u_new)
      results = OSQP.solve!(m)
      osqp_times[i] = results.info.solve_time*1e3

      # @infiltrate
      # error()

end

X_osqp,U_osqp = OSQP_postprocess_mpc(results,N,nx,nu)

return osqp_times, altro_times, X_osqp, states(solver)
end

sim_trials = 10
osqp_time_mat = zeros(45,sim_trials)
altro_time_mat = zeros(45,sim_trials)

for i = 1:sim_trials
    osqp_times, altro_times, X_osqp, X_altro = test_mpc()
    altro_times = altro_times[11:end]
    osqp_times = osqp_times[11:end]
    altro_time_mat[:,i] = altro_times
    osqp_time_mat[:,i] = osqp_times
end


osqp_quants = zeros(size(osqp_time_mat,1),3)
lo = zeros(size(osqp_time_mat,1))
up = zeros(size(osqp_time_mat,1))
osqp_med = zeros(size(osqp_time_mat,1))
altro_med = zeros(size(osqp_time_mat,1))
for i = 1:size(osqp_time_mat,1)
    # osqp_quants[i,1],osqp_quants[i,2],osqp_quants[i,3] =
    # quantile(osqp_time_mat[i,:],[0.25 0.5 0.75])
    # lo[i] = minimum(osqp_time_mat[i,:])
    # up[i] = maximum(osqp_time_mat[i,:])
    osqp_med[i] = median(osqp_time_mat[i,:])
    altro_med[i] = median(altro_time_mat[i,:])
end

function vec_from_mat(mat)
    #vector of vectors from matrix of column vectors

    s = size(mat)
    if length(s) == 3
        a,b,c = size(mat)

        V = fill(zeros(a,b),c)

        for i = 1:c
            V[i] = mat[:,:,i]
        end
    else
        a,b = size(mat)

        V = fill(zeros(a),b)

        for i = 1:b
            V[i] = mat[:,i]
        end
    end


    return V
end

# q1 = osqp_quants[:,1]
# q2 = osqp_quants[:,2]
# q3 = osqp_quants[:,3]
# coords = Coordinates([],[])
# opts=(@pgf {})
# y = 0
width = 1

Q = vec_from_mat(osqp_time_mat')
A = vec_from_mat(altro_time_mat')
idx_range = 1:4:45
Ns = 1:1:45
osqp = map(zip(Ns[idx_range],Q[idx_range])) do (N,q)
    PGFBoxPlot(q,N,2*std(q);width=width,opts=@pgf {color=colors.osqp})
end
altro = map(zip(Ns[idx_range],A[idx_range])) do (N,a)
    PGFBoxPlot(a,N,2*std(a);width=width,opts=@pgf {color=colors.altro})
end
xlabel = "MPC Steps"
# ymode = "linear"
ymode = "log"

p = @pgf TikzPicture(
        Axis(
        {
            # width="8in",
            "ymajorgrids",
            "xmajorgrids",
            xlabel=xlabel,
            ymode=ymode,
            ylabel="computation time (ms)",
            xtick=Ns,
            "legend style"={
                at={"(0.1,0.9)"},
                anchor="north west"
            }
        },
        altro...,
        osqp...,
        PlotInc({"red","dashed","no marks", "very thick"}, Coordinates(Ns, altro_med)),
        PlotInc({"blue","dashed","no marks", "very thick"}, Coordinates(Ns, osqp_med)),
        Legend("ALTRO","OSQP")
    ))

print("ya")

mat"
figure
hold on
plot($altro_time_mat,'r')
plot($osqp_time_mat,'b')
hold off
"

mat"
figure
hold on
boxplot($altro_time_mat')
hold off
"
# p =@pgf PlotInc(
#     {
#         opts...,
#         "solid",
#         "line width"="2pt",
#         "forget plot",
#         "boxplot prepared" = {
#             "draw direction"="y",
#             "draw position"=y,
#             "lower whisker"=lo,
#             "lower quartile"=q1,
#             "median"=q2,
#             "upper quartile"=q3,
#             "upper whisker"=up,
#             "box extend"=width,
#         },
#         "mark options"={scale=.5},
#     }
# )
# altro_avg = zeros(length(altro_times))
# osqp_avg = zeros(length(osqp_times))
#
# for i = 1:length(altro_avg)
#       altro_avg[i] = mean(altro_times[max(1,i-9):i])
#       osqp_avg[i] = mean(osqp_times[max(1,i-9):i])
# end


# function comp_plot(xs, times_altro, times_osqp; kwargs...)
#     times_altro *= 1000
#     times_osqp *= 1000
#     avg_altro = mean.(eachrow(times_altro))
#     std_altro = std.(eachrow(times_altro))
#     avg_osqp = mean.(eachrow(times_osqp))
#     std_osqp = std.(eachrow(times_osqp))
#     p = plot(ylabel="time (ms)"; kwargs...)
#     plot!(xs, avg_altro, yerr=std_altro, markerstrokecolor=:auto, label="ALTRO")
#     plot!(xs, avg_osqp, yerr=std_osqp, markerstrokecolor=:auto, label="OSQP")
#     return p
# end
#
# comp_plot(1:45, altro_time_mat, osqp_time_mat, xlabel="MPC Steps")
#



# using JLD2
# @save joinpath(dirname(@__FILE__),"flexible_satellite_data.jld2") altro_times altro_avg osqp_times osqp_avg
