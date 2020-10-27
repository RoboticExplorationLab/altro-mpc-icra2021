include("src/grasp_model.jl")       # defines model
include("src/grasp_problem.jl")     # defines problem
include("src/new_constraints.jl")   # defines additional constraints
include("src/grasp_mpc_helpers.jl") # defines mpc_update! and gen_ECOS
include("../mpc.jl")                # defines gen_tracking_problem

# Generate and solve reference problem
o = SquareObject()
prob_cold = GraspProblem(o)
opts = SolverOptions(
    verbose = 0,
    projected_newton_tolerance=1e-5,
    cost_tolerance_intermediate=1e-5,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-5
)
altro = ALTROSolver(prob_cold, opts)
Altro.solve!(altro)

# Extract Results
Z_track = get_trajectory(altro)

# MPC Setup
num_iters = 1 # number of MPC iterations
N_mpc = 10 # length of the MPC horizon in number of steps

Q = prob_cold.obj[1].Q[1]
R = prob_cold.obj[1].R[1]
Qf = prob_cold.obj[end].Q[1]

prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)

n, m, N_mpc = size(prob_mpc)
x0 = state(Z_track[1])
u0 = control(Z_track[1])

# fake loop
opts_mpc = SolverOptions(
    verbose = 0,
    projected_newton_tolerance=1e-8,
    cost_tolerance_intermediate=1e-8,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-8
)
iter = 1

# altro solver
altro = ALTROSolver(prob_mpc, opts_mpc)
set_options!(altro, show_summary=false, verbose=0)

# ecos solver
ecos = ECOS.Optimizer(verbose=0,
                    feastol=1e-8,
                        abstol=1e-8,
                        reltol=1e-8)
                        # feastol=opts_mpc.constraint_tolerance,
                        # abstol=opts_mpc.cost_tolerance,
                        # reltol=opts_mpc.cost_tolerance)

# Updates prob_mpc in place, returns an equivalent ecos problem
prob_mpc_ecos, X_ecos, U_ecos = mpc_update!(prob_mpc, o, iter, Z_track)

# Solve Altro
Altro.solve!(altro)

# Solve Ecos
# Convex.solve!(prob_mpc_ecos, ecos)
JuMP.set_optimizer(prob_mpc_ecos, ()->ecos)
JuMP.optimize!(prob_mpc_ecos)

## DUMB TEST
n,m,N_mpc=size(prob_mpc)
dt = prob_mpc.Z[1].dt
# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n

# variables
F1 = Variable(Int(m/2), N_mpc-1)
F2 = Variable(Int(m/2), N_mpc-1)
Z = Variable(n, N_mpc)
X = Z
U = [F1; F2]

# objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
k_mpc=1
# Get Tracking Trajectory
X_ref = states(Z_track)[(k_mpc-1) .+ (1:N_mpc)]
X_ref = hcat(Vector.(X_ref)...)

U_ref = controls(Z_track)[(k_mpc-1) .+ (1:N_mpc-1)]
U_ref = hcat(Vector.(U_ref)...)
objective = 10^9Q*sumsquares(X[:,1:N_mpc-1]-X_ref[:, 1:N_mpc-1]) + Qf*sumsquares(X[:,N_mpc]-X_ref[:,N_mpc]) + R*sumsquares(U-U_ref)
prob = minimize(objective)

# start and goal constraint
prob.constraints += Z[:,1] == x0

# stage constraints
for t = 1:N_mpc-1
    global o, F1, F2, prob, Z, dt

    # torque balance
    prob.constraints += [o.θdd[t],0,0] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]

    # max grasp force
    prob.constraints += o.v[1][t]'*F1[:, t] <= o.f
    prob.constraints += o.v[2][t]'*F2[:, t] <= o.f

    # friction cone
    prob.constraints += norm((I - o.v[1][t]*o.v[1][t]')*F1[:, t]) <= o.mu*o.v[1][t]'*F1[:, t] # friction cone
    prob.constraints += norm((I - o.v[2][t]*o.v[2][t]')*F2[:, t]) <= o.mu*o.v[2][t]'*F2[:, t] # friction cone

    # dynamics
    u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
    prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
    prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2

end

# solve
# bm = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0))
Convex.solve!(prob, ECOS.Optimizer)
plot([X.value' X_ref'],legend=false)
## END OF DUMB TEST

# Compute max infinity norm diff
diffs = []
X = [value.(X_ecos[:,i]) for i in 1:N_mpc]
U = [value.(U_ecos[:,i]) for i in 1:N_mpc-1]
xdiff = maximum(norm.(X - states(altro), Inf))
udiff = maximum(norm.(U - controls(altro), Inf))

# Printouts
println("Timestep $iter")
print("ALTRO runtime: $(round(altro.stats.tsolve, digits=2)) ms")
println("\t Max violation: $(TrajectoryOptimization.max_violation(altro))")
print("ECOS runtime: $(round(1000*ecos.sol.solve_time, digits=2)) ms")
println("\tStatus: ", termination_status(prob_mpc_ecos)) # prob_mpc_ecos.status)
println("State diff = ", round(xdiff, digits=2), "\tControl diff = ", round(udiff, digits=2))

using Plots

X_cold = hcat(Vector.(states(prob_cold)[2:11])...)
U_cold = hcat(Vector.(controls(prob_cold)[2:10])...)
X_ecos = value.(X_ecos)
U_ecos = value.(U_ecos)
X_altro = hcat(Vector.(states(altro))...)
U_altro = hcat(Vector.(controls(altro))...)
plot([X_cold' X_ecos' X_altro'])
plot([U_cold' U_ecos' U_altro'])
