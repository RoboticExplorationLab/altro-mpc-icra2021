
function eye(n::Int)
      return diagm(ones(n))
end
function speye(n::Int)
      sparse(eye(n))
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


N = 50

n = 4
m = 2
A = [1 0 .1 0;
     0 1 0 .1;
     0 0 1  0;
     0 0 0  1]
B = [0.005 0;
     0     0.005;
     .1    0;
     0     .1]


Q = Diagonal(10*rand(n))
R = Diagonal(0.1*ones(m))

Qf = copy(Q)

x̄ = rand(n) .+ 1
ū = 0.1*rand(m)
x0 = (rand(n) .- 1) .* x̄ * 0.5 .+3
xf = zeros(n)

dt = 0.1 # doesn't matter, just needs to be non-zero



N = N-1
Ad = copy(A)
Bd = copy(B)
nx, nu = size(Bd)
umin = -ū
umax = ū

xmin = -x̄
xmax = x̄
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
ueq = leq
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
OSQP.setup!(m; P=P, q=q, A=A, l=l, u=u, eps_abs = 1e-3, eps_rel = 1e-3, eps_prim_inf = 1e-3, eps_dual_inf = 1e-3)
# OSQP.setup!(m; P=P, q=q, A=A, l=l, u=u, eps_abs = 1e-7)
# solve
results = OSQP.solve!(m)

X_osqp,U_osqp = OSQP_postprocess_mpc(results,N,nx,nu)
