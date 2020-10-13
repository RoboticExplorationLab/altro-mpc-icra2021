using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StaticArrays, LinearAlgebra, Plots, BenchmarkTools
using RobotDynamics
using Convex, ECOS

include("src/visualize.jl")
include("src/grasp_model.jl")

ex = 2
include("example$ex.jl")

# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n

# variables
F1 = Variable(Int(m/2), N-1)
F2 = Variable(Int(m/2), N-1)
Z = Variable(n, N)

# objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
objective = Q*sumsquares(Z[:,1:N-1]) + Qf*sumsquares(Z[:,N]) + R*sumsquares([F1;F2])
prob = minimize(objective)

# start and goal constraint
prob.constraints += Z[:,1] == x0
prob.constraints += Z[:,N] == xf

# stage constraints
for t = 1:N-1
    global o, F1, F2, prob, Z, dt, θdd

    # torque balance
    prob.constraints += [θdd[t],0,0] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]

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
Convex.solve!(prob, ECOS.Optimizer)

X_cold = Z.value
U_cold = [F1.value; F2.value]

## RUN MPC
include("src/mpc.jl")

# Initial Conditions
x_curr = X_cold[:, 1] # Current state (x, v)
u_curr = X_cold[:, 1] # Controls (u) applied at this instant
hor = 10 # length of the MPC horizon in number of steps

benchmark_arr = []
MPC_Iteration_States = []
num_iters = N - hor - 1
push!(MPC_Iteration_States, x_curr)

for iter in 1:num_iters
    global x_curr, u_curr
    # Propagate the physics forward to the next timestep
    x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt)

    # Compute next control
    u_curr, b = ecos_mpc_step(o, x_curr, X_cold[:, iter + hor], hor, iter)
    println("Iter = $iter & Median Time = $(round(median(b.times) / 1e6, digits=2)) ms")

    # push values
    push!(MPC_Iteration_States, x_curr)
    push!(benchmark_arr, b)
end

# plot trajectory
x = [MPC_Iteration_States[t][1] for t = 1:num_iters]
y = [MPC_Iteration_States[t][2] for t = 1:num_iters]
z = [MPC_Iteration_States[t][3] for t = 1:num_iters]
xd = [MPC_Iteration_States[t][4] for t = 1:num_iters]
yd = [MPC_Iteration_States[t][5] for t = 1:num_iters]
zd = [MPC_Iteration_States[t][6] for t = 1:num_iters]
plot([x y z xd yd zd], label = ["x" "y" "z" "xd" "yd" "zd"])
