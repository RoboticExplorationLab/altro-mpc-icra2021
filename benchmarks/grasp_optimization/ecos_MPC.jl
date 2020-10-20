using BenchmarkTools
using Convex, ECOS

include("src/grasp_model.jl") # generates model o, x0, xf

## GENERATE REFERENCE TRAJECTORY OVER FULL HORIZON
ex = 2
include("settings$ex.jl") # creates o, x0, xf

# Indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n

# Variables
F1 = Variable(Int(m/2), N-1)
F2 = Variable(Int(m/2), N-1)
Z = Variable(n, N)

# Objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
objective = Q*sumsquares(Z[:,1:N-1]) + Qf*sumsquares(Z[:,N]) + R*sumsquares([F1;F2])
prob = minimize(objective)

# Start and Goal Constraint
prob.constraints += Z[:,1] == x0
prob.constraints += Z[:,N] == xf

# Stage Constraints
for t = 1:N-1
    global o, F1, F2, prob, Z, dt, θdd

    # Torque Balance
    prob.constraints += [θdd[t],0,0] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]

    # Max Grasp Force
    prob.constraints += o.v[1][t]'*F1[:, t] <= o.f
    prob.constraints += o.v[2][t]'*F2[:, t] <= o.f

    # SOCP friction cone
    prob.constraints += norm((I - o.v[1][t]*o.v[1][t]')*F1[:, t]) <= o.mu*o.v[1][t]'*F1[:, t] # friction cone
    prob.constraints += norm((I - o.v[2][t]*o.v[2][t]')*F2[:, t]) <= o.mu*o.v[2][t]'*F2[:, t] # friction cone

    # Dynamics
    u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
    prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
    prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2
end

# Solve
Convex.solve!(prob, ECOS.Optimizer(feastol=1e-4, abstol=1e-4, reltol=1e-4))

# Extract Results
X_cold = Z.value # for setting reference X trajectory
U_cold = [F1.value; F2.value] # for setting reference U trajectory

## RUN MPC AND STORE SOLVE TIMES
include("src/mpc_helpers_ecos.jl")

# Setup
num_iters = 15 # number of MPC iterations
hor = 10 # length of the MPC horizon in number of steps
x_curr = X_cold[:, 1] # Current state (x, v)
u_curr = U_cold[:, 1] # Controls (u) applied at this instant
X_warm = X_cold[:, 1:hor] # X trajectory for warmstarting
U_warm = U_cold[:, 1:hor-1] # U trajectory for warmstarting

# Arrays for results
ecos_times = []
ecos_states = []
ecos_controls = []
push!(ecos_states, x_curr)

for iter in 1:num_iters
    global x_curr, u_curr, X_warm, U_warm
    local prob, Z, F1, F2

    # Propagate the physics forward to the next timestep
    x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt)

    # Set trajectory to track
    X_ref = [x_curr X_cold[:, iter .+ (2:hor)]]
    U_ref = U_cold[:, iter .+ (1:hor-1)]

    # Setup problem
    prob, Z, F1, F2 = ecos_mpc_setup(o, X_ref, U_ref, hor, iter, X_warm, U_warm)

    # Solve
    b = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0, feastol=1e-4, abstol=1e-4, reltol=1e-4)) samples=3 evals=1
    Convex.solve!(prob, ECOS.Optimizer(verbose=0))
    println(prob.status)

    # Extract trajectory for warmstarting
    X_warm = Z.value
    U_warm = [F1.value; F2.value]
    u_curr = U_warm[:, 1]

    # Printouts
    println("Iter $iter: $(round(median(b).time/1e6, digits=2)) ms")

    # Push values
    append!(ecos_times, b.times/1e6)
    push!(ecos_states, x_curr)
    push!(ecos_controls, u_curr)
end
