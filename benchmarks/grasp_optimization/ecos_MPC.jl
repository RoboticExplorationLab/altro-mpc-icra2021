using BenchmarkTools
using Convex, ECOS

include("src/grasp_model.jl") # generates model o, x0, xf

## GENERATE REFERENCE TRAJECTORY OVER FULL HORIZON
ex = 2
include("settings$ex.jl") # creates o, x0, xf

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
Convex.solve!(prob, ECOS.Optimizer(feastol=1e-4, abstol=1e-4, reltol=1e-4))

X_cold = Z.value
U_cold = [F1.value; F2.value]

## RUN MPC AND STORE SOLVE TIMES
include("src/mpc_helpers_ecos.jl")

x_curr = X_cold[:, 1] # Current state (x, v)
u_curr = U_cold[:, 1] # Controls (u) applied at this instant
hor = 10 # length of the MPC horizon in number of steps

ecos_times = []
ecos_states = []
ecos_controls = []
num_iters = N - hor - 1
push!(ecos_states, x_curr)

for iter in 1:num_iters
    global x_curr, u_curr
    local prob, b, F1, F2
    # Propagate the physics forward to the next timestep
    x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt)

    # Setup and solve
    prob, F1, F2 = ecos_mpc_setup(o, x_curr, X_cold[:, iter + hor], hor, iter)
    b = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0)) samples=3 evals=1
    Convex.solve!(prob, ECOS.Optimizer(verbose=0))
    println(prob.status)

    # extract control
    u_curr = [F1.value[:, 1]; F2.value[:, 1]]

    # printouts
    println("Iter $iter: $(round(median(b).time/1e6, digits=2)) ms")

    # push values
    append!(ecos_times, b.times/1e6)
    push!(ecos_states, x_curr)
    push!(ecos_controls, u_curr)
end
