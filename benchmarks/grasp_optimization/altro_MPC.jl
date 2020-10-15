using BenchmarkTools
using TrajectoryOptimization
const TO = TrajectoryOptimization

include("src/grasp_model.jl") # generates model o, x0, xf

## GENERATE REFERENCE TRAJECTORY OVER FULL HORIZON
ex = 2
include("settings$ex.jl") # creates o, x0, xf

# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n
u_ind = n .+ (1:m)
F1_ind = n .+ (1:Int(m/2))
F2_ind = n .+ (1+Int(m/2):m)

# objective
Q = 1.0e-3*Diagonal(@SVector ones(n))
Qf = 10.0*Diagonal(@SVector ones(n))
R = 1.0*Diagonal(@SVector ones(m))
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Goal Constraint
goal = GoalConstraint(SVector{n}(xf))
add_constraint!(conSet, goal, N)

include("src/new_constraints.jl")
for i = 1:N-1
    global model, u_ind, F1_ind, F2_ind, n, m
    local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2

    # Torque Balance
    B = [o.B[1][i] o.B[2][i]]
    t_bal = LinearConstraint(n, m, B, [θdd[i],0,0], Equality(), u_ind)
    add_constraint!(conSet, t_bal, i:i)

    # Max Grasp Force
    A = zeros(2, m)
    A[1,1:Int(m/2)] = o.v[1][i]
    A[2,1+Int(m/2):end] = o.v[2][i]
    max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
    add_constraint!(conSet, max_f, i:i)

    # SOCP friction cone
    v1_i = o.v[1][i]
    A1 = (I - v1_i*v1_i')
    c1 = o.mu*v1_i
    nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
    add_constraint!(conSet, nc1, i:i)

    v2_i = o.v[2][i]
    A2 = (I - v2_i*v2_i')
    c2 = o.mu*v2_i
    nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
    add_constraint!(conSet, nc2, i:i)
end

# Problem
prob = Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet);

u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

using Altro
# for example 2, solves in 50 ms
opts = SolverOptions(
    verbose = 1,
    projected_newton_tolerance=1e-5,
    cost_tolerance_intermediate=1e-5,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)

# normal solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
solve!(altro)

# extract results
X_cold = states(altro)
U_cold = controls(altro)

## RUN MPC
include("src/mpc_helpers_altro.jl")

# Initial Conditions
x_curr = X_cold[1] # Current state (x, v)
u_curr = U_cold[1] # Controls (u) applied at this instant
hor = 10 # length of the MPC horizon in number of steps

altro_times = []
altro_states = []
altro_controls = []
num_iters = N - hor - 1

for iter in 1:num_iters
    global x_curr, u_curr
    local X_warm, U_warm, obj, conSet, prob, altro, b

    # Propagate the physics forward to the next timestep
    x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt)

    # Construct the warm-started state and control arrays
    X_warm = [[x_curr]; X_cold[2 + iter:hor + iter]]
    U_warm = U_cold[1 + iter:hor + iter]

    # Set up problem
    obj, conSet = altro_mpc_setup(o, X_warm, U_warm, hor, iter)
    prob = Problem(o, obj, X_warm[end], dt*hor, x0=SVector{n}(X_warm[1]), constraints=conSet);

    initial_controls!(prob, U_warm)
    initial_states!(prob, X_warm)
    rollout!(prob)

    # Solve
    altro = ALTROSolver(prob, opts)
    set_options!(altro, show_summary=false, verbose=0)
    solve!(altro)
    b = benchmark_solve!(altro, samples=3, evals=1)

    # extract control
    U_new = controls(prob)
    u_curr = U_new[1]

    # printouts
    println("Iter $iter: $(round(median(b).time/1e6, digits=2)) ms")
    println("Max violation: $(TrajectoryOptimization.max_violation(altro))")

    # push values
    append!(altro_times, b.times/1e6)
    push!(altro_states, x_curr)
    push!(altro_controls, u_curr)
end
