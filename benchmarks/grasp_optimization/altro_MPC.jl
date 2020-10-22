using BenchmarkTools
using TrajectoryOptimization
const TO = TrajectoryOptimization

include("src/grasp_model.jl") # generates model o, x0, xf

## GENERATE REFERENCE TRAJECTORY OVER FULL HORIZON
ex = 2
include("settings$ex.jl") # creates o, x0, xf

# Indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n
u_ind = n .+ (1:m)
F1_ind = n .+ (1:Int(m/2))
F2_ind = n .+ (1+Int(m/2):m)

# Objective
Q = 1.0e-3*Diagonal(@SVector ones(n))
Qf = 10.0*Diagonal(@SVector ones(n))
R = 1.0*Diagonal(@SVector ones(m))
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Goal Constraint
goal = GoalConstraint(SVector{n}(xf))
add_constraint!(conSet, goal, N)

# Stage Constraints
include("src/new_constraints.jl")
for i = 1:N-1
    global model, u_ind, F1_ind, F2_ind, n, m
    local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2

    # Torque Balance
    B = [o.B[1][i] o.B[2][i]]
    t_bal = LinearConstraint2(n, m, B, [θdd[i],0,0], Equality(), u_ind)
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

# Solver
using Altro
opts = SolverOptions(
    verbose = 1,
    projected_newton_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)

# Solve
altro = ALTROSolver(prob, opts)
solve!(altro)

# Extract Results
X_cold = states(altro) # for setting reference X trajectory
U_cold = controls(altro) # for setting reference U trajectory

## RUN MPC
include("src/mpc_helpers_altro.jl")

# Setup
num_iters = 15 # number of MPC iterations
hor = 10 # length of the MPC horizon in number of steps
x_curr = X_cold[1] # Current state (x, v)
u_curr = U_cold[1] # Controls (u) applied at this instant
X_warm = X_cold[1:hor] # X trajectory for warmstarting
U_warm = U_cold[1:hor-1] # U trajectory for warmstarting

obj, conSet = altro_mpc_setup(o, X_cold[1:hor], U_cold[1:hor-1], hor)
prob = Problem(o, obj, X_cold[hor], dt*(hor-1), x0=X_cold[1], constraints=conSet)
initial_controls!(prob, U_cold[1:hor-1])
rollout!(prob);

# Arrays for results
altro_times = []
altro_states = []
altro_controls = []
push!(altro_states, x_curr)

for iter in 1:num_iters
    global x_curr, u_curr, prob, altro, opts, X_warm, U_warm, noise

    # Propagate the physics forward to the next timestep
    x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt, noise[iter])

    # Set trajectory to track
    X_ref = [[x_curr]; X_cold[iter .+ (2:hor)]]
    U_ref = U_cold[iter .+ (1:hor-1)]

    # Update problem
    altro_mpc_update!(prob, o, X_ref, U_ref, hor, iter, X_warm, U_warm)

    # Solve
    altro = ALTROSolver(prob, opts)
    set_options!(altro, show_summary=false, verbose=0)
    solve!(altro)
    b = benchmark_solve!(altro, samples=3, evals=1)

    # Extract trajectory for warmstarting
    X_warm=states(altro)
    U_warm=controls(altro)
    u_curr = controls(altro)[1]

    # Printouts
    println("Iter $iter: $(round(median(b).time/1e6, digits=2)) ms")
    println("Max violation: $(TrajectoryOptimization.max_violation(altro))")

    # Push values
    append!(altro_times, b.times/1e6)
    push!(altro_states, x_curr)
    push!(altro_controls, u_curr)
end
