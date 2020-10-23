using BenchmarkTools
using TrajectoryOptimization
const TO = TrajectoryOptimization
using Altro

include("src/grasp_model.jl")      # defines model
include("src/new_constraints.jl")  # new constraints
include("src/mpc_helpers_altro.jl")

## GENERATE REFERENCE TRAJECTORY OVER FULL HORIZON
ex = 2
include("settings$ex.jl") # creates o, x0, xf


# Solver
opts = SolverOptions(
    verbose = 1,
    projected_newton_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)

# Solve reference problem
prob_cold = GraspProblem()
altro = ALTROSolver(prob_cold, opts)
solve!(altro)

# Extract Results
X_cold = states(altro) # for setting reference X trajectory
U_cold = controls(altro) # for setting reference U trajectory
Z_track = get_trajectory(altro)

# Setup
num_iters = 15 # number of MPC iterations
hor = 10 # length of the MPC horizon in number of steps

# MPC problem
obj, conSet = altro_mpc_setup(prob_cold.model, X_cold[1:hor], U_cold[1:hor-1], hor)
prob_mpc = Problem(o, obj, X_cold[hor], dt*(hor-1), x0=X_cold[1], constraints=conSet)
initial_controls!(prob_mpc, U_cold[1:hor-1])
rollout!(prob_mpc);

function run_grasp_mpc(prob_mpc, opts_mpc, Z_track, num_iters = 15;
        noise = [0.03*randn(6) for i = 1:num_iters] 
    )
    N_mpc = prob_mpc.N
    x_curr = state(Z_track[1])
    u_curr = control(Z_track[1])
    X_warm = states(Z_track)[1:N_mpc]
    U_warm = controls(Z_track)[1:N_mpc]

    # Arrays for results
    altro_times = zeros(num_iters) 
    iters = zeros(num_iters)
    altro_states = [zero(x_curr) for i = 1:num_iters+1] 
    altro_controls = [zero(u_curr) for i = 1:num_iters] 
    altro_states[1] = x_curr

    for iter in 1:num_iters

        # Propagate the physics forward to the next timestep
        x_curr = noisy_discrete_dynamics(o, x_curr, u_curr, dt, noise[iter])

        # Set trajectory to track
        X_ref = [[x_curr]; X_cold[iter .+ (2:hor)]]
        U_ref = U_cold[iter .+ (1:hor-1)]

        # Update problem
        altro_mpc_update!(prob_mpc, o, X_ref, U_ref, hor, iter, X_warm, U_warm)

        # Solve
        altro = ALTROSolver(prob_mpc, opts_mpc)
        set_options!(altro, show_summary=true, verbose=0)
        solve!(altro)
        # b = benchmark_solve!(altro, samples=3, evals=1)

        # Extract trajectory for warmstarting
        X_warm=states(altro)
        U_warm=controls(altro)
        u_curr = controls(altro)[1]

        # Printouts
        println("Iter $iter: $(altro.stats.tsolve), digits=2)) ms")
        println("Max violation: $(TrajectoryOptimization.max_violation(altro))")

        # Push values
        altro_times[iter] = altro.stats.tsolve
        iters[iter] = iterations(altro)
        altro_states[iter+1] = x_curr
        altro_controls[iter] = u_curr
        # append!(altro_times, b.times/1e6)
        # push!(altro_states, x_curr)
        # push!(altro_controls, u_curr)
    end
    return altro_states, Dict(:time=>altro_times, :iter=>iters)
end

run_grasp_mpc(prob_mpc, opts, Z_track)