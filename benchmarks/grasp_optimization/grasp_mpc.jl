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

function run_grasp_mpc(prob_mpc, opts_mpc, Z_track, num_iters = 20)
    n, m, N_mpc = size(prob_mpc)
    x0 = state(Z_track[1])
    u0 = control(Z_track[1])

    # Arrays for results
    altro_times = zeros(num_iters)
    altro_iters = zeros(num_iters)
    altro_states = [zero(x0) for i = 1:num_iters+1]
    altro_states[1] = x0
    altro_controls = [zero(u0) for i = 1:num_iters]
    ecos_times = zeros(num_iters)
    ecos_controls = copy(altro_controls)

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

    for iter in 1:num_iters
        # Updates prob_mpc in place, returns an equivalent ecos problem
        prob_mpc_ecos, X_ecos, U_ecos = mpc_update!(prob_mpc, o, iter, Z_track)

        # TODO Shift the multipliers and penalties
        # Altro.shift_fill!(TO.get_constraints(altro))

        # Solve Altro
        Altro.solve!(altro)

        # Solve Ecos
        # Convex.solve!(prob_mpc_ecos, ecos)
        JuMP.set_optimizer(prob_mpc_ecos, ()->ecos)
        JuMP.optimize!(prob_mpc_ecos)

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

        # Update arrays
        altro_times[iter] = altro.stats.tsolve
        altro_iters[iter] = iterations(altro)
        altro_states[iter+1] = state(prob_mpc.Z[1])
        altro_controls[iter] = control(prob_mpc.Z[1])
        ecos_times[iter] = 1000*ecos.sol.solve_time
        ecos_controls[iter] = value.(U_ecos)[:, 1] #U_ecos.value[:, 1]
    end

    altro_traj = Dict(:states=>altro_states, :controls=>altro_controls)
    res = Dict(:time=>[altro_times ecos_times], :iter=>altro_iters)

    # Print Solve Time Difference
    ave_diff = (sum(ecos_times) - sum(altro_times))/length(altro_times)
    println("\n Average ALTRO solve time was $(round(ave_diff, digits=2)) ms faster than that of ECOS")

    return res, altro_traj, ecos_controls
end

## Test single run
res, altro_traj, ecos_controls = run_grasp_mpc(prob_mpc, opts, Z_track, num_iters)

## Plot Timing Results
# using Plots
# altro_times = res[:time][:,1]
# ecos_times = res[:time][:,2]
# bounds = extrema([altro_times; ecos_times])
# bin_min = floor(Int, bounds[1]) - 1
# bin_max = ceil(Int, bounds[2]) + 1
# bins = collect(bin_min:bin_max)
# histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
# histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
# xlabel!("Solve Time (ms)")
# ylabel!("Counts")
