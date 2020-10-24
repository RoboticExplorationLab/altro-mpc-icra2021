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
    projected_newton_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)
altro = ALTROSolver(prob_cold, opts)
Altro.solve!(altro)

# Extract Results
Z_track = get_trajectory(altro)

# MPC Setup
num_iters = 15 # number of MPC iterations
hor = 10 # length of the MPC horizon in number of steps

Q = prob_cold.obj[1].Q[1]
R = prob_cold.obj[1].R[1]
Qf = prob_cold.obj[end].Q[1]

prob_mpc = gen_tracking_problem(prob_cold, hor, Qk = Q, Rk = R, Qfk = Qf)

function run_grasp_mpc(prob_mpc, opts_mpc, Z_track, num_iters = 15)
    x0 = state(Z_track[1])
    u0 = control(Z_track[1])

    # Arrays for results
    altro_times = zeros(num_iters)
    altro_iters = zeros(num_iters)
    altro_states = [zero(x0) for i = 1:num_iters+1]
    altro_states[1] = x0
    altro_controls = [zero(u0) for i = 1:num_iters]
    ecos_times = zeros(num_iters)

    # altro solver
    altro = ALTROSolver(prob_mpc, opts_mpc)
    set_options!(altro, show_summary=false, verbose=0)

    # ecos solver
    ecos = ECOS.Optimizer(verbose=0,
                        feastol=opts_mpc.constraint_tolerance,
                        abstol=opts_mpc.cost_tolerance,
                        reltol=opts_mpc.cost_tolerance)

    for iter in 1:num_iters
        # Updates prob_mpc in place, returns an equivalent ecos problem
        prob_mpc_ecos, U_ecos = mpc_update!(prob_mpc, o, iter, Z_track)

        # TODO Shift the multipliers and penalties
        # Altro.shift_fill!(TO.get_constraints(altro))

        # Solve Altro
        Altro.solve!(altro)

        # Solve Ecos
        b_ecos = @benchmark Convex.solve!($prob_mpc_ecos, $ecos) samples=1 evals=1
        Convex.solve!(prob_mpc_ecos, ecos)

        # Printouts
        println("Timestep $iter")
        print("ALTRO runtime: $(round(altro.stats.tsolve, digits=2)) ms")
        println("\t Max violation: $(TrajectoryOptimization.max_violation(altro))")
        print("ECOS runtime: $(round(median(b_ecos).time/1e6, digits=2)) ms")
        println("\tStatus: ", prob_mpc_ecos.status)
        println("Control diff = ", round(norm(control(prob_mpc.Z[1]) - U_ecos.value[:, 1]), digits=2))

        # Update arrays
        altro_times[iter] = altro.stats.tsolve
        altro_iters[iter] = iterations(altro)
        altro_states[iter+1] = state(prob_mpc.Z[1])
        altro_controls[iter] = control(prob_mpc.Z[1])
        ecos_times[iter] = median(b_ecos).time/1e6
    end

    altro_traj = Dict(:states=>altro_states, :controls=>altro_controls)
    altro_res = Dict(:time=>altro_times, :iter=>altro_iters)
    return altro_traj, altro_res, ecos_times
end
