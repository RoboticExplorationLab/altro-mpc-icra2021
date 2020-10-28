include("src/grasp_model.jl")       # defines model
include("src/grasp_problem.jl")     # defines LQR problem
include("src/new_constraints.jl")   # defines additional constraints
include("src/grasp_mpc_helpers.jl") # defines mpc_update!() and gen_ECOS()
include("../mpc.jl")                # defines tracking problem

function run_grasp_mpc(prob_mpc, opts_mpc, Z_track, num_iters = 20; print_all=true)
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

    # ALTRO solver
    altro = ALTROSolver(prob_mpc, opts_mpc)
    set_options!(altro, show_summary=false, verbose=0)

    # ECOS solver
    ecos = ECOS.Optimizer(verbose=0,
                        feastol=opts_mpc.constraint_tolerance,
                        abstol=opts_mpc.cost_tolerance,
                        reltol=opts_mpc.cost_tolerance)

    for iter in 1:num_iters
        # Updates prob_mpc in place, returns an equivalent ecos problem
        prob_mpc_ecos, X_ecos, U_ecos = mpc_update!(prob_mpc, o, iter, Z_track)

        # Solve Altro
        Altro.solve!(altro)

        # Solve Ecos
        JuMP.set_optimizer(prob_mpc_ecos, ()->ecos)
        JuMP.optimize!(prob_mpc_ecos)

        # Compute max infinity norm diff
        diffs = []
        X = [value.(X_ecos[:,i]) for i in 1:N_mpc]
        U = [value.(U_ecos[:,i]) for i in 1:N_mpc-1]
        xdiff = maximum(norm.(X - states(altro), Inf))
        udiff = maximum(norm.(U - controls(altro), Inf))

        # Printouts
        if print_all
            println("Timestep $iter")
            print("ALTRO runtime: $(round(altro.stats.tsolve, digits=2)) ms")
            println("\t Max violation: $(TrajectoryOptimization.max_violation(altro))")
            print("ECOS runtime: $(round(1000*ecos.sol.solve_time, digits=2)) ms")
            println("\tStatus: ", termination_status(prob_mpc_ecos))
            println("State diff = ", round(xdiff, digits=2), "\tControl diff = ", round(udiff, digits=2))
        end

        # Update arrays
        altro_times[iter] = altro.stats.tsolve
        altro_iters[iter] = iterations(altro)
        altro_states[iter+1] = state(prob_mpc.Z[1])
        altro_controls[iter] = control(prob_mpc.Z[1])
        ecos_times[iter] = 1000*ecos.sol.solve_time
        ecos_controls[iter] = value.(U_ecos)[:, 1]
    end

    altro_traj = Dict(:states=>altro_states, :controls=>altro_controls)
    res = Dict(:time=>[altro_times ecos_times], :iter=>altro_iters)

    # Print Solve Time Difference
    ave_diff = mean(ecos_times) - mean(altro_times)
    println("\nAverage ALTRO solve time was $(round(ave_diff, digits=2)) ms faster than that of ECOS\n")

    return res, altro_traj, ecos_controls
end

## Generate and solve reference problem
o = SquareObject()
prob_cold = GraspProblem(o,251)
opts = SolverOptions(
    verbose = 0,
    projected_newton=false,
    cost_tolerance=1e-6,
    cost_tolerance_intermediate=1e-4,
    constraint_tolerance=1e-6
)
altro = ALTROSolver(prob_cold, opts, show_summary=true)
# benchmark_solve!(altro)
Altro.solve!(altro)
norm(states(altro) - X_ref,Inf)
states(altro) ≈ X_ref

# Extract Results
Z_track = get_trajectory(altro)

## MPC Setup
Random.seed!(1)
opts_mpc = SolverOptions(
    cost_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    constraint_tolerance=1e-4,
    projected_newton=false,
    penalty_initial=10.,
    penalty_scaling=100.,
    reset_duals=false,
)
num_iters = 20 # number of MPC iterations
N_mpc = 21 # length of the MPC horizon in number of steps

Q = prob_cold.obj[1].Q[1]
R = prob_cold.obj[1].R[1]
Qf = prob_cold.obj[end].Q[1]

prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)
altro = ALTROSolver(prob_mpc, opts_mpc, show_summary=true)
solve!(altro)

mpc_update!(prob_mpc, o, 1, Z_track)
solve!(altro)
X_ref_mpc ≈ states(altro)

# Test single run
num_iters = 20
print_all = true
res, altro_traj, ecos_controls = run_grasp_mpc(prob_mpc, opts_mpc, Z_track,
                                            num_iters, print_all=print_all)
mean(res[:iter])

# Histogram of timing results
altro_times = res[:time][:,1]
ecos_times = res[:time][:,2]
bounds = extrema([altro_times; ecos_times])
bin_min = floor(Int, bounds[1]) - 1
bin_max = ceil(Int, bounds[2]) + 1
bins = collect(bin_min:bin_max)
histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
xlabel!("Solve Time (ms)")
ylabel!("Counts")

# Plot of tangential to normal force ratio
U = altro_traj[:controls]
normal1 = [dot(o.v[1][i+1], U[i][1:3]) for i = 1:num_iters]
normal2 = [dot(o.v[2][i+1], U[i][4:6]) for i = 1:num_iters]
tangent1 = [norm((I - o.v[1][i+1]*o.v[1][i+1]')*U[i][1:3]) for i = 1:num_iters]
tangent2 = [norm((I - o.v[2][i+1]*o.v[2][i+1]')*U[i][4:6]) for i = 1:num_iters]
friction1 = tangent1 ./ normal1
friction2 = tangent2 ./ normal2
plot([friction1 friction2 o.mu*ones(num_iters)],
    xlabel="Time Step",
    ylabel="Tangential Force/Normal Force",
    linestyle = [:solid :solid :dash],
    label = ["F_T1/F_N1" "F_T2/F_N2" "mu"])

# Compare F1
U = altro_traj[:controls]
Ue = ecos_controls
Uc = controls(prob_cold)
N = num_iters

# altro
u1 = [U[t][1] for t = 1:N-1]
u2 = [U[t][2] for t = 1:N-1]
u3 = [U[t][3] for t = 1:N-1]
# ecos
u1e = [Ue[t][1] for t = 1:N-1]
u2e = [Ue[t][2] for t = 1:N-1]
u3e = [Ue[t][3] for t = 1:N-1]
# cold solve
u1c = [Uc[t][1] for t = 1:N-1]
u2c = [Uc[t][2] for t = 1:N-1]
u3c = [Uc[t][3] for t = 1:N-1]

plot([u1 u2 u3 u1e u2e u3e u1c u2c u3c],
    xlabel="Time Step",
    ylabel="Coordinate",
    label = ["u1" "u2" "u3" "u1e" "u2e" "u3e" "u1c" "u2c" "u3c"])
