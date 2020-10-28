import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using BenchmarkTools
using JuMP, ECOS
using TrajectoryOptimization, Altro
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using JLD2
using Random
using Plots

include("grasp_mpc.jl") # Run cold solve and set up MPC tracking problem
include("../plotting.jl") # Function for comparison box plots

## Run Benchmark
println("Starting MPC Benchmark...")
Random.seed!(1)
opts_cold = SolverOptions(
    verbose = 0,
    projected_newton=false,
    cost_tolerance=1e-6,
    cost_tolerance_intermediate=1e-4,
    constraint_tolerance=1e-6
)
opts_mpc = SolverOptions(
    cost_tolerance=1e-4,
    cost_tolerance_intermediate=1e-3,
    constraint_tolerance=1e-4,
    projected_newton=false,
    penalty_initial=10_000.,
    penalty_scaling=100.,
    # reset_duals = false,
)
Ns = [11, 21, 31, 41, 51]
results = map(Ns) do N_mpc
    println("Running with $N_mpc knot points...")
    # Solve Cold Start
    o = SquareObject()
    prob_cold = GraspProblem(o,251)
    altro = ALTROSolver(prob_cold, opts_cold)
    solve!(altro)
    Z_track = get_trajectory(altro)

    # Run MPC
    Q,R,Qf = 1e3,1e0,1e1
    prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)
    run_grasp_mpc(prob_mpc, opts_mpc, Z_track, print_all=false)
end

# Save Results
@save joinpath(@__DIR__, "grasp_benchmark_data.jld2") results Ns
@load joinpath(@__DIR__, "grasp_benchmark_data.jld2") results Ns

# Generate Plots
timing_results = [results[i][1] for i=1:length(Ns)]
p = comparison_plot(timing_results, Ns, "knot points (N)", shift=0, width=4)
pgfsave(joinpath(IMAGE_DIR,"grasp_horizon_comp.tikz"), p, include_preamble=false)
