import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using BenchmarkTools
using Convex, ECOS
using TrajectoryOptimization, Altro
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using JLD2

include("grasp_mpc.jl")
include("src/plotting.jl")

# Run Benchmark
println("Starting MPC Benchmark...")
Ns = [10, 15, 25]
results = map(Ns) do N_mpc
    println("Running with $N_mpc knot points...")
    prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)
    run_grasp_mpc(prob_mpc, opts, Z_track)
end

# Save Results
@save string(@__DIR__,"/grasp_benchmark_data.jld2") results Ns

# Generate Plots
# timing_results = [results[i][1] for i=1:length(Ns)]
# p = comparison_plot(timing_results, Ns, "knot points (N)")
# pgfsave(string(@__DIR__,"/horizon_comp.tikz"), p, include_preamble=false)
