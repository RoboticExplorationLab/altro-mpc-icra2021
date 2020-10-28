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
include("src/plotting.jl") # Function for comparison box plots

## Run Benchmark
println("Starting MPC Benchmark...")
Ns = [10, 15, 20, 25, 30, 35]
results = map(Ns) do N_mpc
    println("Running with $N_mpc knot points...")
    prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)
    run_grasp_mpc(prob_mpc, opts, Z_track, print_all=false)
end

# Save Results
@save string(@__DIR__,"/grasp_benchmark_data.jld2") results Ns
@load string(@__DIR__,"/grasp_benchmark_data.jld2") results Ns

# Generate Plots
timing_results = [results[i][1] for i=1:length(Ns)]
p = comparison_plot(timing_results, Ns, "knot points (N)")
pgfsave(string(@__DIR__,"/horizon_comp.tikz"), p, include_preamble=false)
