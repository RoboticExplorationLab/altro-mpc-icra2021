import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using BenchmarkTools
using Convex, ECOS
using TrajectoryOptimization, Altro
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics

# Setup Benchmark
include("grasp_mpc.jl")

# Run Benchmark
println("Starting MPC Benchmark:")
altro_traj, altro_res, ecos_times = run_grasp_mpc(prob_mpc, opts, Z_track)

# Solve Time Difference
altro_times = altro_res[:time]
ave_diff = (sum(ecos_times) - sum(altro_times))/length(altro_times)
println("\n Average ALTRO solve time was $(round(ave_diff, digits=2)) ms faster than that of ECOS")

## Uncomment to save results

# # Save Results
# using JLD2
# @save string(@__DIR__,"/grasp_benchmark_data.jld2") altro_times ecos_times
#
# # Plot Timing Results
# using Plots
# bounds = extrema([altro_times; ecos_times])
# bin_min = floor(Int, bounds[1]) - 1
# bin_max = ceil(Int, bounds[2]) + 1
# bins = collect(bin_min:bin_max)
# histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
# histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
# xlabel!("Solve Time (ms)")
# ylabel!("Counts")
# png(string(@__DIR__,"/grasp_benchmark_hist.png"))
