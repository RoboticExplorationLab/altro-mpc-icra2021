using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Run MPC
noise = [.03*randn(6) for i=1:20] # use the same dynamics noise for both runs
include("altro_MPC.jl") # generates altro_times, altro_states, altro_controls
include("ecos_MPC.jl") # generates ecos_times, ecos_states, ecos_controls

# Verify Matching Trajectories
for i=1:num_iters
    print("Timestep $i: ")
    print("\tState diff = ", round(norm(altro_states[i+1] - ecos_states[i+1]), digits=2))
    println("\tControl diff = ", round(norm(altro_controls[i] - ecos_controls[i]), digits=2))
end

# Solve Time Difference
ave_diff = (sum(ecos_times) - sum(altro_times))/length(altro_times)
println("\n Average ALTRO solve time was $(round(ave_diff, digits=2)) ms faster than that of ECOS")

# Save Results
using JLD2
@save string(@__DIR__,"/grasp_benchmark_data.jld2") altro_times altro_states altro_controls ecos_times ecos_states ecos_controls

# Plot Timing Results
using Plots
bounds = extrema([altro_times; ecos_times])
bin_min = floor(Int, bounds[1]) - 1
bin_max = ceil(Int, bounds[2]) + 1
bins = collect(bin_min:2:bin_max)
histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
xlabel!("Solve Time (ms)")
ylabel!("Counts")
png(string(@__DIR__,"/grasp_hist.png"))
