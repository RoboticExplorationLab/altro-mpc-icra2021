using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("altro_MPC.jl") # generates altro_times
include("ecos_MPC.jl") # generates ecos_times

using JLD2
@save string(@__DIR__,"/grasp_timing_data.jld2") altro_times ecos_times

using Plots
bins = collect(4:2:22)
histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
xlabel!("Solve Time (ms)")
ylabel!("Counts")
png(string(@__DIR__,"/grasp_hist.png"))
