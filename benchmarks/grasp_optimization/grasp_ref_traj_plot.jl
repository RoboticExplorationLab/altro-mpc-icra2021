using JLD2
using Plots
using PGFPlots

# plotting functions
include(joinpath(@__DIR__,"src","utils.jl"))
include(joinpath(@__DIR__,"src","visualize.jl"))

# load trajectory
@load string(@__DIR__,"/grasp_ref_traj.jld2") y z F1 F2 θ o_p o_mass o_g

## plot traj with Plots
# plot([])
# for t = 1:2:length(y)-1
#     global y, z, F1, F2, o_p, θ
#     local p, F
#     p = [o_p[1][t][2:3], o_p[2][t][2:3]]
#     F = [F1[t], F2[t]]
#     visualize_square([y[t],z[t]], θ[t], p, F, o_mass*o_g[2:3], fa=t/(length(y)-1))
# end
# plot!([])

## plot traj with PGFPlots
plots = []
for t = 1:2:length(y)-1
    global y, z, F1, F2, o_p, θ
    local p, F
    p = [o_p[1][t][2:3], o_p[2][t][2:3]]
    F = [F1[t], F2[t]]
    sq = pgf_square!(plots, [y[t],z[t]], θ[t], p, F, o_mass*o_g[2:3], fa=t/(length(y)-1))
end
p = Axis([plots[i] for i in 1:length(plots)], axisEqual=true)
