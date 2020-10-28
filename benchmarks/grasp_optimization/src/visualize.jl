# plot square
function visualize_square(cg, θ, p, F, Fg; r = 1, xlims=(-3,7), ylims=(-4,4.5),fa=1.)
    # plot square
    corners = [[r, r], [-r, r], [-r, -r], [r, -r], [r, r]]
    rotated = [rot(θ)*corners[i]+cg for i = 1:5]
    x = [rotated[i][1] for i = 1:5]
    y = [rotated[i][2] for i = 1:5]
    Plots.plot!(x, y, seriescolor=:blue, seriesalpha=fa)

    # plot gravity
    Plots.plot!([cg[1]; cg[1]], cg[2] .+ Fg, seriescolor=:green, seriesalpha=fa)

    # plot forces
    for i = 1:length(p)
        pi = p[i]
        Fi = F[i]
        Plots.plot!((cg[1]+pi[1]) .+ [0; -Fi[1]], (cg[2]+pi[2]) .+ [0; -Fi[2]], seriescolor=:red, seriesalpha=fa)
    end

    return Plots.plot!([], legend = false, aspect_ratio=:equal, xlims=xlims, ylims=ylims)
end

function pgf_square!(plots, cg, θ, p, F, Fg; r = 1, fa=1.)
    # plot square
    corners = [[r, r], [-r, r], [-r, -r], [r, -r], [r, r]]
    rotated = [rot(θ)*corners[i]+cg for i = 1:5]
    x = [rotated[i][1] for i = 1:5]
    y = [rotated[i][2] for i = 1:5]
    push!(plots, PGFPlots.Plots.Linear(x, y, mark="none", style="solid, blue, opacity=$fa"))

    # plot gravity
    # plot!([cg[1]; cg[1]], cg[2] .+ Fg, seriescolor=:green, seriesalpha=fa)
    push!(plots, PGFPlots.Plots.Linear([cg[1]; cg[1]], cg[2] .+ Fg, mark="none", style="solid, green, opacity=$fa"))

    # plot forces
    for i = 1:length(p)
        pi = p[i]
        Fi = F[i]
        push!(plots, PGFPlots.Plots.Linear((cg[1]+pi[1]) .+ [0; -Fi[1]], (cg[2]+pi[2]) .+ [0; -Fi[2]], mark="none", style="solid, red, opacity=$fa"))

        # plot!((cg[1]+pi[1]) .+ [0; -Fi[1]], (cg[2]+pi[2]) .+ [0; -Fi[2]], seriescolor=:red, seriesalpha=fa)
    end

    # return plot!([], legend = false, aspect_ratio=:equal, xlims=xlims, ylims=ylims)
end
