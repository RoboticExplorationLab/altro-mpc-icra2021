# plot square
function visualize_square(cg, θ, p, F, Fg; r = 1, xlims=(-4,4), ylims=(-4,4))
    # plot square
    corners = [[r, r], [-r, r], [-r, -r], [r, -r], [r, r]]
    rotated = [rot(θ)*corners[i]+cg for i = 1:5]
    x = [rotated[i][1] for i = 1:5]
    y = [rotated[i][2] for i = 1:5]
    plot(x, y, seriescolor=:blue)

    # plot gravity
    plot!([cg[1]; cg[1]], cg[2] .+ Fg, seriescolor=:green)

    # plot forces
    for i = 1:length(p)
        pi = p[i]
        Fi = F[i]
        plot!((cg[1]+pi[1]) .+ [0; -Fi[1]], (cg[2]+pi[2]) .+ [0; -Fi[2]], seriescolor=:red)
    end

    return plot!([], legend = false, aspect_ratio=:equal, xlims=xlims, ylims=ylims)
end
