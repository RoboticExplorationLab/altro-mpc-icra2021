# Comparison Plotting

get_err(a, e) = log10.(max.(10^(-20), norm.(a - e)))

function plot_3set(a, e, offset = 0;
                            title = "Position between ALTRO and ECOS",
                            show_plot = true,
                            show_error = true,
                            theta = 70)

    var_a1 = getArrAtInd(a, 1 + offset)
    var_a2 = getArrAtInd(a, 2 + offset)
    var_a3 = getArrAtInd(a, 3 + offset)

    var_e1 = evaluate(e)[1 + offset,1:end]
    var_e2 = evaluate(e)[2 + offset,1:end]
    var_e3 = evaluate(e)[3 + offset,1:end]

    plt_compare3d = plot3d(var_a1, var_a2, var_a3, label = "ALTRO Trajectory",
                                xlabel = "X", ylabel = "Y", zlabel = "Z")
    plot3d!(var_e1, var_e2, var_e3, label = "ECOS Trajectory")
    title!(title)

    if show_plot
        display(plt_compare3d)
    end

    if show_error
        err_xs = get_err(var_a1, var_e1)
        err_ys = get_err(var_a2, var_e2)
        err_zs = get_err(var_a3, var_e3)

        plt_err = plot(err_xs, label = "x error", legend = :bottomright,
                                xlabel = "Timestep", ylabel = "Log10 of Error")
        plot!(err_ys, label = "y error")
        plot!(err_zs, label = "z error")
        title!(title)

        if show_plot
            display(plt_err)
        end

        return plt_compare3d, plt_err
    else
        return plt_compare3d
    end

end


function plot_3setRef(a, e, ref, offset = 0;
                            title = "Position between ALTRO and ECOS",
                            show_plot = true,
                            show_error = true)

    var_a1 = getArrAtInd(a, 1 + offset)
    var_a2 = getArrAtInd(a, 2 + offset)
    var_a3 = getArrAtInd(a, 3 + offset)

    var_e1 = evaluate(e)[1 + offset,1:end]
    var_e2 = evaluate(e)[2 + offset,1:end]
    var_e3 = evaluate(e)[3 + offset,1:end]

    var_ref1 = getArrAtInd(ref, 1 + offset)
    var_ref2 = getArrAtInd(ref, 2 + offset)
    var_ref3 = getArrAtInd(ref, 3 + offset)

    plt_compare3d = plot3d(var_ref1, var_ref2, var_ref3,
                                label = "Reference Trajectory",
                                xlabel = "X", ylabel = "Y", zlabel = "Z")
    plot3d!(var_a1, var_a2, var_a3, label = "ALTRO Trajectory")
    plot3d!(var_e1, var_e2, var_e3, label = "ECOS Trajectory")
    title!(title)

    if show_plot
        display(plt_compare3d)
    end

    if show_error
        err_xs = get_err(var_a1, var_e1)
        err_ys = get_err(var_a2, var_e2)
        err_zs = get_err(var_a3, var_e3)

        plt_err = plot(err_xs, label = "x error", legend = :bottomright,
                                xlabel = "Timestep", ylabel = "Log10 of Error")
        plot!(err_ys, label = "y error")
        plot!(err_zs, label = "z error")
        title!(title)

        if show_plot
            display(plt_err)
        end

        return plt_compare3d, plt_err
    else
        return plt_compare3d
    end

end

slope(θ, deg=true) = deg ? (sind(90 - θ) / cosd(90 - θ)) : (sin(pi/2 - θ) / cos(pi/2 - θ))

function plot_glide_angle(a, e, ref, offset = 0;
                            title = "Glidescope for ALTRO and ECOS",
                            show_plot = true,
                            theta = 70)
    var_a1 = getArrAtInd(a, 1 + offset)
    var_a2 = getArrAtInd(a, 2 + offset)
    var_a3 = getArrAtInd(a, 3 + offset)

    var_e1 = evaluate(e)[1 + offset,1:end]
    var_e2 = evaluate(e)[2 + offset,1:end]
    var_e3 = evaluate(e)[3 + offset,1:end]

    var_ref1 = getArrAtInd(ref, 1 + offset)
    var_ref2 = getArrAtInd(ref, 2 + offset)
    var_ref3 = getArrAtInd(ref, 3 + offset)

    var_a_lat = [norm([var_a1[i]; var_a2[i]]) for i in 1:length(var_a1)]
    var_e_lat = [norm([var_e1[i]; var_e2[i]]) for i in 1:length(var_e1)]
    var_ref_lat = [norm([var_ref1[i]; var_ref2[i]]) for i in 1:length(var_ref1)]

    println("Size a -> $(size(var_a_lat))")
    println("Size e -> $(size(var_e_lat))")
    println("Size ref -> $(size(var_ref_lat))")

    println("Size a -> $(size(var_a3))")
    println("Size e -> $(size(var_e3))")
    println("Size ref -> $(size(var_ref3))")

    m1 = slope(theta)
    max_x = max(maximum(var_a_lat), maximum(var_e_lat), maximum(var_ref_lat))

    plt_glide = plot(var_ref_lat, var_ref3,
                                label = "Reference Trajectory",
                                xlabel = "|XY|",
                                ylabel = "Z",
                                legend = :outerright,
                                linestyle = :dot)
    plot!(var_a_lat, var_a3, label = "ALTRO Trajectory", linestyle = :dash)
    plot!(var_e_lat, var_e3, label = "ECOS Trajectory", linestyle = :dash)
    plot!(0:max_x, x -> m1 * x, label = "Glidescope Edge")
    title!(title)

    if show_plot
        display(plt_glide)
    end

    return plt_glide
end
