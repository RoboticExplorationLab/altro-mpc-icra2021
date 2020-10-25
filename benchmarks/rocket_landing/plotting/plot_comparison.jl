# Comparison Plotting

get_err(a, e) = log10.(max.(10^(-20), norm.(a - e)))

function plot_3set(a, e, offset = 0;
                            title = "Position between ALTRO and ECOS",
                            show_plot = true,
                            show_error = true)

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
