# Plot a Trajectory from ALTRO

include("..\\src\\utils.jl")
include("..\\src\\struct_setup.jl")

using Plots
gr()

function plotALTRO_traj(solver::ALTROSolver, t_opts::TrajectoryOptions;
                            show_plot = true)
    X = states(solver)

    gr()

    plt3D = plot3d(getArrAtInd(X, 1), getArrAtInd(X, 2), getArrAtInd(X, 3),
                    #markershape = :circle, markersize = 2,
                    label = "Reference Trajectory", legend = :topleft,
                    xlabel = "X", ylabel = "Y", zlabel = "Z")

    scatter!([t_opts.x0[1]], [t_opts.x0[2]], [t_opts.x0[3]],
                            markershape = :hexagon, label = "Initial Location")
    scatter!([t_opts.xf[1]], [t_opts.xf[2]], [t_opts.xf[3]],
                            markershape = :hexagon, label = "Landing Site")

    if show_plot
        display(plt3D)
    end

    plt2D_XZ = plot(getArrAtInd(X, 1), getArrAtInd(X, 3),
                    label = "Reference Trajectory", legend = :topleft,
                    xlabel = "X", ylabel = "Z")

    if show_plot
        display(plt2D_XZ)
    end

    plt2D_YZ = plot(getArrAtInd(X, 2), getArrAtInd(X, 3),
                    label = "Reference Trajectory", legend = :topleft,
                    xlabel = "Y", ylabel = "Z")

    if show_plot
        display(plt2D_YZ)
    end

    return plt3D, plt2D_XZ, plt2D_YZ

end

function plotALTRO_controls(solver::ALTROSolver, r1::Rocket; show_plot = true)

    U = controls(solver)

    plt_ux = plot(getAngle3D.(U), label = "U Control Angle",
                        legend = :outerright)
    ylabel!("Thrust Angle (deg)")
    title!("Controls over Time")

    plt_umag = plot(norm.(U), label = "U Magnitude", legend = :outerright)
    hline!([r1.u_max], linecolor = :grey, linestyle = :dashdot,
                label = "Max Thrust")
    hline!([r1.mass * r1.grav], linecolor = :orange, linestyle = :dash,
                label = "Gravity Balance")
    xlabel!("time (s)")
    ylabel!("control (N)")

    plt_u = plot(plt_ux, plt_umag, layout = (2, 1))

    if show_plot
        display(plt_u)
    end

    theta = getAngleXY.(U)
    rad = getAngle3D.(U)

    plt_polar = plot(theta, rad, proj = :polar, legend = :none,
                    title = "Thrust Angle Over Time \n (Angles in DEG)")

    if show_plot
        display(plt_polar)
    end

    return plt_u, plt_polar
end

function plotALTRO_trajMPCLOOP(X_States, t_opts::TrajectoryOptionsWARM;
                            show_plot::Bool = true, cut_traj::Bool = false,
                            k_start::Int64 = 1)
    x = getArrAtInd(X_States, 1)
    y = getArrAtInd(X_States, 2)
    z = getArrAtInd(X_States, 3)

    if cut_traj
        N = k_start + t_opts.Horizon
    else
        N = size(t_opts.ref_traj_x, 1)
    end

    ref_x = getArrAtInd(t_opts.ref_traj_x[k_start:N], 1)
    ref_y = getArrAtInd(t_opts.ref_traj_x[k_start:N], 2)
    ref_z = getArrAtInd(t_opts.ref_traj_x[k_start:N], 3)

    gr()

    plt3D = plot3d(ref_x, ref_y, ref_z, #markershape = :circle, markersize = 2,
                    label = "Reference Trajectory", legend = :topleft,
                    xlabel = "X", ylabel = "Y", zlabel = "Z")
    plot3d!(x, y, z, label = "MPC Trajectory")


    # scatter!([t_opts.x0[1]], [t_opts.x0[2]], [t_opts.x0[3]],
    #                         markershape = :hexagon, label = "Initial Location")
    # scatter!([t_opts.xf[1]], [t_opts.xf[2]], [t_opts.xf[3]],
    #                         markershape = :hexagon, label = "Landing Site")

    if show_plot
        display(plt3D)
    end

    plt2D_XZ = plot(ref_x, ref_z,
                    label = "Reference Trajectory", legend = :topleft,
                    xlabel = "X", ylabel = "Z")
    plot!(x, z, label = "MPC Trajectory")


    if show_plot
        display(plt2D_XZ)
    end

    plt2D_YZ = plot(y, z,
                    label = "Reference Trajectory", legend = :topright,
                    xlabel = "Y", ylabel = "Z")
    plot!(y, z, label = "MPC Trajectory")

    if show_plot
        display(plt2D_YZ)
    end

    return plt3D, plt2D_XZ, plt2D_YZ

end


function plotALTRO_trajMPC(solver::ALTROSolver, t_opts::TrajectoryOptionsWARM;
                            show_plot::Bool = true, cut_traj::Bool = true,
                            k_start::Int64 = 1)
    X = states(solver)

    plotALTRO_trajMPCLOOP(X, t_opts, show_plot = show_plot, cut_traj = cut_traj,
                                k_start = k_start)

end
