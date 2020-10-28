using PGFPlotsX, Colors, Statistics
colors = (altro=colorant"red", osqp=colorant"blue")
cd(joinpath(@__DIR__,".."))
IMAGE_DIR = joinpath("figures")
function PGFBoxPlot(x, y::Real=0, thresh=3*std(x);
        opts=(@pgf {}),
        plot_outliers=true,
        width=6
    )
    q1,q2,q3 = quantile(x,[0.25, 0.5, 0.75])
    μ = mean(x)
    is_inlier = μ - thresh .< abs.(x) .< μ + thresh
    inliers = x[is_inlier]
    outliers = x[.!is_inlier]
    lo = minimum(inliers)
    up = maximum(inliers)
    if plot_outliers
        coords = Coordinates(fill(y,length(outliers)), outliers)
    else
        coords = Coordinates([],[])
    end

    @pgf PlotInc(
        {
            opts...,
            "solid",
            "line width"="2pt",
            "forget plot",
            "boxplot prepared" = {
                "draw direction"="y",
                "draw position"=y,
                "lower whisker"=lo,
                "lower quartile"=q1,
                "median"=q2,
                "upper quartile"=q3,
                "upper whisker"=up,
                "box extend"=width,
            },
            "mark options"={scale=.5},
        },
        coords
    )
end

function comparison_plot(results, Ns, xlabel;
        shift=4,
        width=6,
        ymode="linear"
    )
    altro = map(zip(Ns,results)) do (N,res)
        times = res[:time][:,1]
        PGFBoxPlot(times, N-shift, plot_outliers=false, width=width,
            opts=@pgf {color=colors.altro}
        )
    end
    osqp = map(zip(Ns,results)) do (N,res)
        times = res[:time][:,2]
        PGFBoxPlot(times, N+shift, plot_outliers=false, width=width,
            opts=@pgf {color=colors.osqp}
        )
    end
    altro_avg = [mean(res[:time][:,1]) for res in results] 
    osqp_avg = [mean(res[:time][:,2]) for res in results] 
    p = @pgf TikzPicture(
        Axis(
        {
            # width="8in",
            "ymajorgrids",
            "xmajorgrids",
            xlabel=xlabel,
            ymode=ymode,
            ylabel="computation time (ms)",
            xtick=Ns,
            "legend style"={
                at={"(0.1,0.9)"},
                anchor="north west"
            }
        },
        altro...,
        osqp...,
        PlotInc({"red","dashed","no marks", "very thick"}, Coordinates(Ns .- shift, altro_avg)),
        PlotInc({"blue","dashed","no marks", "very thick"}, Coordinates(Ns .+ shift, osqp_avg)),
        Legend("ALTRO","OSQP")
    ))
end