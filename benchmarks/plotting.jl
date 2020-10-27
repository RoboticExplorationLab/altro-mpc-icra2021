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
