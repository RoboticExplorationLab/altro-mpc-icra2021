include(joinpath(@__DIR__,"..","..","plotting.jl"))
colors = (altro=colorant"red", ecos=colorant"blue")

## Generate Plots
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
    ecos = map(zip(Ns,results)) do (N,res)
        times = res[:time][:,2]
        PGFBoxPlot(times, N+shift, plot_outliers=false, width=width,
            opts=@pgf {color=colors.ecos}
        )
    end
    altro_avg = [mean(res[:time][:,1]) for res in results]
    ecos_avg = [mean(res[:time][:,2]) for res in results]
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
        ecos...,
        PlotInc({"red","dashed","no marks", "very thick"}, Coordinates(Ns .- shift, altro_avg)),
        PlotInc({"blue","dashed","no marks", "very thick"}, Coordinates(Ns .+ shift, ecos_avg)),
        Legend("ALTRO","OSQP")
    ))
end
