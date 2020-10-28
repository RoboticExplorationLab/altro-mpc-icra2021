import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using JLD2
using PGFPlotsX

include(joinpath("..","plotting.jl"))
@load "rocket.jld2" res tol_comp tols

@pgf Axis(
    {
        "xmajorgrids",
        "ymajorgrids",
        ylabel="computation time (ms)",
    },
    PGFBoxPlot(res[:time][:,1], 1, opts ={color=colors.altro}, plot_outliers=false, width=1),
    PGFBoxPlot(res[:time][:,2], 2, opts = {color=colors.osqp}, plot_outliers=false, width=1),
)

p = @pgf Axis(
    {
        ylabel="computation time (ms)",
        xlabel="optimality tolerance",
        xmode="log",
        "legend style"={
            at={"(0.1,0.5)"},
            anchor="west"
        }
    },
    PlotInc({"very thick"}, Coordinates(tols, tol_comp[3,:])),
    PlotInc({"very thick"}, Coordinates(tols, tol_comp[4,:])),
    Legend("ALTRO","ECOS")
)
pgfsave(joinpath(IMAGE_DIR, "rocket_tol_comp.tikz"), p, include_preamble=false)

p = @pgf Axis(
    {
        ylabel="solver tolerance",
        xlabel="optimality tolerance",
        xmode="log",
        ymode="log",
        "legend style"={
            at={"(0.9,0.1)"},
            anchor="south east"
        }
    },
    PlotInc({"very thick"}, Coordinates(tols, tol_comp[1,:])),
    PlotInc({"very thick"}, Coordinates(tols, tol_comp[2,:])),
    Legend("ALTRO","ECOS")
)
pgfsave(joinpath(IMAGE_DIR, "rocket_solver_tol.tikz"), p, include_preamble=false)