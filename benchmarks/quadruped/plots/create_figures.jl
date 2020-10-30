import Pkg; Pkg.activate(joinpath(@__DIR__,"..","..")); Pkg.instantiate()

# using MATLAB
using JLD2
using LaTeXTabulars
using PGFPlotsX
using LaTeXStrings
include(joinpath("..","..","plotting.jl"))
# using StatsBase

@load joinpath(dirname(@__FILE__),"timing_data.jld2") altro_times osqp_times altro_times_socp ecos_times

w = 0.7
lw="1.25pt"
p = @pgf Axis(
   {
      "ymajorgrids",
      ylabel="computation time (ms)",
      xtick=[1,2,4,5],
      xticklabels=["ALTRO","OSQP","ALTRO","ECOS"],
      "x tick label style"={rotate=-45, anchor="north west"}
   },
   PGFBoxPlot(altro_times,      1, width=w, linewidth=lw, plot_outliers=false, opts={color=colors["ALTRO"]}),
   PGFBoxPlot(osqp_times,       2, width=w, linewidth=lw, plot_outliers=false, opts={color=colors["OSQP"]}),
   PGFBoxPlot(altro_times_socp, 4, width=w, linewidth=lw, plot_outliers=false, opts={color=colors["ALTRO"]}),
   PGFBoxPlot(ecos_times,       5, width=w, linewidth=lw, plot_outliers=false, opts={color=colors["ECOS"]}),
)
pgfsave(joinpath(IMAGE_DIR,"quadruped_times.tikz"), p, include_preamble=false)

xs = [1,2,4,5]
p = @pgf Axis(
   {
      "ybar",
      "ymajorgrids",
      "enlargelimits"=0.25,
      ylabel="computation time (ms)",
      xtick="data",
      "symbolic x coords"=["QP","SOCP"],
      "nodes near coords align"="{vertical}"
   },
   PlotInc({fill=colors["OSQP"]}, Coordinates(["QP","SOCP"],[mean(osqp_times),0])),
   PlotInc({fill=colors["ALTRO"]},Coordinates(["QP","SOCP"],[mean(altro_times),mean(altro_times_socp)])),
   PlotInc({fill=colors["ECOS"]}, Coordinates(["SOCP"],[mean(ecos_times)])),
)
pgfsave(joinpath(IMAGE_DIR,"quadruped_times.tikz"), p, include_preamble=false)
print_tex(p)

xs = [0,0.3]
w = 0.1 
quantile(altro_times, 0.25)
errbars = @pgf {
   "error bars/.cd", 
   "y dir=both, y explicit", 
   "error bar style"={"black", "line width"="1.5pt"},
   "error mark"="|",
   "error mark options"={scale=5, "line width"="1.5pt"}
}
p = @pgf Axis(
   {
      "width=3.5in",
      "height=4cm",
      "ybar",
      "ymajorgrids",
      "enlarge x limits"=1.0,
      ylabel="computation time (ms)",
      xtick=xs,
      xticklabels=["QP","SOCP"],
      # "x tick label style"={rotate=0, anchor="west"},
      "bar width"=w,
      "legend style"={
            
            at={"(0.1,0.9)"},
            anchor="north west"
      }
   },
   PlotInc({color=colors["OSQP"], fill=colors["OSQP"], errbars...}, 
      Coordinates(xs[1:1],[mean(osqp_times)], yerror=std(osqp_times))),
   PlotInc({color=colors["ALTRO"], fill=colors["ALTRO"], errbars...}, 
      Coordinates(xs[1:2],[mean(altro_times),mean(altro_times_socp)], yerror=std(altro_times))),
   PlotInc({color=colors["ECOS"], fill=colors["ECOS"], errbars...}, 
      Coordinates(xs[2:2],[mean(ecos_times)], yerror=std(ecos_times))),
   Legend("OSQP","ALTRO","ECOS")
)
         # "error mark options={line width=1pt, mark size=1pt}"
pgfsave(joinpath(IMAGE_DIR,"quadruped_times.tikz"), p, include_preamble=false)
print_tex(p)

disp(x) = string(round(mean(x), digits=3)) * " ms"

latex_tabular(joinpath(dirname(@__FILE__), "table.tex"),
              Tabular("lcl"),
              [Rule(:top),
               ["Full Friction Cone", "ECOS", "ALTRO"],
               Rule(:mid),
            #    ["N=10", mean(ecos_times), mean(altro_times_socp)],
               ["", disp(ecos_times), disp(altro_times_socp)],
               Rule(:top),
               ["Linearized Friction Cone", "OSQP", "ALTRO"],
               Rule(:mid),
            #    ["N=10", 2, 3],
               ["", disp(osqp_times), disp(altro_times)],
               Rule(:bottom)])