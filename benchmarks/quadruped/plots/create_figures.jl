import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
# using MATLAB
using JLD2
using LaTeXStrings
using LaTeXTabulars
using StatsBase

@load joinpath(dirname(@__FILE__),"timing_data.jld2") altro_times osqp_times altro_times_socp ecos_times

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