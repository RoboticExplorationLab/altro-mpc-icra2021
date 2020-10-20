import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
# using MATLAB
using JLD2

@load joinpath(dirname(@__FILE__),"timing_data.jld2") altro_times osqp_times altro_times_socp

# t = 1:1:36
# mat"
# figure
# hold on
# plot($osqp_times,'*')
# plot($altro_times,'o')
# plot($osqp_avg,'b')
# plot($altro_avg,'r')
# ylabel('Solve Time (ms)')
# xlabel('MPC Steps')
# legend('OSQP','ALTRO','OSQP 10 Step Average','ALTRO 10 Step Average')
# %ylim([0,10])
# hold off
# "

using Plots
using PGFPlotsX
pgfplotsx()

bins = collect(0:100:5000)

histogram(altro_times ./ 1000, bins=bins, fillalpha=0.33, xlabel="Solve Times (microsecond)", title="Altro vs OSQP Solve Times", label="Altro")
histogram!(osqp_times ./ 1000, bins=bins, fillalpha=0.33, label="OSQP")
histogram!(altro_times_socp ./ 1000, bins=bins, fillalpha=0.33, label="Altro w/ SOCP")
pgfsave("histogram.tikz", h)


using StatsBase

figure = @pgf Axis(
    {
        "ybar interval",
        xmajorgrids = false,
        xtick_distance="{50}",
    },
    Plot(Table(fit(Histogram, altro_times ./ 1000)))
)
pgfsave("altro_hist.tikz", figure, include_preamble=false)
