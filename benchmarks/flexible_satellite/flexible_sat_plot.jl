import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using MATLAB
using JLD2

@load joinpath(dirname(@__FILE__),"flexible_satellite_data.jld2") altro_times altro_avg osqp_times osqp_avg

t = 1:1:36
mat"
figure
hold on
plot($osqp_times,'*')
plot($altro_times,'o')
plot($osqp_avg,'b')
plot($altro_avg,'r')
ylabel('Solve Time (ms)')
xlabel('MPC Steps')
legend('OSQP','ALTRO','OSQP 10 Step Average','ALTRO 10 Step Average')
%ylim([0,10])
hold off
"
using PGFPlotsX
pgfplotsx()

p = @pgf Axis(
    {
        xlabel="MPC Steps",
        ylabel="Solve Time (ms)"},
    Plot(
        {
            color="cyan",
            no_marks,
            "very thick"
        },
        Coordinates(t,osqp_times)
    ),
    PlotInc(
        {
            color="orange",
            no_marks,
            "very thick"
        },
        Coordinates(t, altro_times)
    ),
    Legend("OSQP","ALTRO")
)

pgfsave("paper/figures/c_max_convergence.tikz", p, include_preamble=false)
