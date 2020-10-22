import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()
using Altro
using TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
const TO = TrajectoryOptimization

using JuMP
using OSQP
using SparseArrays
using BenchmarkTools
using TimerOutputs
using Random
using Profile
using Statistics
using LinearAlgebra
using StaticArrays
# using ProfileView: @profview
using StatProfilerHTML
using JLD2
using Plots

include("../plotting.jl")
include("random_linear.jl")
include("random_linear_problem.jl")

##
function gen_trajectory(n,m,N,dt)
    prob = gen_random_linear(n,m,N, dt)
    U = [@SVector randn(m) for k = 1:N-1]
    X = [@SVector zeros(n) for k = 1:N]
    for k = 1:N-1
        X[k+1] = discrete_dynamics(TO.integration(prob), prob.model, X[k], U[k], (k-1)*dt, dt)
    end
    Z = Traj(X,U,fill(dt,N-1))
    initial_trajectory!(prob, Z)
    prob
end

opts = SolverOptions(
    cost_tolerance = 1e-6,
    constraint_tolerance = 1e-6,
    reset_duals = false,
    projected_newton = false
)
n,m = 12,6
dt = 0.1
N = 1101
prob = gen_trajectory(n,m,N,dt)
Z_track = prob.Z

## Generate and solve Initial problem
Random.seed!(10)

## Generate the (tracking) MPC problem
prob_mpc = gen_tracking_problem(prob, N_mpc)

## Run MPC
opts = SolverOptions(
    cost_tolerance = 1e-4,
    cost_tolerance_intermediate = 1e-4,
    constraint_tolerance = 1e-4,
    penalty_initial = 1_000.,
    penalty_scaling = 100.,
    reset_duals = false,
    projected_newton = false
)
N_mpc = 21
Random.seed!(1)
res = run_MPC(prob_mpc, opts, Z_track)
res[:time]
@show mean(res[:iter], dims=1)
norm(res[:err_x0],Inf)
norm(res[:err_traj],Inf)
res[:err_traj]

## Time Horizon
Ns = [11,31,51,71,101]
results = map(Ns) do N_mpc
    println("Running with $N_mpc knot points...")
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    run_MPC(prob_mpc, opts, Z_track, length(Z_track) - Ns[end])
end
@save "horizon_comp.jld2" results Ns
@load "horizon_comp.jld2" results Ns

## State dimension
ns = [2,10,15,20,30,50]
m = 2
results_n = map(ns) do n 
    println("Running with state dimension $n...")
    prob = gen_trajectory(n, m, N, dt)
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    opts.static_bp = (n <= 30)
    run_MPC(prob_mpc, opts, prob.Z, 100)
end
@save "state_dim_comp.jld2" results_n ns
@load "state_dim_comp.jld2" results_n ns

## Control dimension
ms = [2,6,10,15,20]
n = 20 
results_m = map(ms) do m 
    println("Running with control dimension $m...")
    prob = gen_trajectory(n, m, N, dt)
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    opts.static_bp = true 
    run_MPC(prob_mpc, opts, prob.Z, 100)
end
@save "control_dim_comp.jld2" results_m ms
@load "control_dim_comp.jld2" results_m ms

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
p = comparison_plot(results, Ns, "knot points (N)")
pgfsave(joinpath(IMAGE_DIR, "horizon_comp.tikz"), p, include_preamble=false)

p = comparison_plot(results_n, ns, "state dimension (n)", shift=1.5,width=2,ymode="log")
pgfsave(joinpath(IMAGE_DIR, "state_dim_comp.tikz"), p, include_preamble=false)

p = comparison_plot(results_m, ms, "control dimension (m)", shift=1,width=1.5,ymode="log")
pgfsave(joinpath(IMAGE_DIR, "control_dim_comp.tikz"), p, include_preamble=false)