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

include("random_linear.jl")
include("random_linear_problem.jl")

##
opts = SolverOptions(
    cost_tolerance = 1e-6,
    constraint_tolerance = 1e-6,
    reset_duals = false,
    projected_newton = false
)

## Generate and solve Initial problem
Random.seed!(10)
n,m = 12,6
dt = 0.1
N = 1001
prob = gen_random_linear(n,m,N, dt)
U = [@SVector randn(m) for k = 1:N-1]
X = [@SVector zeros(n) for k = 1:N]
for k = 1:N-1
    X[k+1] = discrete_dynamics(TO.integration(prob), prob.model, X[k], U[k], (k-1)*dt, dt)
end
Z = Traj(X,U,fill(dt,N-1))
initial_trajectory!(prob, Z)
Z_track = prob.Z 

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

using PGFPlotsX
function PGFBoxPlot(x, y::Real=0, thresh=3*std(x); opts=(@pgf {}), plot_outliers=true)
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
            "very thick",
            "forget plot",
            "boxplot prepared" = {
                "draw direction"="y",
                "draw position"=y,
                "lower whisker"=lo,
                "lower quartile"=q1,
                "median"=q2,
                "upper quartile"=q3,
                "upper whisker"=up, 
                "box extend"=6,
            },
            "mark options"={scale=.5}
        }, 
        coords
    )
end
altro_avg = [mean(res[:time][:,1]) for res in results] 
osqp_avg = [mean(res[:time][:,2]) for res in results] 

shift = 4 
altro = map(zip(Ns,results)) do (N,res)
    times = res[:time][:,1]
    PGFBoxPlot(times, N-shift, plot_outliers=false, opts=@pgf {"red"})
end
osqp = map(zip(Ns,results)) do (N,res)
    times = res[:time][:,2]
    PGFBoxPlot(times, N+shift, plot_outliers=false, opts=@pgf {"blue"})
end
p = @pgf TikzPicture(
    Axis(
    {
        width="8in",
        grid="both",
        xlabel="Knot Points (N)"
    },
    altro...,
    osqp...,
    PlotInc({"red","dashed","no marks"}, Coordinates(Ns .- shift, altro_avg)),
    PlotInc({"blue","dashed","no marks"}, Coordinates(Ns .+ shift, osqp_avg)),
    Legend("ALTRO","OSQP")
))
print_tex(p)