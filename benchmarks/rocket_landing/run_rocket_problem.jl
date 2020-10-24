import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using Altro
using TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
const TO = TrajectoryOptimization

# using JuMP
# using OSQP
using Convex
using ECOS
using SparseArrays
using BenchmarkTools
using TimerOutputs
using Random
using Profile
using Statistics
using LinearAlgebra
using StaticArrays
using StatsPlots
# using ProfileView: @profview
using StatProfilerHTML
using JLD2
using Plots
include("rocket_landing_problem.jl")
# include(joinpath("src", "struct_setup.jl"))
# include(joinpath("src","make_problem.jl"))
include(joinpath("src","convert_Altro_to_Convex.jl"))


## Set up initial "cold start" problem
x0_new = @SVector [4.0, 2.0, 20.0, -3.0, 2.0, -5.0]
xf_new = @SVector zeros(6)
N = 301
dt = 0.05
theta = 20

opts = SolverOptions(
    cost_tolerance_intermediate=1e-2,
    penalty_scaling=10.,
    penalty_initial=1.0,
    verbose = 1,
    projected_newton = false,
    constraint_tolerance = 1.0e-8,
)
prob = RocketProblem(N, (N-1)*dt, x0=x0_new, θ_max=theta, integration=Exponential)
solver = ALTROSolver(prob, opts, show_summary=true)
Altro.solve!(solver)

## Set up MPC problem
Random.seed!(1)
opts_mpc = SolverOptions(
    cost_tolerance = 1e-4,
    cost_tolerance_intermediate = 1e-4,
    constraint_tolerance = 1e-4,
    projected_newton = false,
    penalty_initial = 10.,
    penalty_scaling = 100.,
    verbose = 0,
    reset_duals = false,
    show_summary = false
)
N_mpc = 51
prob_mpc = RocketProblemMPC(prob, N_mpc,
    Qk = 10.,
    Qfk = 100.0,
    Rk = 1e-1,
)
prob_mpc.x0
Z_track = prob.Z
X, res = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track)
mean(res[:time],dims=1)[1]
mean(res[:iter],dims=1)[1]

using Plots
plot(X, inds=1:3)
plot!(states(Z_track), inds=1:3, color=[1 2 3], linestyle=:dash)
