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
using StatsPlots
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
N_mpc = 21
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
Random.seed!(1)
res = run_MPC(prob_mpc, opts, Z_track)
res[:time]
@show mean(res[:iter], dims=1)
norm(res[:err_x0],Inf)
norm(res[:err_traj],Inf)
res[:err_traj]

boxplot(res[:time][:,1], label="ALTRO", ylabel="solve time (ms)")
boxplot!(res[:time][:,2], label="OSQP")
