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
include("rocket_landing_problem.jl")

##
x0_new = @SVector [4.0, 2.0, 20.0, -3.0, 2.0, -5.0]
xf_new = @SVector zeros(6)
N = 301
dt = 0.05
theta = 20

r1 = Rocket(10.0, [0; 0; 9.81], 20, t1.dt)
s1 = selection(USE_ALTRO, COLD)
obj1 = ObjectiveOptions(1e-2, 100.0, 1e-1)
out1 = OutputOptions(true, true, true, true, true)

prob = RocketProblem(N, (N-1)*dt, x0=x0_new, θ_max=theta, integration=RK3)
prob0 = make_problem_ALTRO_COLD(r1, obj1, t1, out1)
cost(prob)
cost(prob0)

cons = prob.constraints
cons0 = prob0.constraints
cons[1].xf ≈ cons0[1].xf
cons[2].z_min ≈ cons0[2].z_min
cons[2].z_max ≈ cons0[2].z_max
cons[3].val ≈ cons0[3].val
cons[4].A ≈ cons0[4].A
cons[4].c ≈ cons0[4].c

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
solve!(solver)

##
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
