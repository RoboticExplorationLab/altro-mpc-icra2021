# Test Comparison between ALTRO and ECOS

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
include(joinpath("src","utils.jl"))
include(joinpath("plotting","plot_comparison.jl"))


##
## Set up initial "cold start" problem
x0_new = @SVector [4.0, 2.0, 20.0, -3.0, 2.0, -5.0]
xf_new = @SVector zeros(6)
N = 301
dt = 0.05
theta = 20 # deg
glide = 50 # deg

opts = SolverOptions(
    cost_tolerance_intermediate=1e-2,
    penalty_scaling=10.,
    penalty_initial=1.0,
    verbose = 1,
    projected_newton = false,
    constraint_tolerance = 1.0e-8,
)
prob_altro_cold = RocketProblem(N, (N-1)*dt,
                                        x0=x0_new,
                                        θ_thrust_max=theta,
                                        θ_glidescope=glide,
                                        integration=Exponential)

solver = ALTROSolver(prob_altro_cold, opts, show_summary=true)
Altro.solve!(solver)

X_altro = states(solver)
xs_altro = getArrAtInd(X_altro, 1)
ys_altro = getArrAtInd(X_altro, 2)
zs_altro = getArrAtInd(X_altro, 3)

plt_altro3d = plot3d(xs_altro, ys_altro, zs_altro, label = "ALTRO Trajectory")
display(plt_altro3d)

Random.seed!(1)
opts_mpc = SolverOptions(
    cost_tolerance = 1e-8,
    cost_tolerance_intermediate = 1e-8,
    constraint_tolerance = 1e-8,
    projected_newton = false,
    penalty_initial = 10.,
    penalty_scaling = 100.,
    verbose = 0,
    reset_duals = false,
    show_summary = false
)

ecos_optimizer = ECOS.Optimizer(
    verbose=opts_mpc.verbose > 0,
    feastol=opts_mpc.constraint_tolerance,
    abstol=opts_mpc.cost_tolerance,
    reltol=opts_mpc.cost_tolerance
)

N_mpc = 51
prob_altro_mpc = RocketProblemMPC(prob_altro_cold, N_mpc,
    Qk = 10.,
    Qfk = 100.0,
    Rk = 1e-1,
)

Z_track = prob_altro_cold.Z
ecos, X_ecos, U_ecos = gen_ECOS_Rocket(prob_altro_mpc, Z_track, 1,
                                                    verbose = true)

altro_mpc = ALTROSolver(prob_altro_mpc, opts_mpc)
Altro.solve!(altro_mpc)
Convex.solve!(ecos, ecos_optimizer)

# Show error in position
plot_3set(states(altro_mpc), X_ecos,
                        title = "Position between ALTRO and ECOS")
# Show error in velocity
plot_3set(states(altro_mpc), X_ecos, 3,
                        title = "Velocity between ALTRO and ECOS")
# # Show error in controls
plot_3set(controls(altro_mpc), U_ecos,
                        title = "Controls between ALTRO and ECOS")

# Show error in position with reference trajectory
plot_3setRef(states(altro_mpc), X_ecos, states(Z_track)[1:N_mpc],
                        title = "Position between ALTRO and ECOS")

plot_glide_angle(states(altro_mpc), X_ecos, states(Z_track))
