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

prob = gen_trajectory(n,m,N,dt)
Z_track = prob.Z
maximum(norm.(states(Z_track),Inf))
maximum(norm.(controls(Z_track),Inf))

## Generate and solve Initial problem
Random.seed!(10)

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
n,m = 30,20 
dt = 0.1
N = 1101
opts.static_bp = false
prob = gen_trajectory(n,m,N,dt)
Z_track = prob.Z
prob_mpc = gen_tracking_problem(prob, N_mpc)
res = run_MPC(prob_mpc, opts, Z_track)
@show mean(res[:iter], dims=1)
@show median(res[:time], dims=1)

##
Random.seed!(1)
n,m = 12,6 
prob = gen_trajectory(n,m,N,dt)
Z_track = prob.Z
prob_mpc = gen_tracking_problem(prob, N_mpc)
altro = ALTROSolver(prob_mpc, opts, show_summary=true, static_bp = false, save_S=true)
max_violation(altro)
cost(altro)
Z0 = deepcopy(get_trajectory(altro))
λ0 = deepcopy(Altro.get_duals(altro))
b = benchmark_solve!(altro)

##
# Update initial state by using 1st control, and adding some noise 
k_mpc = 2
x0 = discrete_dynamics(TO.integration(prob_mpc), prob_mpc.model, prob_mpc.Z[1])
x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
TO.set_initial_state!(prob_mpc, x0)

# Update tracking cost
TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

# Shift the initial trajectory
RD.shift_fill!(prob_mpc.Z)

# Shift the multipliers and penalties
Altro.shift_fill!(TO.get_constraints(altro))

solve!(solver)

## Time Horizon
Random.seed!(1)
Ns = [11,31,51,71,101]
n = 12
m = 6
prob = gen_trajectory(n, m, N, dt)
prob_mpc = gen_tracking_problem(prob, N_mpc)
Z_track = prob.Z
results = map(Ns) do N_mpc
    println("Running with $N_mpc knot points...")
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    run_MPC(prob_mpc, opts, Z_track, 100) 
end
@show [median(res[:iter], dims=1) for res in results]
@show [median(res[:time], dims=1) for res in results]
@save "horizon_comp.jld2" results Ns

## State dimension
Random.seed!(10)
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

## Control dimension
Random.seed!(15)
ms = [2,6,10,15,20,25]
n = 30 
results_m = map(ms) do m 
    println("Running with control dimension $m...")
    prob = gen_trajectory(n, m, N, dt)
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    opts.static_bp = (m < 14) 
    opts.static_bp = false
    run_MPC(prob_mpc, opts, prob.Z, 100)
end
@save "control_dim_comp.jld2" results_m ms

@show [median(res[:iter], dims=1) for res in results_m]
@show [median(res[:time], dims=1) for res in results_m]


## Generate Plots
p = comparison_plot(results, Ns, "knot points (N)")
pgfsave(joinpath(IMAGE_DIR, "horizon_comp.tikz"), p, include_preamble=false)

@load "state_dim_comp.jld2" results_n ns
p = comparison_plot(results_n, ns, "state dimension (n)", shift=1.5,width=2,ymode="log")
pgfsave(joinpath(IMAGE_DIR, "state_dim_comp.tikz"), p, include_preamble=false)

@load "control_dim_comp.jld2" results_m ms
p = comparison_plot(results_m, ms, "control dimension (m)", shift=1,width=1.5,ymode="log")
pgfsave(joinpath(IMAGE_DIR, "control_dim_comp.tikz"), p, include_preamble=false)
