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

## Generate a tracking problem
N_mpc = 21
prob_mpc = gen_tracking_problem(prob, N_mpc)
err_traj, err_x0 = run_MPC(prob_mpc, opts, Z_track)
err_x0
err_traj

altro = ALTROSolver(prob_mpc, opts, show_summary=true)
solve!(altro)
osqp,l,u = gen_OSQP(prob_mpc, opts)
res = OSQP.solve!(osqp)

res.info.obj_val
cost(solver)

xinds = [(k-1)*(n+m) .+ (1:n) for k = 1:N_mpc]
uinds = [(k-1)*(n+m) + n .+ (1:m) for k = 1:N_mpc-1]
xi = vcat(xinds...)
ui = vcat(uinds...)
x0_l = view(l, (N_mpc-1)*n .+ (1:n))
x0_u = view(u, (N_mpc-1)*n .+ (1:n))

X_altro = vcat(Vector.(states(altro))...)
X_osqp = res.x[xi] 
norm(X_altro - X_osqp, Inf)

U_altro = vcat(Vector.(controls(altro))...)
U_osqp = res.x[ui] 
norm(U_altro - U_osqp, Inf)

# Update
t0 = prob_mpc.t0
k_mpc = 1

# Shift the problem by 1 time-step
t0 += dt
k_mpc += 1
TO.set_initial_time!(prob_mpc, t0)
x0 = discrete_dynamics(TO.integration(prob), prob_mpc.model, prob_mpc.Z[1])
x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
TO.set_initial_state!(prob_mpc, x0)

# Update the tracking cost
TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

# Shift the trajectory
RD.shift_fill!(prob_mpc.Z)

# Shift the multipliers and penalties
Altro.shift_fill!(TO.get_constraints(altro))

x0_l .= x0
x0_u .= x0
NN = RD.num_vars(prob_mpc.Z)
q = zeros(NN)
obj = prob_mpc.obj
for k = 1:length(obj) - 1
    q[xinds[k]] .= obj[k].q * dt
    q[uinds[k]] .= obj[k].r * dt
end
q[xinds[end]] .= obj[end].q
OSQP.update!(osqp, q=q, l=l, u=u)

solve!(altro)
res = OSQP.solve!(osqp)

X_altro = vcat(Vector.(states(altro))...)
X_osqp = res.x[xi] 
norm(X_altro - X_osqp, Inf)
X_altro[1:n] ≈ x0
norm(X_osqp[1:n] - x0)

U_altro = vcat(Vector.(controls(altro))...)
U_osqp = res.x[ui] 
norm(U_altro - U_osqp, Inf)
 
## Horizon Length comparison 
t_altro, t_osqp = run_comparison(12,4,101, 1000, opts=opts)

Ns = [11,21,31,41,51,101]
n_runs = 5
times_altro = zeros(length(Ns), n_runs)
times_osqp = zeros(length(Ns), n_runs)
for (i,N) in enumerate(Ns)
    for j = 1:n_runs
        times_altro[i,j], times_osqp[i,j] = run_comparison(12,4,N, opts=opts)
    end
end
@save joinpath(@__DIR__,"horizon_comp.jld2") times_altro times_osqp
@load joinpath(@__DIR__,"horizon_comp.jld2") times_altro times_osqp
comp_plot(Ns, times_altro, times_osqp, xlabel="horizon length")
savefig(joinpath(@__DIR__, "horizon_comp.png"))

## State Dimension
m = 2
ns = [2,5,10,15,20,30,40,50]
N = 51
n_runs = 5
times_altro = zeros(length(ns), n_runs)
times_osqp = zeros(length(ns), n_runs)
for (i,n) in enumerate(ns)
    println("Size dim $n")
    opts.static_bp = n < 40
    for j = 1:n_runs
        times_altro[i,j], times_osqp[i,j] = run_comparison(n,m,N, opts=opts)
    end
end
@save joinpath(@__DIR__, "state_dim_comp.jld2") times_altro times_osqp
@load joinpath(@__DIR__, "state_dim_comp.jld2") times_altro times_osqp
comp_plot(ns, times_altro, times_osqp, xlabel="state dimension")
savefig(joinpath(@__DIR__, "state_dim_comp.png"))

## Control Dimension
n = 20 
ms = [2,5,10,15,20]
N = 51
n_runs = 5
times_altro = zeros(length(ms), n_runs)
times_osqp = zeros(length(ms), n_runs)
for (i,m) in enumerate(ms)
    println("Size dim $m")
    opts.static_bp = n < 40
    for j = 1:n_runs
        times_altro[i,j], times_osqp[i,j] = run_comparison(n,m,N, opts=opts)
    end
end
@save joinpath(@__DIR__, "control_dim_comp.jld2") times_altro times_osqp
@load joinpath(@__DIR__, "control_dim_comp.jld2") times_altro times_osqp
comp_plot(ms, times_altro, times_osqp, xlabel="control dimension")
savefig(joinpath(@__DIR__, "control_dim_comp.png"))
