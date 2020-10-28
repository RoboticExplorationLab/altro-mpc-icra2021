import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using Altro
using TrajectoryOptimization
using RobotDynamics
import RobotZoo.LinearModels: DoubleIntegrator
const RD = RobotDynamics
const TO = TrajectoryOptimization

using JuMP
using OSQP, ECOS, SCS, COSMO
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
using Printf
using StatProfilerHTML
using JLD2
using Plots

include("../mpc.jl")
include("simple_rocket.jl")
include("rocket_landing_problem.jl")

## Cold Start
x0_new = @SVector [4.0, 2.0, 20.0, -3.0, 2.0, -5.0]
xf_new = @SVector zeros(6)
N = 301
dt = 0.05
theta = 5 # deg
glide = 45 # deg

opts = SolverOptions(
    cost_tolerance_intermediate=1e-2,
    penalty_scaling=10.,
    penalty_initial=1e-2,
    # verbose = 1,
    projected_newton = false,
    constraint_tolerance = 1.0e-5,
    iterations = 5000,
    iterations_inner = 100,
    iterations_linesearch = 100,
    iterations_outer = 500,
)
prob = RocketProblem(N, (N-1)*dt,
    Qfk=1e4,
    Rk=1e-0,
    x0=x0_new,
    θ_thrust_max=theta,
    θ_glideslope=glide,
    integration=Exponential,
    gravity=SA[0,0,-9.81],
    include_goal=true,
    include_thrust_angle=true,
    include_glideslope=true,
)
# initial_controls!(prob, zero.(controls(prob)))
rollout!(prob)
plot(states(prob), inds=1:3)

altro = ALTROSolver(prob, opts, show_summary=true, verbose=1)
Altro.solve!(altro)
Z_track = TO.get_trajectory(altro)

x_altro = vcat([Vector(RD.get_z(z)) for z in get_trajectory(altro)]...)
x = [x[1] for x in states(altro)]
y = [x[2] for x in states(altro)]
z = [x[3] for x in states(altro)]
plot(x,y,z, aspect_ratio=:equal)
plot(controls(altro), inds=1:3)
U = controls(altro)
maximum(norm.(controls(altro)))
maximum([atand(norm(u[1:2])/u[3]) for u in controls(altro)])
maximum([atand(norm(x[1:2])/x[3]) for x in states(altro)])


## Convert to tracking problem
Random.seed!(1)
opts_mpc = SolverOptions(
    cost_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    constraint_tolerance=1e-4,
    reset_duals=false,
    penalty_initial=1000.0,
    penalty_scaling=10.0,
    projected_newton=false
)
N_mpc = 21
prob_mpc = gen_tracking_problem(prob, N_mpc)
X_traj, res = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, ecos_tol=1e-4)
println("Average number of iterations: ", mean(res[:iter][:,1]))
times = mean(res[:time], dims=1)
@printf("Average times: ALTRO-%0.2f ms, ECOS-%0.2f ms", times[1], times[2])
@save "rocket.jld2" res=res

# Show actual vs reference trajectory 
plot(states(Z_track), inds=1:3, linestyle=:dash, label="reference")
plot!(X_traj, inds=1:3, color=[1 2 3], linewidth=2, label="actual")

## Compare convergence
function find_max_tolerance(prob, N_mpc, X_ref, U_ref, eps)
    Z_track = get_trajectory(prob)
    tols = 10.0 .^ range(-1, -14, length=131)

    # Find ALTRO tolerance
    tol_altro = tols[1]
    t_altro = Inf
    tol_min = 1
    tol_max = 14
    iter = 1 
    mid_prev = Inf
    while iter < 20
        mid = (tol_max + tol_min)/2
        tol = 10^(-mid)
        @show mid
        opts_mpc.cost_tolerance = tol 
        opts_mpc.cost_tolerance_intermediate = tol 
        opts_mpc.constraint_tolerance = tol
        Random.seed!(1)
        prob_mpc = gen_tracking_problem(prob, N_mpc)
        _,res, = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, num_iters=1, benchmark=false)
        err_X = norm(X_ref - states(prob_mpc), Inf)
        err_U = norm(U_ref - controls(prob_mpc), Inf)
        if err_X < eps && err_U < eps  # pass, bump up lower lower tolerance
            tol_max = mid 
        else
            tol_min = mid
        end
        if abs(mid - mid_prev) < 1e-2
            println("Number of iterations: ", res[:iter][1,1])
            tol_altro = tol
            t_altro = res[:time][1,1]
            break
        end
        mid_prev = mid
        iter += 1
    end

    # Find ECOS tolerance
    tol_ecos = tols[1] 
    t_ecos = Inf
    tol_min = 1
    tol_max = 14
    iter = 1 
    mid_prev = Inf
    while iter < 20
        mid = (tol_max + tol_min)/2
        tol = 10^(-mid)
        opts_mpc.cost_tolerance = tol 
        opts_mpc.cost_tolerance_intermediate = tol 
        opts_mpc.constraint_tolerance = tol
        Random.seed!(1)
        prob_mpc = gen_tracking_problem(prob, N_mpc)
        _,res,X,U = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, num_iters=1)
        err_X = norm(X_ref - X, Inf)
        err_U = norm(U_ref - U, Inf)
        if err_X < eps && err_U < eps  # pass, bump up lower lower tolerance
            tol_max = mid 
        else
            tol_min = mid
        end
        if abs(mid - mid_prev) < 1e-2
            tol_ecos = tol
            t_ecos = res[:time][1,2]
            break
        end
        mid_prev = mid
        iter += 1
    end
    return tol_altro, tol_ecos, t_altro, t_ecos
end

opts_mpc = SolverOptions(
    cost_tolerance=1e-8,
    cost_tolerance_intermediate=1e-8,
    constraint_tolerance=1e-8,
    reset_duals=false,
    penalty_initial=1000.0,
    penalty_scaling=10.0,
    projected_newton=false,
    show_summary=false
)
optimizer = optimizer_with_attributes(COSMO.Optimizer, 
    "verbose"=>false, "eps_abs"=>1e-12, "eps_rel"=>1e-12)

# Get reference trajectory
Random.seed!(1)
prob_mpc = gen_tracking_problem(prob, N_mpc)
X_traj, _, X_ref, U_ref = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, num_iters=1, optimizer=optimizer) 
norm(X_ref - states(prob_mpc), Inf)
norm(U_ref - controls(prob_mpc), Inf)

find_max_tolerance(prob, N_mpc, X_ref, U_ref, 1e-4)
tols = [1e-0, 1e-1, 1e-2, 1e-3, 1e-4]
tol_comp = map(tols) do tol
    SVector{4}(find_max_tolerance(prob, N_mpc, X_ref, U_ref, tol))
end
tol_comp_mat = hcat(tol_comp...)
@save "rocket.jld2" tol_comp=tol_comp_mat res=res tols=tols

plot(tols, tol_comp_mat[3,:], xscale=:log10)
plot!(tols, tol_comp_mat[4,:])


