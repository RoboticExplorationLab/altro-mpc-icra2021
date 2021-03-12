import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using Altro
using TrajectoryOptimization
using RobotDynamics
import RobotZoo.LinearModels: DoubleIntegrator
const RD = RobotDynamics
const TO = TrajectoryOptimization

using JuMP
using OSQP, ECOS, SCS, COSMO, Mosek, MosekTools
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
    cost_tolerance_intermediate=1e-4,
    penalty_scaling=500.,
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
rollout!(prob)

# Solve Cold Start
altro = ALTROSolver(prob, opts, show_summary=true, verbose=1)
b1 = benchmark_solve!(altro)

# Extract reference trajectory
Z_track = TO.get_trajectory(altro)

##  Compare to ALTRO
Random.seed!(1)
prob_naive = RocketProblem(N, (N-1)*dt,
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
    contype=TrajectoryOptimization.Inequality()
)
rollout!(prob)
opts2 = SolverOptions(
    cost_tolerance_intermediate=1e-4,
    penalty_scaling=2.,
    penalty_initial=1e-2,
    verbose = 1,
    projected_newton = false,
    constraint_tolerance = 1.0e-5,
    iterations = 5000,
    iterations_inner = 100,
    iterations_linesearch = 100,
    iterations_outer = 1000,
    show_summary=true,
)
altro_naive = ALTROSolver(prob_naive, opts2)
b2 = benchmark_solve!(altro_naive)
Z_track2 = TO.get_trajectory(altro)

## Plot the solution
x_altro = vcat([Vector(RD.get_z(z)) for z in get_trajectory(altro)]...)
x = [x[1] for x in states(altro)]
y = [x[2] for x in states(altro)]
z = [x[3] for x in states(altro)]
plot(x,y,z, aspect_ratio=:equal)

# Check if the solution lies on the SOCP boundaries
U = controls(altro)
maxthrust = maximum(norm.(controls(altro)))
thrustangle = maximum([atand(norm(u[1:2])/u[3]) for u in controls(altro)])  # less than 5
glideslope = maximum([atand(norm(x[1:2])/x[3]) for x in states(altro)])     # less than 45


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
X_traj, res = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, ecos_tol=1e-4, benchmark=false)
println("Average number of iterations: ", mean(res[:iter][:,1]))
times = mean(res[:time], dims=1)
@printf("Average times: ALTRO-%0.2f ms, ECOS-%0.2f ms", times[1], times[2])
std(res[:time], dims=1)
@save joinpath(@__DIR__, "rocket.jld2") res=res

# Show actual vs reference trajectory 
@load joinpath(@__DIR__, "rocket.jld2") res
plot(states(Z_track), inds=1:3, linestyle=:dash, label="reference")
plot!(X_traj, inds=1:3, color=[1 2 3], linewidth=2, label="actual")


## Compare convergence
tols = [1e-2,1e-4,1e-6,1e-8,1e-10,1e-12]
results = map(tols) do tol
    println("Running at tolerance $tol")
    opts_mpc.constraint_tolerance = tol
    opts_mpc.cost_tolerance = tol
    opts_mpc.cost_tolerance_intermediate = tol
    opts_mpc.gradient_tolerance = tol
    opts_mpc.gradient_tolerance_intermediate = tol
    opts_mpc.show_summary = false 
    opts_mpc.show_summary = true
    optimizers = (
        JuMP.optimizer_with_attributes(ECOS.Optimizer, 
            "verbose"=>false,
            "feastol"=>opts_mpc.constraint_tolerance,
            "abstol"=>opts_mpc.cost_tolerance,
            "reltol"=>opts_mpc.cost_tolerance
        ),
        optimizer_with_attributes(COSMO.Optimizer,
            "verbose"=>false,
            "eps_abs"=>tol,
            "eps_rel"=>tol,
            "rho"=>1e-2,
            "scaling"=>0,
            "alpha"=>1.0
        ),
        optimizer_with_attributes(Mosek.Optimizer,
            "QUIET"=>true,
            "INTPNT_CO_TOL_DFEAS"=>tol,
            "INTPNT_CO_TOL_PFEAS"=>tol,
            "INTPNT_CO_TOL_MU_RED"=>tol,
            "INTPNT_CO_TOL_REL_GAP"=>tol,
            "INTPNT_QO_TOL_DFEAS"=>tol,
            "INTPNT_QO_TOL_PFEAS"=>tol,
            "INTPNT_QO_TOL_MU_RED"=>tol,
            "INTPNT_QO_TOL_REL_GAP"=>tol,
            # "INTPNT_QO_TOL_DFEAS"=>1e-4,
            # "INTPNT_QO_TOL_PFEAS"=>1e-4,
            "INTPNT_TOL_DFEAS"=>tol,
            "INTPNT_TOL_PFEAS"=>tol,
            "INTPNT_TOL_MU_RED"=>tol,
            "INTPNT_TOL_REL_GAP"=>tol,
        )
    )
    map(optimizers) do optimizer
        solvername = MOI.get(MOI.instantiate(optimizer), MOI.SolverName())
        println("  Running with $solvername")
        Random.seed!(1)
        prob_mpc = gen_tracking_problem(prob, N_mpc)
        X_traj, res, X, U = run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, num_iters=1, optimizer=optimizer, tol0=1e-6) 
        Xerr_altro = norm(X_ref - states(prob_mpc), Inf)
        Uerr_altro = norm(U_ref - controls(prob_mpc), Inf)
        Xerr_jump  = norm(X_ref - X, Inf) 
        Uerr_jump  = norm(U_ref - U, Inf)
        SA[max(Xerr_altro, Uerr_altro), max(Xerr_jump, Uerr_jump)]
    end
end
tol_comp = hcat(map(results) do res
    m = hcat(res...)
    insert(m[2,:],1,m[1])
end...)
@save joinpath(@__DIR__,"rocket.jld2") res=res tols=tols tol_comp=tol_comp

## Generate plot
using PGFPlotsX
include(joinpath("..","plotting.jl"))
@load joinpath(@__DIR__,"rocket.jld2") res tols tol_comp
legend = ("ALTRO","ECOS","COSMO","Mosek")

plots = map(enumerate(legend)) do (i,label)
    @pgf PlotInc(
        {
            color=colors[label],
            "very thick",
            mark="*",
            "mark options"={fill=colors[label]}
        },
        Coordinates(tols, tol_comp[i,:])
    )
end
p = @pgf Axis(
    {
        xlabel="solver tolerance",
        ylabel="trajectory error",
        ymode="log",
        xmode="log",
        "xmajorgrids",
        "ymajorgrids",
        "x dir"="reverse",
        "legend style"={
            "legend columns"=-1,
            at={"(0.5,-0.35)"},
            anchor="north"
        }
    },
    plots...,
    Legend(legend...)
)
pgfsave(joinpath(IMAGE_DIR,"rocket_solver_tol.tikz"),p,include_preamble=false)