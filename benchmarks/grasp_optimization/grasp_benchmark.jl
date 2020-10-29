import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using BenchmarkTools
using JuMP, ECOS, COSMO, Mosek, MosekTools
using TrajectoryOptimization, Altro
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using JLD2
using Random
using Plots

include("grasp_mpc.jl") # Run cold solve and set up MPC tracking problem
include("../plotting.jl") # Function for comparison box plots

## Run Benchmark
println("Starting MPC Benchmark...")
Random.seed!(1)
opts_cold = SolverOptions(
    verbose = 0,
    projected_newton=false,
    cost_tolerance=1e-6,
    cost_tolerance_intermediate=1e-4,
    constraint_tolerance=1e-6
)
opts_mpc = SolverOptions(
    cost_tolerance=1e-4,
    cost_tolerance_intermediate=1e-3,
    constraint_tolerance=1e-4,
    projected_newton=false,
    penalty_initial=10_000.,
    penalty_scaling=100.,
    # reset_duals = false,
)
Ns = [11, 21, 31, 41, 51]
optimizers = (
    JuMP.optimizer_with_attributes(ECOS.Optimizer, 
        "verbose"=>false,
        "feastol"=>opts_mpc.constraint_tolerance,
        "abstol"=>opts_mpc.cost_tolerance,
        "reltol"=>opts_mpc.cost_tolerance
    ),
    optimizer_with_attributes(COSMO.Optimizer,
        "verbose"=>false,
        "eps_abs"=>1e-4,
        "eps_rel"=>1e-4,
        "rho"=>1e-2,
        "scaling"=>0,
        "alpha"=>1.0
    ),
    optimizer_with_attributes(Mosek.Optimizer,
        "QUIET"=>true,
        "INTPNT_CO_TOL_DFEAS"=>1e-4,
        "INTPNT_CO_TOL_PFEAS"=>1e-4,
        "INTPNT_CO_TOL_MU_RED"=>1e-4,
        # "INTPNT_QO_TOL_DFEAS"=>1e-4,
        # "INTPNT_QO_TOL_PFEAS"=>1e-4,
        "INTPNT_TOL_DFEAS"=>1e-4,
        "INTPNT_TOL_PFEAS"=>1e-4,
        "INTPNT_TOL_MU_RED"=>1e-4,
    )
)
warmstart = Dict("Mosek"=>true, "COSMO"=>true, "ECOS"=>false)
warmstart_dual = Dict("Mosek"=>false, "COSMO"=>true, "ECOS"=>false)

results = map(optimizers) do optimizer
    solvername = MOI.get(MOI.instantiate(optimizer), MOI.SolverName())
    println("Running with $solvername")
    map(Ns) do N_mpc
        println("  Running with $N_mpc knot points...")
        # Solve Cold Start
        o = SquareObject()
        prob_cold = GraspProblem(o,251)
        altro = ALTROSolver(prob_cold, opts_cold)
        solve!(altro)
        Z_track = get_trajectory(altro)

        # Run MPC
        Q,R,Qf = 1e3,1e0,1e1
        prob_mpc = gen_tracking_problem(prob_cold, N_mpc, Qk = Q, Rk = R, Qfk = Qf)
        run_grasp_mpc(prob_mpc, opts_mpc, Z_track, print_all=false, optimizer=optimizer, 
            warmstart=warmstart[solvername], warmstart_dual=warmstart_dual[solvername]
        )
    end
end

# Save Results
@save joinpath(@__DIR__, "grasp_benchmark_data.jld2") results Ns optimizers
@load joinpath(@__DIR__, "grasp_benchmark_data.jld2") results Ns

# Generate Plots
timing_results = [results[i][1] for i=1:length(Ns)]
p = comparison_plot(timing_results, Ns, "knot points (N)", shift=0, width=4, 
    legend=("ALTRO","ECOS","COSMO","Mosek"))
pgfsave(joinpath(IMAGE_DIR,"grasp_horizon_comp.tikz"), p, include_preamble=false)


timing_results = [copy(res[1]) for res in results[1]]
for i = 2:length(optimizers)
    for j = 1:length(Ns) 
        res = results[i][j][1]
        timing_results[j][:time] = [timing_results[j][:time] res[:time][:,2]]
    end
end
timing_results[1][:time]
legend=("ALTRO","ECOS","COSMO","Mosek")

avg = hcat(map(timing_results) do res 
    vec(mean(res[:time], dims=1))
end...)
q1 = hcat(map(timing_results) do res 
    [quantile(res[:time][:,i], 0.25) for i = 1:length(legend)]
end...)
q3 = hcat(map(timing_results) do res 
    [quantile(res[:time][:,i], 0.75) for i = 1:length(legend)]
end...)

avgs = map(1:length(legend)) do i
    @pgf PlotInc(
        {
            color=colors[legend[i]],
            "very thick",
            mark="*",
            "mark options"={fill=colors[legend[i]]},
        }, 
        Coordinates(Ns, avg[i,:])
    )
end
q1s = map(1:length(legend)) do i
    @pgf PlotInc(
        {
            color=colors[legend[i]],
            "no markers",
            "dashed",
            "forget plot"
        }, 
        Coordinates(Ns, q1[i,:])
    )
end
q3s = map(1:length(legend)) do i
    @pgf PlotInc(
        {
            color=colors[legend[i]],
            "no markers",
            "dashed",
            "forget plot"
        }, 
        Coordinates(Ns, q3[i,:])
    )
end

p = @pgf Axis(
    {
        "legend style"={
            at={"(0.1,0.9)"},
            anchor="north west"
        }
    },
    q1s...,
    q3s...,
    avgs...,
    Legend(legend...)
)
pgfsave(joinpath(IMAGE_DIR,"grasp_horizon_comp.tikz"), p, include_preamble=false)