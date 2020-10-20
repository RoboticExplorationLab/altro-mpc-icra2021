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
# using ProfileView: @profview
using StatProfilerHTML
using JLD2
using Plots

include("random_linear.jl")

##
function gen_random_linear(n,m,N)
    # Create model
    dt = 0.1 # doesn't matter, just needs to be non-zero
    A,B = gendiscrete(n,m)
    model = RD.LinearModel(A, B; dt=dt)

    # Objective
    Q = Diagonal(10*rand(n))
    R = Diagonal(0.1*ones(m))
    Qf = Q * (N-1)

    x̄ = rand(n) .+ 1
    ū = rand(m) * 10 / (N-1)
    # ū = rand(m) * 0.1 
    x0 = (rand(n) .- 1) .* x̄ * 0.5
    xf = zeros(n)
    obj = LQRObjective(Q, R, Qf, xf, N)

    # Constraints
    constraints = ConstraintList(n, m, N)
    bound = BoundConstraint(n, m, x_min=-x̄, x_max=x̄, u_min=-ū, u_max=ū)
    add_constraint!(constraints, bound, 1:N-1)
    add_constraint!(constraints, GoalConstraint(xf), N)

    # Problem
    tf = (N-1)*dt
    prob = Problem(model, obj, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)

end

function gen_OSQP_JuMP(prob::Problem)
    n,m,N = size(prob)
    A,B = Matrix(RD.get_A(prob.model)), Matrix(RD.get_B(prob.model))
    x0 = prob.x0
    xf = prob.xf
    x̄ = prob.constraints[1].z_max[1:n]
    ū = prob.constraints[1].z_max[n+1:n+m]
    dt = prob.Z[1].dt
    Q = prob.obj[1].Q * dt
    R = prob.obj[1].R * dt
    Qf = prob.obj[end].Q

    select(i, n) = (n*(i-1)+1):(n*(i-1)+n)

    jump_model = Model(OSQP.Optimizer)
    set_silent(jump_model)

    @variable(jump_model, x[1:N*n])
    @variable(jump_model, u[1:(N-1)*m])

    objective_exp = @expression(jump_model, 0.5*transpose(x[select(N, n)]) * Qf * x[select(N, n)])

    @constraint(jump_model, x[select(1, n)] .== x0)

    for i=1:N-1
        # dynamics constraints
        @constraint(jump_model, A*x[select(i, n)] + B*u[select(i, m)] .== x[select(i+1, n)])

        # stagewise state cost
        add_to_expression!(objective_exp, 0.5*transpose(x[select(i, n)]) * Q * x[select(i, n)])

        # stagewise control cost
        add_to_expression!(objective_exp, 0.5*transpose(u[select(i, m)]) * R * u[select(i, m)])

        # control/state bound constraints
        @constraint(jump_model, x[select(i, n)] .<= x̄)
        @constraint(jump_model, x[select(i, n)] .>= -x̄)
        @constraint(jump_model, u[select(i, m)] .<= ū)
        @constraint(jump_model, u[select(i, m)] .>= -ū)
    end

    @objective(jump_model, Min, objective_exp)
    return jump_model
end

"""
Generate an OSQP problem from a TrajectoryOptimization Problem
"""
function gen_OSQP(prob0::Problem, opts::SolverOptions)
    prob = copy(prob0)
    TO.add_dynamics_constraints!(prob)
    n,m,N = size(prob)
    nlp = TrajOptNLP(prob, remove_bounds=true)
    NN = N*n + (N-1)*m

    # Cost function
    TO.hess_f!(nlp)
    P = nlp.data.G
    q = zeros(NN)

    # Constraints
    TO.jac_c!(nlp)
    A = nlp.data.D
    gL, gU = TO.constraint_bounds(nlp)
    zL, zU = TO.primal_bounds!(nlp)

    # Put dynamics on top of bound constraints
    A = [A; I(NN)]
    u = [gU; zU]
    l = [gL; zL]

    model = OSQP.Model()
    OSQP.setup!(model, P=P, q=q, A=A, l=l, u=u;
        verbose=opts.verbose>0,
        eps_abs=opts.cost_tolerance,
        eps_rel=opts.cost_tolerance,
        eps_prim_inf=opts.constraint_tolerance,
        eps_dual_inf=opts.constraint_tolerance,
    )
    return model, l,u
end

"""
Generate a set of random initial conditions
"""
function gen_ICs(prob, iters=10)
    n,m,N = size(prob)
    x̄ = prob.constraints[1].z_max[1:n]
    [(rand(n) .- 1) .* x̄ * 0.5 for i = 1:iters]
end

function MPC_Altro(solver, ICs)
    times = Float64[]
    iters = Int[]
    for ic in ICs
        t = @elapsed begin
            TO.set_initial_state!(solver, ic)
            solve!(solver)
        end
        push!(times, t)
        push!(iters, iterations(solver))
    end
    println("Average iters: ", mean(iters))
    return times
end

function MPC_OSQP(model, l, u, x0_l, x0_u, ICs)
    times = Float64[]
    for ic in ICs
        t = @elapsed begin
            x0_l .= ic 
            x0_u .= ic 
            OSQP.update!(model, u=u, l=l)
            results = OSQP.solve!(model)
        end
        push!(times, t)
    end
    return times
end

"""
Compare ALTRO and OSQP on a randomly-generated linear problem of the given size.

Run an MPC controller that generates `steps` initial conditions and solves the problem,
warm-starting from the previous initial condition. Reports the median times of the MPC 
iterations.
"""
function run_comparison(n,m,N,steps=100; opts=SolverOptions())
    # Generate the problem
    prob = gen_random_linear(n,m,N)

    # Convert to OSQP
    osqp,l,u = gen_OSQP(prob, opts) 
    NN = N*n + (N-1)*m
    xi = vcat([(k-1)*(n+m) .+ (1:n) for k = 1:N]...)
    ui = vcat([(k-1)*(n+m) + n .+ (1:m) for k = 1:N-1]...)
    x0_l = view(l, (N-1)*n .+ (1:n))
    x0_u = view(u, (N-1)*n .+ (1:n))

    # Solve the first time
    altro = ALTROSolver(prob, opts)
    solve!(altro)
    res = OSQP.solve!(osqp)

    # Compare the results
    # println("Difference in Cost: ", abs(res.info.obj_val - cost(altro)))

    X_altro = vcat(Vector.(states(altro))...)
    X_osqp = res.x[xi] 
    # println("Difference in states: ", norm(X_altro - X_osqp, Inf))

    U_altro = vcat(Vector.(controls(altro))...)
    U_osqp = res.x[ui] 
    # println("Difference in controls: ", norm(U_altro - U_osqp, Inf))

    # Generate the initial conditions 
    ICs = gen_ICs(prob, steps)

    # Change the ALTRO solver options for MPC
    set_options!(altro, reset_duals=false, penalty_initial=1e-1, penalty_scaling=1000., 
        show_summary=false, verbose=0)

    # Run MPC
    t_altro = MPC_Altro(altro, ICs)
    t_osqp = MPC_OSQP(osqp, l, u, x0_l, x0_u, ICs)

    println("Median ALTRO time: ", median(t_altro))
    println("Median OSQP time:  ", median(t_osqp))
    return median(t_altro), median(t_osqp)
end

function comp_plot(xs, times_altro, times_osqp; kwargs...)
    times_altro *= 1000
    times_osqp *= 1000
    avg_altro = mean.(eachrow(times_altro))
    std_altro = std.(eachrow(times_altro))
    avg_osqp = mean.(eachrow(times_osqp))
    std_osqp = std.(eachrow(times_osqp))
    p = plot(ylabel="time (ms)"; kwargs...) 
    plot!(xs, avg_altro, yerr=std_altro, markerstrokecolor=:auto, label="ALTRO")
    plot!(xs, avg_osqp, yerr=std_osqp, markerstrokecolor=:auto, label="OSQP")
    return p
end

opts = SolverOptions(
    cost_tolerance = 1e-4,
    constraint_tolerance = 1e-4,
    projected_newton = false
)

## Horizon Length comparison 
Random.seed!(10)
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
