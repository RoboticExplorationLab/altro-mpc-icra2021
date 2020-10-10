include("random_linear.jl")

n = 10
m = Int(round(n/2))
horizon_lengths = collect(10:10:100)
A,B = gendiscrete(n,m)

Q = Diagonal(10*rand(n))
R = Diagonal(0.1*ones(m))

Qf = Q #dare(A, B, Q, R)

x̄ = rand(n) .+ 1
ū = 0.1*rand(m)
x0 = (rand(n) .- 1) .* x̄ * 0.5
x0_warm = (rand(n) .- 1) .* x̄ * 0.5
xf = zeros(n)

tol = 1e-6

###############################################################################################################################
# ALTRO Benchmark Function:
using Altro
using TrajectoryOptimization
using RobotDynamics

const RD = RobotDynamics

# returns the benchmarks for ALTRO
function get_altro_benchmark(N)
    dt = 0.1 # doesn't matter, just needs to be non-zero
    model = RD.LinearModel(A, B; dt=dt)
    objective = LQRObjective(Q, R, Qf, xf, N)

    constraints = ConstraintList(n, m, N)
    bound = BoundConstraint(n, m, x_min=-x̄, x_max=x̄, u_min=-ū, u_max=ū)
    add_constraint!(constraints, bound, 1:N)

    tf = (N-1)*dt

    problem = Problem(model, objective, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)
    solver = ALTROSolver(problem)
    set_options!(solver, projected_newton=false, constraint_tolerance = tol, static_bp=(n<=14))
    solve!(solver)

    X_altro = states(solver)
    U_altro = controls(solver)

    b_no_ws = benchmark_solve!(solver)

    initial_controls!(solver, U_altro)
    initial_states!(solver, X_altro)
    set_options!(solver, reset_duals=false, penalty_initial=10.0)

    # update x0 to see warmstarted solution:
    solver.solver_al.solver_uncon.x0 .= x0_warm

    b_ws = benchmark_solve!(solver)
    @assert states(solver)[1] ≈ x0_warm

    return (b_no_ws, b_ws)
end

###############################################################################################################################
# OSQP Benchmark Function:
using ParameterJuMP, JuMP
using OSQP
using BenchmarkTools

select(i, n) = (n*(i-1)+1):(n*(i-1)+n)

function get_osqp_benchmark(N)
    jump_model = ModelWithParams(
    optimizer_with_attributes(
            OSQP.Optimizer, "eps_abs" => tol, "eps_rel" => tol, "eps_prim_inf" => tol, "eps_dual_inf" => tol 
        )
    )
    set_silent(jump_model)

    x0_param = [add_parameter(jump_model, x0[i]) for i in 1:n]

    @variable(jump_model, x[1:((N)*n)])
    @variable(jump_model, u[1:((N-1)*m)])

    objective_exp = @expression(jump_model, 0.5*transpose(x[select(N, n)]) * Qf * x[select(N, n)])

    @constraint(jump_model, initial_value_constraint, x[select(1, n)] .== x0_param)

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
        @constraint(jump_model, u[select(i, m)] .<= ū)
        @constraint(jump_model, u[select(i, m)] .>= -ū)
    end

    @objective(jump_model, Min, objective_exp)

    optimize!(jump_model)

    X_osqp = [value.(x)[select(i, n)] for i=1:N]
    U_osqp = [value.(u)[select(i, m)] for i=1:N-1]

    b_no_ws = @benchmark optimize!($jump_model) samples=10 evals=10

    # update x0 to see warmstart
    fix.(x0_param, x0_warm)

    set_start_value.(x, value.(x))
    set_start_value.(u, value.(u))

    b_ws = @benchmark optimize!($jump_model) samples=10 evals=10

    return (b_no_ws, b_ws)
end

function get_benchmarks(horizon_lengths)
    altro_no_ws = zeros(length(horizon_lengths))
    osqp_no_ws = zeros(length(horizon_lengths))

    altro_ws = zeros(length(horizon_lengths))
    osqp_ws = zeros(length(horizon_lengths))

    for i=1:length(horizon_lengths)
        N = horizon_lengths[i]
        b_altro = get_altro_benchmark(N)
        altro_no_ws[i] = median(b_altro[1].times)
        altro_ws[i] = median(b_altro[2].times)

        b_osqp = get_osqp_benchmark(N)
        osqp_no_ws[i] = median(b_osqp[1].times)
        osqp_ws[i] = median(b_osqp[2].times)
    end

    return (altro_no_ws, osqp_no_ws, altro_ws, osqp_ws)
end

(altro_no_ws, osqp_no_ws, altro_ws, osqp_ws) = get_benchmarks(horizon_lengths)

using Plots

scale = 1e-6

scatter(horizon_lengths, altro_no_ws .* scale, label="Altro", title="Random Linear Problem (no warmstart)", xlabel="Horizon Length", ylabel="Time (milliseconds)")
scatter!(horizon_lengths, osqp_no_ws .* scale, label="OSQP")
png("benchmarks/random_linear_mpc/no_warmstart_horizon")

scatter(horizon_lengths, altro_ws .* scale, label="Altro", title="Random Linear Problem (warmstart)", xlabel="Horizon Length", ylabel="Time (milliseconds)")
scatter!(horizon_lengths, osqp_ws .* scale, label="OSQP")
png("benchmarks/random_linear_mpc/warmstart_horizon")