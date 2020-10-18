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

include(joinpath(dirname(@__FILE__),"generate_AB.jl"))

##
function gen_random_linear(n,m,N)
    # Create model
    dt = 0.5 # doesn't matter, just needs to be non-zero
    # A,B = gendiscrete(n,m)
    A,B = generate_AB()
    n,m = size(B)
    model = RD.LinearModel(A, B; dt=dt)

    Q = 10*I(n)
    R = .1*I(m)

    Qf = copy(Q)

    # x̄ = rand(n) .+ 1
    ū = .01*ones(m)
    # x0 = (rand(n) .- 1) .* x̄ * 0.5
    xf = zeros(n)
    x0 = [.1;.1;.1;zeros(n-3)]

    # Objective
    # Q = Diagonal(10*rand(n))
    # R = Diagonal(0.1*ones(m))
    # Qf = Q * (N-1)

    # x̄ = rand(n) .+ 1
    # ū = rand(m) .+ 0.5
    # x0 = (rand(n) .- 1) .* x̄ * 0.5
    # xf = zeros(n)
    obj = LQRObjective(Q, R, Qf, xf, N)

    # Constraints
    constraints = ConstraintList(n, m, N)
    bound = BoundConstraint(n, m, u_min=-ū, u_max=ū)
    add_constraint!(constraints, bound, 1:N-1)
    # add_constraint!(constraints, GoalConstraint(xf), N)

    # Problem
    tf = (N-1)*dt
    prob = Problem(model, obj, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)

end

function gen_OSQP_JuMP(prob::Problem)
    n,m,N = size(prob)
    A,B = Matrix(RD.get_A(prob.model)), Matrix(RD.get_B(prob.model))
    x0 = prob.x0
    xf = prob.xf
    @infiltrate
    error()
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
        # @constraint(jump_model, x[select(i, n)] .<= x̄)
        # @constraint(jump_model, x[select(i, n)] .>= -x̄)
        @constraint(jump_model, u[select(i, m)] .<= ū)
        @constraint(jump_model, u[select(i, m)] .>= -ū)
    end

    @objective(jump_model, Min, objective_exp)
    return jump_model
end

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

function gen_ICs(prob, iters=10)
    n,m,N = size(prob)
    x̄ = prob.constraints[1].z_max[1:n]
    [(rand(n) .- 1) .* x̄ * 0.5 for i = 1:iters]
end

opts = SolverOptions(
    cost_tolerance = 1e-4,
    constraint_tolerance = 1e-4
)

## Solve with ALTRO
prob = gen_random_linear(12,6,11)
solver = ALTROSolver(prob, show_summary=true)
solve!(solver)
cost(solver)

# Solve with OSQP / JuMP
jump_model = gen_OSQP_JuMP(prob)
optimize!(jump_model)
termination_status(jump_model)

# Solve with OSQP.jl
model,l,u = gen_OSQP(prob, opts)
results = OSQP.solve!(model)
results.prim_inf_cert
results.info.status
n,m,N = size(prob)
NN = N*n + (N-1)*m
xi = vcat([(k-1)*(n+m) .+ (1:n) for k = 1:N]...)
ui = vcat([(k-1)*(n+m) + n .+ (1:m) for k = 1:N-1]...)

## Check solutions
abs(objective_value(jump_model) - cost(solver))
abs(results.info.obj_val - cost(solver))

X_altro = vcat(Vector.(states(solver))...)
X_jump = value.(jump_model.obj_dict[:x])
X_osqp = results.x[xi]
norm(X_altro - X_jump)
norm(X_altro - X_osqp)

U_altro = vcat(Vector.(controls(solver))...)
U_jump = value.(jump_model.obj_dict[:u])
U_osqp = results.x[ui]
norm(U_altro - U_jump)
norm(U_altro - U_osqp)


## Update initial condition
ICs = gen_ICs(prob)
TO.set_initial_state!(solver, ICs[1])
solve!(solver)

x0_l = view(l, (N-1)*n .+ (1:n))
x0_u = view(u, (N-1)*n .+ (1:n))
x0_l .= ICs[1]
x0_u .= ICs[1]
OSQP.update!(model, u=u, l=l)
results = OSQP.solve!(model)

results.info.obj_val
norm(states(solver)[1] - ICs[1])
norm(results.x[1:n] - ICs[1])

function MPC_Altro(solver, ICs)
    times = Float64[]
    for ic in ICs
        t = @elapsed begin
            TO.set_initial_state!(solver, ic)
            solve!(solver)
        end
        push!(times, t)
    end
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
ICs = gen_ICs(prob)
set_options!(solver, show_summary = false)
t_altro = MPC_Altro(solver, ICs)
t_osqp = MPC_OSQP(model, l, u, x0_l, x0_u, ICs)
