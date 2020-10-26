import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

using Altro
using TrajectoryOptimization
using RobotDynamics
import RobotZoo.LinearModels: DoubleIntegrator
const RD = RobotDynamics
const TO = TrajectoryOptimization

using JuMP
using OSQP, ECOS
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

##
function RocketModel(mass, grav, dt, ωPlanet = [0.0; 0.0; 0.0];
                                                integration=RD.Exponential)
    A = [
        zeros(3,3)      I(3);
        zeros(3,3)      zeros(3,3)
    ]
    B = [
        zeros(3,3);
        -1/mass * I(3)
    ]
    d = [
        zeros(3);
        grav
    ]

    # Continuous Model
    cmodel =  LinearModel(A,B,d)
    if dt == zero(dt)
        return cmodel
    end

    # Discrete Model
    model = LinearizedModel(cmodel, dt=dt, is_affine=true, integration=integration)
end

function gen_rocket_problem(N, dt=0.1)
    model = RocketModel(10, SA[0,0,-9.81], dt) 
    n,m = size(model)
    tf = (N-1)*dt

    # Initial and Final conditions
    x0 = SA_F64[0,0,100, 5,5,-10]
    xf = @SVector zeros(n)

    # Objective
    Q = Diagonal(@SVector fill(1.0, n))
    R = Diagonal(@SVector fill(0.1, m))
    Qf = Diagonal(@SVector fill(100.0, n))
    obj = LQRObjective(Q,R,Qf,xf,N)

    # Constraints
    cons = ConstraintList(n,m,N)
    xmin = SA[-Inf,-Inf,0,-Inf,-Inf,-Inf]
    add_constraint!(cons, BoundConstraint(n, m, x_min=xmin), 2:N-1)
    add_constraint!(cons, GoalConstraint(xf), N)

    thrust = NormConstraint(n,m, 1000.0, TO.SecondOrderCone(), :control)
    add_constraint!(cons, thrust, 1:N-1)

    prob = Problem(model, obj, xf, tf, x0=x0, constraints=cons, integration=PassThrough)
    return prob
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
    q = zeros(n+m, N)
    for k = 1:N
        q[1:n,k] .= prob.obj[k].q
        q[n+1:end,k] .= prob.obj[k].r
    end
    dt = prob.Z[1].dt
    q[:,1:N-1] .*= dt
    q = q[1:end-m]

    # Constraints
    TO.jac_c!(nlp)
    A = nlp.data.D
    d = RD.get_d(prob.model.linmodel,1)
    gL, gU = TO.constraint_bounds(nlp)
    gL = -repeat(d,N-1)
    gU = copy(gL)
    zL, zU = TO.primal_bounds!(nlp)
    active = isfinite.(zL) .| isfinite.(zU)  # filter out the bounds at infinity

    # Put dynamics on top of bound constraints
    A = [A; I(NN)[active,:]]
    u = [gU; zU[active]]
    l = [gL; zL[active]]

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

function gen_JuMP(prob0::Problem, opts::SolverOptions)
    prob = copy(prob0)
    TO.add_dynamics_constraints!(prob)
    n,m,N = size(prob)
    nlp = TrajOptNLP(prob, remove_bounds=true)
    NN = N*n + (N-1)*m
    dt = prob.Z[1].dt
    xinds = nlp.Z.Zdata.xinds
    uinds = nlp.Z.Zdata.uinds

    # Put dynamics on top of bound constraints
    model = Model(JuMP.optimizer_with_attributes(ECOS.Optimizer, 
        "verbose"=>false,
        "feastol"=>opts.constraint_tolerance,
        "abstol"=>opts.cost_tolerance,
        "reltol"=>opts.cost_tolerance
    ))
    @variable(model, z[1:NN]) 

    # Cost function
    TO.hess_f!(nlp)
    P = nlp.data.G
    q = zeros(n+m, N)
    for k = 1:N
        q[1:n,k] .= prob.obj[k].q
        q[n+1:end,k] .= prob.obj[k].r
    end
    dt = prob.Z[1].dt
    q[:,1:N-1] .*= dt
    q = q[1:end-m]
    @objective(model, Min, dot(z,P,z) + dot(q,z))

    # Dynamics Constraints
    TO.jac_c!(nlp)
    A = Array(RD.get_A(prob.model.linmodel,1))
    B = Array(RD.get_B(prob.model.linmodel,1))
    d = Array(RD.get_d(prob.model.linmodel,1))
    b = -repeat(d,N-1)
    for k = 1:N-1
        @constraint(model, A*z[xinds[k]] .+ B*z[uinds[k]] .+ d .== z[xinds[k+1]])
    end

    zL, zU = TO.primal_bounds!(nlp)

    # @constraint(model, A * z .== b)
    @constraint(model, z[xinds[1]] .== prob.x0)
    @constraint(model, z[xinds[N]] .== prob.xf)
    for k = 2:N-1
        @constraint(model, z[xinds[1][3]] >= 0)
    end
    for k = 1:N-1
        @constraint(model, [1000.0, z[uinds[k]]...] in JuMP.SecondOrderCone()) 
    end
    return model, z
end

opts = SolverOptions(
    constraint_tolerance = 1e-8,
    cost_tolerance=1e-8,
    cost_tolerance_intermediate=1e-8
)

##
prob = gen_rocket_problem(21)
altro = ALTROSolver(prob, opts, show_summary=true, verbose=1, projected_newton=false)
rollout!(prob)
Altro.solve!(altro)
TO.findmax_violation(altro)
x_altro = vcat([RD.get_z(z) for z in get_trajectory(altro)]...)
norm.(controls(altro))

##
using JuMP
n,m,N = size(prob)
model, z = gen_JuMP(prob, altro.opts)
optimize!(model)
norm(value.(z) - x_altro, Inf)
inds = reshape(1:(n+m)*N,n+m,N)
xinds = [z[1:n] for z in eachcol(inds)]
uinds = [z[n+1:end] for z in eachcol(inds)][1:N-1]
X = [value.(z[ind]) for ind in xinds]
U = [value.(z[ind]) for ind in uinds]
maximum(norm.(X - states(altro), Inf))
maximum(norm.(U - controls(altro), Inf))
norm.(U)

# ## Solve using OSQP
# osqp, = gen_OSQP(prob, altro.opts) 
# res = OSQP.solve!(osqp)
# res.x
# norm(res.x - x_altro, Inf)

## Solve using Convex.jl