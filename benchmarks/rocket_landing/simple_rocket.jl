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
using StatProfilerHTML
using JLD2
using Plots

include("../mpc.jl")

##
function gen_JuMP_rocket(prob0::Problem, opts::SolverOptions, optimizer;
        goal_constraint = false
    )
    prob = copy(prob0)
    TO.add_dynamics_constraints!(prob)
    n,m,N = size(prob)
    nlp = TrajOptNLP(prob, remove_bounds=true)
    NN = N*n + (N-1)*m
    dt = prob.Z[1].dt
    xinds = nlp.Z.Zdata.xinds
    uinds = nlp.Z.Zdata.uinds

    # Create the model 
    model = Model(optimizer)
    @variable(model, z[1:NN]) 
    z0 = vcat([Vector(RD.get_z(z)) for z in get_trajectory(prob0)]...)
    JuMP.set_start_value.(z, z0)    

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
    dts = fill(dt,N_mpc)
    dts[end] = 1
    c = [c.c for c in prob.obj]'dts

    @objective(model, Min, 0.5*dot(z,P,z) + dot(q,z) + c)

    # Dynamics Constraints
    A = Array(RD.get_A(prob.model.linmodel,1))
    B = Array(RD.get_B(prob.model.linmodel,1))
    d = Array(RD.get_d(prob.model.linmodel,1))
    for k = 1:N-1
        @constraint(model, A*z[xinds[k]] .+ B*z[uinds[k]] .+ d .== z[xinds[k+1]])
    end

    # Initial condition 
    @constraint(model, z[xinds[1]] .== prob.x0)

    # Goal constraint
    if goal_constraint 
        @constraint(model, z[xinds[N]] .== prob.xf)
    end

    # Ground constraint
    # for k = 2:N-1
    #     @constraint(model, z[xinds[1][3]] >= 0)
    # end

    # Thrust cone constraint
    con = prob_mpc.constraints[findfirst(x-> x isa NormConstraint, prob_mpc.constraints)]
    maxThrust = con.val
    for k = 1:N-1
        @constraint(model, [maxThrust, z[uinds[k]]...] in JuMP.SecondOrderCone()) 
    end

    # Thrust angle constraint
    cones = prob_mpc.constraints[findall(x-> x isa NormConstraint2, prob_mpc.constraints)]
    α_max = cones[1].c[3]
    for k = 1:N-1
        u1,u2,u3 = z[uinds[k]]
        @constraint(model, [α_max * u3, u1, u2] in JuMP.SecondOrderCone())
    end

    return model, z
end

function mpc_update(altro, prob_mpc, Z_track, t0, k_mpc)
    TO.set_initial_time!(prob_mpc, t0)

    # Propagate the system forward w/ noise
    x0 = discrete_dynamics(TO.integration(prob_mpc),
                                prob_mpc.model, prob_mpc.Z[1])
    noise = [(@SVector randn(3)) / 100; (@SVector randn(3)) / 1e6]
    x0 += noise
    TO.set_initial_state!(prob_mpc, x0)

    # Update tracking cost
    TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

    # Shift the initial trajectory
    RD.shift_fill!(prob_mpc.Z)

    # Shift the multipliers and penalties
    Altro.shift_fill!(TO.get_constraints(altro))
end

function dynamics_violation(prob,X,U)
    err = [zero(X[1]) for u in U]
    for k = 1:length(U)
        t = prob.Z[k].t
        dt = prob.Z[k].dt
        err = discrete_dynamics(TO.integration(prob), prob.model, X[k], U[k], t, dt) - X[k+1]
    end
    return maximum(norm.(err,Inf))
end

opts = SolverOptions(
    constraint_tolerance = 1e-4,
    cost_tolerance=1e-4,
    cost_tolerance_intermediate=1e-4,
    penalty_initial=1e-4,
    penalty_scaling=2.0
)

##
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
opts_mpc = SolverOptions(
    cost_tolerance=1e-6,
    constraint_tolerance=1e-6,
    projected_newton=false
)
optimizer = JuMP.optimizer_with_attributes(ECOS.Optimizer, 
    "verbose"=>false,
    "feastol"=>1e-12,
    "abstol"=>opts_mpc.cost_tolerance,
    "reltol"=>opts_mpc.cost_tolerance
)
optimizer = JuMP.optimizer_with_attributes(COSMO.Optimizer, "eps_abs"=>1e-12, "eps_rel"=>1e-12)
# optimizer = JuMP.optimizer_with_attributes(SCS.Optimizer, "eps"=>1e-2, "max_iters"=>10e3)

function test_mpc(prob_mpc, opts_mpc, optimizer)
    N_mpc = 21
    prob_mpc = gen_tracking_problem(prob, N_mpc)
    altro = ALTROSolver(prob_mpc, opts_mpc, show_summary=true)
    solve!(altro)

    t0 = 2prob.Z[2].t 
    k_mpc = 3
    mpc_update(altro, prob_mpc, Z_track, t0, k_mpc)
    model,z = gen_JuMP_rocket(prob_mpc, altro.opts, optimizer)
    solve!(altro)
    optimize!(model)

    # Compare results
    prob_ = copy(prob_mpc)
    TO.add_dynamics_constraints!(prob_)
    n,m,N = size(prob_)
    nlp = TrajOptNLP(prob_, remove_bounds=true)
    
    TO.hess_f!(nlp)
    P = nlp.data.G
    q = zeros(n+m, N)
    for k = 1:N
        q[1:n,k] .= prob_.obj[k].q
        q[n+1:end,k] .= prob_.obj[k].r
    end
    dt = prob_.Z[1].dt
    q[:,1:N-1] .*= dt
    q = q[1:end-m]
    dts = fill(prob_mpc.Z[1].dt,N_mpc)
    dts[end] = 1
    c = [c.c for c in prob_mpc.obj]'dts

    n,m,N = size(prob_mpc)
    inds = reshape(1:(n+m)*N,n+m,N)
    xinds = [z[1:n] for z in eachcol(inds)]
    uinds = [z[n+1:end] for z in eachcol(inds)][1:N-1]

    x_altro = vcat([Vector(RD.get_z(z)) for z in get_trajectory(altro)]...)
    norm(value.(z) - x_altro, Inf)
    X = [value.(z[ind]) for ind in xinds]
    U = [value.(z[ind]) for ind in uinds]
    err_X = maximum(norm.(X - states(altro), Inf))
    err_U = maximum(norm.(U - controls(altro), Inf))
    err_dyn = dynamics_violation(prob_mpc, X, U)

    dts[end] = 0
    Z = Traj(X, U, dts)
    x_ecos = value.(z)
    @show cost(prob_mpc)
    @show 0.5*dot(x_altro,P,x_altro) + q'x_altro + c
    @show 0.5*dot(x_ecos,P,x_ecos) + q'x_ecos + c
    @show cost(prob_mpc.obj, Z)
    @show objective_value(model) 
    return SA[err_X, err_U, err_dyn]
end

function get_opts(optimizer, opts)
    if optimizer === COSMO.Optimizer
        return ("eps_abs"=>opts.cost_tolerance, "eps_rel"=>opts.cost_tolerance)
    elseif optimizer == ECOS.Optimizer
        return ("feastol"=>opts.constraint_tolerance, "abstol"=>opts.cost_tolerance, "reltol"=>opts.cost_tolerance)
    elseif optimizer == SCS.Optimizer
        return ("eps"=>opt.cost_tolerance, "verbose"=>false)
    end
end

test_mpc(prob_mpc, opts_mpc, optimizer)
dts = fill(prob_mpc.Z[1].dt,N_mpc)
dts[end] = 1
c = [c.c for c in prob_mpc.obj]'dts



##
using DataFrames
res = Dict{Symbol,Vector{Union{Float64,Symbol}}}(
    :err_X => Float64[],
    :err_U => Float64[],
    :err_dyn => Float64[],
    :tol => Float64[],
    :solver => Symbol[]
)

for tol in (1e-4,1e-6,1e-8,1e-10)
    opts = SolverOptions(
        constraint_tolerance = tol,
        cost_tolerance=tol
    )
    for solver in (COSMO, ECOS)
        optimizer = JuMP.optimizer_with_attributes(solver.Optimizer, get_opts(solver.Optimizer, opts)...)
        err_X, err_U, err_dyn = test_mpc(opts_mpc, optimizer) 
        push!(res[:err_X], err_X)
        push!(res[:err_U], err_U)
        push!(res[:err_dyn], err_dyn)
        push!(res[:tol], tol)
        push!(res[:solver], Symbol(solver))
    end
end

df = DataFrame(res)
using Plots
cosmo_res = df[df.solver .== :COSMO,:]
ecos_res = df[df.solver .== :ECOS,:]
plot(title="COSMO", legend=:topleft, xlabel="tolerance", ylabel="error")
plot!(cosmo_res.tol, cosmo_res.err_X, label="err_X", xscale=:log10, yscale=:log10)
plot!(cosmo_res.tol, cosmo_res.err_U, label="err_U", xscale=:log10, yscale=:log10)
plot!(cosmo_res.tol, cosmo_res.err_dyn, label="err_dyn", xscale=:log10, yscale=:log10)

plot(title="ECOS", legend=:topleft, xlabel="tolerance", ylabel="error")
plot!(ecos_res.tol, ecos_res.err_X, label="err_X", xscale=:log10, yscale=:log10)
plot!(ecos_res.tol, ecos_res.err_U, label="err_U", xscale=:log10, yscale=:log10)
plot!(ecos_res.tol, ecos_res.err_dyn, label="err_dyn", xscale=:log10, yscale=:log10)