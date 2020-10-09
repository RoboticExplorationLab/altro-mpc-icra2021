using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StaticArrays, LinearAlgebra, Plots, BenchmarkTools
using TrajectoryOptimization, RobotDynamics
const TO = TrajectoryOptimization
using Convex, ECOS

include("src/grasp_model.jl")
include("src/visualize.jl")

ex = 2
include("example$ex.jl")

## ALTRO
# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n
u_ind = n .+ (1:m)
F1_ind = n .+ (1:Int(m/2))
F2_ind = n .+ (1+Int(m/2):m)

# objective
Q = 1.0e-3*Diagonal(@SVector ones(n))
Qf = 10.0*Diagonal(@SVector ones(n))
R = 1.0*Diagonal(@SVector ones(m))
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# Goal Constraint
goal = GoalConstraint(SVector{n}(xf))
TO.add_constraint!(conSet, goal, N)

include("src/new_constraints.jl")
for i = 1:N-1
    global model, u_ind, F1_ind, F2_ind, n, m
    local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2

    # Torque Balance
    B = [o.B[1][i] o.B[2][i]]
    t_bal = LinearConstraint(n, m, B, [θdd[i],0,0], Equality(), u_ind)
    TO.add_constraint!(conSet, t_bal, i:i)

    # Max Grasp Force
    A = zeros(2, m)
    A[1,1:Int(m/2)] = o.v[1][i]
    A[2,1+Int(m/2):end] = o.v[2][i]
    max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
    TO.add_constraint!(conSet, max_f, i:i)

    # SOCP friction cone
    v1_i = o.v[1][i]
    A1 = (I - v1_i*v1_i')
    c1 = o.mu*v1_i
    nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
    TO.add_constraint!(conSet, nc1, i:i)

    v2_i = o.v[2][i]
    A2 = (I - v2_i*v2_i')
    c2 = o.mu*v2_i
    nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
    TO.add_constraint!(conSet, nc2, i:i)
end

# Problem
prob = TO.Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet);

u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

using Altro
# for example 2, solves in 50 ms
opts = SolverOptions(
    verbose = 0,
    projected_newton_tolerance=1e-5,
    cost_tolerance_intermediate=1e-5,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)

# normal solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
Altro.solve!(altro)

# benchmark solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=false)
bm_altro = benchmark_solve!(altro, samples=10, evals=10)

## ALTRO W WARMSTART
X_altro = states(altro)
U_altro = controls(altro)

initial_controls!(altro, U_altro)
initial_states!(altro, X_altro)

x0_warm = SVector{n}(x0 + .1randn(n))
altro.solver_al.solver_uncon.x0 .= x0_warm
bm_altro_warm = benchmark_solve!(altro, samples=10, evals=10)

@assert states(altro)[1] ≈ x0_warm

## ECOS

# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n

# variables
F1 = Variable(Int(m/2), N-1)
F2 = Variable(Int(m/2), N-1)
Z = Variable(n, N)

# objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
objective = Q*sumsquares(Z[:,1:N-1]) + Qf*sumsquares(Z[:,N]) + R*sumsquares([F1;F2])
prob = minimize(objective)

# start and goal constraint
prob.constraints += Z[:,1] == x0
prob.constraints += Z[:,N] == xf

# stage constraints
for t = 1:N-1
    global o, F1, F2, prob, Z, dt, θdd

    # torque balance
    prob.constraints += [θdd[t],0,0] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]

    # max grasp force
    prob.constraints += o.v[1][t]'*F1[:, t] <= o.f
    prob.constraints += o.v[2][t]'*F2[:, t] <= o.f

    # friction cone
    prob.constraints += norm((I - o.v[1][t]*o.v[1][t]')*F1[:, t]) <= o.mu*o.v[1][t]'*F1[:, t] # friction cone
    prob.constraints += norm((I - o.v[2][t]*o.v[2][t]')*F2[:, t]) <= o.mu*o.v[2][t]'*F2[:, t] # friction cone

    # dynamics
    u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
    prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
    prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2

end

# normal solve
Convex.solve!(prob, ECOS.Optimizer)

# benchmark solve
bm_ecos = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0)) samples=10 evals=10

## ECOS WITH WARMSTART
prob.constraints[1].rhs = x0_warm
bm_ecos_warm = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0), warmstart=true) samples=10 evals=10

## COMPARE

# run times
plotly()
histogram(bm_altro.times/1e6,fillalpha=.5,label="ALTRO")
histogram!(bm_ecos.times/1e6,fillalpha=.5,label="ECOS")
histogram!(bm_altro_warm.times/1e6,fillalpha=.5,label="ALTRO w warm start")
histogram!(bm_ecos_warm.times/1e6,fillalpha=.5,label="ECOS w warm start")
#
# # controls
# Fz_altro = [U_altro[t][3] for t = 1:N-1]
# Fz_ecos = F1.value[3,:]
# plot([Fz_altro Fz_ecos], label=["ALTRO" "ECOS"])
#
# # trajectory
# plot([])
# for t = 1:3:N-1
#     local p, F
#     p = [o.p[1][t][2:3], o.p[2][t][2:3]]
#     F = [F1.value[2:3, t], F2.value[2:3, t]]
#     visualize_square(Z.value[2:3,t], θ[t], p, F, o.mass*o.g[2:3], r=1, xlims=(-3,7), ylims=(-4,4.5), fa=t/(N-1))
# end
# plot!([])
