import Pkg; Pkg.activate(string(@__DIR__,"/..")); Pkg.instantiate();
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays, LinearAlgebra
using RobotDynamics
using Plots
using BenchmarkTools

include("src/grasp_model.jl")
include("src/visualize.jl")

ex = 2
include("example$ex.jl")

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
add_constraint!(conSet, goal, N)

include("src/new_constraints.jl")
for i = 1:N-1
    global model, u_ind, F1_ind, F2_ind, n, m
    local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2

    # Torque Balance
    B = [o.B[1][i] o.B[2][i]]
    t_bal = LinearConstraint(n, m, B, [θdd[i],0,0], Equality(), u_ind)
    add_constraint!(conSet, t_bal, i:i)

    # Max Grasp Force
    A = zeros(2, m)
    A[1,1:Int(m/2)] = o.v[1][i]
    A[2,1+Int(m/2):end] = o.v[2][i]
    max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
    add_constraint!(conSet, max_f, i:i)

    # SOCP friction cone
    v1_i = o.v[1][i]
    A1 = (I - v1_i*v1_i')
    c1 = o.mu*v1_i
    nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
    add_constraint!(conSet, nc1, i:i)

    v2_i = o.v[2][i]
    A2 = (I - v2_i*v2_i')
    c2 = o.mu*v2_i
    nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
    add_constraint!(conSet, nc2, i:i)
end

# Problem
prob = Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet);

u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

using Altro
# for example 2, solves in 50 ms
opts = SolverOptions(
    verbose = 1,
    projected_newton_tolerance=1e-5,
    cost_tolerance_intermediate=1e-5,
    penalty_scaling=10.,
    penalty_initial=1.0,
    constraint_tolerance=1e-4
)

# normal solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
solve!(altro)

# benchmark solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=false)
bm_altro = benchmark_solve!(altro)

# extract results
X = states(altro)
U = controls(altro)

x = [X[t][1] for t = 1:N]
y = [X[t][2] for t = 1:N]
z = [X[t][3] for t = 1:N]
xd = [X[t][4] for t = 1:N]
yd = [X[t][5] for t = 1:N]
zd = [X[t][6] for t = 1:N]

F1 = [U[t][2:3] for t = 1:N-1]
F2 = [U[t][5:6] for t = 1:N-1]

# visualize
anim = Animation()
for t = 1:N-1
    global F1, F2, y, z, model, θ
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1[t], F2[t]]
    plot([])
    visualize_square([y[t],z[t]], θ[t], p, F, o.mass*o.g[2:3], r=1)
    frame(anim)
end
gif(anim, string(@__DIR__,"/altro_3D_w_rotation$ex.gif"), fps=2)

plot([x y z xd yd zd])
png(string(@__DIR__,"/altro_3D_w_rotation$ex.png"))

plot([])
for t = 1:2:N-1
    global F1, F2, y, z, model, θ
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1[t], F2[t]]
    visualize_square([y[t],z[t]], θ[t], p, F, o.mass*o.g[2:3], fa=t/(N-1))
end
plot!([])

# # verify in friction cone
# for i = 1:N-1
#     global X, U, dt, nc1
#     local t, z1, x1
#     t=1
#     z1 = KnotPoint(X[t], U[t], dt)
#     x1 = TO.evaluate(nc1, z1)
#     @show x1
# end
#
# for i = 1:N-1
#     global X, U, dt, nc2
#     local t, z1, x1
#     t=1
#     z1 = KnotPoint(X[t], U[t], dt)
#     x1 = TO.evaluate(nc2, z1)
#     @show x1
# end
#
#
# # verify max grasp force
# for i = 1:N-1
#     global U
#     local A
#     A = zeros(2, m)
#     A[1,1:Int(m/2)] = o.v[1][i]
#     A[2,1+Int(m/2):end] = o.v[2][i]
#     @show A*U[i]
# end
