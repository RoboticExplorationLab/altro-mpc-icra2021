import Pkg; Pkg.activate(string(@__DIR__,"/..")); Pkg.instantiate();
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays, LinearAlgebra
using RobotDynamics
using Plots
using BenchmarkTools

include("src/grasp_model.jl")
include("src/visualize.jl")

N = 15
tf = 1.4
dt = tf/(N-1)

n = 6
m = 6
g = @SVector [0, 0, -9.81] # gravity
mu = .8                # friction constant
mass = .2               # mass
j = 1.                  # inertia
f = 3.                  # max grasp force

# rotational trajectory
θdd = .01*ones(N)
θd = [.01*t for t = 0:N-1]
θ = [.005*t^2 for t = 0:N-1]

# generate p v B matrices
p1_0 = [.0,1, 0]; v1_0 = [.0,-1, 0]
p2_0 = [.0,-1, 0]; v2_0 = [.0,1, 0]
p1, v1, B1 = generate_pvB_3D(p1_0, v1_0, θ)
p2, v2, B2 = generate_pvB_3D(p2_0, v2_0, θ)

model = SquareObject(n, m, mu, mass, j, f, g, [p1, p2], [v1, v2], [B1, B2])

x0 = @SVector [0.,1.,1.,0.,0.,0.]
xf = @SVector zeros(n)

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
goal = GoalConstraint(xf)
add_constraint!(conSet, goal, N)

include("src/new_constraints.jl")
for i = 1:N-1
    global model, u_ind, F1_ind, F2_ind, n, m
    local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2
    # Torque Balance
    B = [model.B[1][i] model.B[2][i]]
    t_bal = LinearConstraint(n, m, B, [θdd[i],0,0], Equality(), u_ind)
    add_constraint!(conSet, t_bal, i:i)

    # Max Grasp Force
    A = zeros(2, m)
    A[1,1:Int(m/2)] = model.v[1][i]
    A[2,1+Int(m/2):end] = model.v[2][i]
    max_f = LinearConstraint(n, m, A, model.f*ones(2), Inequality(), u_ind)
    add_constraint!(conSet, max_f, i:i)

    # SOCP friction Cone
    v1_i = model.v[1][i]
    A1 = (I - v1_i*v1_i')
    c1 = model.mu*v1_i
    nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
    add_constraint!(conSet, nc1, i:i)

    v2_i = model.v[2][i]
    A2 = (I - v2_i*v2_i')
    c2 = model.mu*v2_i
    nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
    add_constraint!(conSet, nc2, i:i)
end

# Problem
prob = Problem(model, obj, xf, tf, x0=x0, constraints=conSet);

u0 = @SVector [0, 1, model.mass*9.81/2, 0, -1, model.mass*9.81/2]
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

using Altro
opts = SolverOptions(
    verbose = 0,
    projected_newton_tolerance=4e-4,
    cost_tolerance_intermediate=4e-3,
    penalty_scaling=10.,
    penalty_initial=1.0,
    projected_newton=true
)

# normal solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
solve!(altro)
as = altro.stats

# benchmark solve
altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=false)
benchmark_solve!(altro)

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
    p = [model.p[1][t][2:3], model.p[2][t][2:3]]
    F = [F1[t], F2[t]]
    visualize_square([y[t],z[t]], θ[t], p, F, model.mass*model.g[2:3], r=1)
    frame(anim)
end
gif(anim, string(@__DIR__,"/altro_3D_w_rotation.gif"), fps=2)

plot([x y z xd yd zd])
png(string(@__DIR__,"/altro_3D_w_rotation.png"))

"""
Solve Statistics
  Total Iterations: 42
  Solve Time: 208.64059899999998 (ms)
"""

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
#     A[1,1:Int(m/2)] = model.v[1][i]
#     A[2,1+Int(m/2):end] = model.v[2][i]
#     @show A*U[i]
# end
