import Pkg; Pkg.activate(string(@__DIR__,"/..")); Pkg.instantiate();
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays, LinearAlgebra
using RobotDynamics
using Plots
using BenchmarkTools

include("src/grasp_model.jl")
include("src/visualize.jl")

n = 6
m = 6
g = @SVector [0, 0, -9.81] # gravity
mu = .8                # friction constant
mass = .2               # mass
j = 1.                  # inertia
f = 3.                  # max grasp force
model = SquareObject(n, m, mu, mass, j, f, g, [0], [0], [0])

N = 15
tf = 1.4
dt = tf/(N-1)

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

# Torque Balance
p1_0 = [0., -1, 0]
p2_0 = [0, 1, 0]
B = [skew(p1_0) skew(p2_0)]
t_bal = LinearConstraint(n, m, B, [0.,0,0], Equality(), u_ind)
add_constraint!(conSet, t_bal, 1:N-1)

# Max Grasp Force
v1_0 = [0., 1, 0]
v2_0 = [0., -1, 0]
A = zeros(2, m)
A[1,1:Int(m/2)] = v1_0
A[2,1+Int(m/2):end] = v2_0
b = model.f*ones(2)
max_f = LinearConstraint(n, m, A, b, Inequality(), u_ind)
add_constraint!(conSet, max_f, 1:N-1)

# Linear friction Cone
# A = [0. -.8 1 0 0 0; 0 0 0 0 .8 1]
# b = zeros(2)
# fc = LinearConstraint(n, m, A, b, Inequality(), u_ind)
# add_constraint!(conSet, fc, 1:N-1)

# SOCP friction Cone
include("src/new_constraints.jl")

A1 = (I - v1_0*v1_0')
c1 = model.mu*v1_0
nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
add_constraint!(conSet, nc1, 1:N-1)

A2 = (I - v2_0*v2_0')
c2 = model.mu*v2_0
nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
add_constraint!(conSet, nc2, 1:N-1)

# Problem
prob = Problem(model, obj, xf, tf, x0=x0, constraints=conSet);

u0 = @SVector [0, 0, model.mass*9.81/2, 0, 0, model.mass*9.81/2]
U0 = [u0 for k = 1:N-1]
initial_controls!(prob, U0)
rollout!(prob);

using Altro
opts = SolverOptions(
    cost_tolerance_intermediate=1e-2,
    penalty_scaling=10.,
    penalty_initial=1.0
)

altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
solve!(altro)

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
    global F1, F2, y, z, model
    local p, F
    p = [[-1, 0], [1, 0]]
    F = [F1[t], F2[t]]
    visualize_square([y[t],z[t]], 0., p, F, model.mass*model.g[2:3], r=1)
    frame(anim)
end
gif(anim, string(@__DIR__,"/altro_3D.gif"), fps=2)

plot([x y z xd yd zd])
png(string(@__DIR__,"/altro_3D.png"))

"""
Solve Statistics
  Total Iterations: 21
  Solve Time: 461.56739999999996 (ms)
"""

# verify in friction cone
# t=1
# z = KnotPoint(X[t], U[t], dt)
# x = TO.evaluate(nc1, z)
# @show in(x, TO.SecondOrderCone())
