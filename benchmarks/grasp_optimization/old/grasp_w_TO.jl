import Pkg; Pkg.activate(string(@__DIR__,"/..")); Pkg.instantiate();
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays, LinearAlgebra
using RobotDynamics
using Plots
include("grasp_model.jl")
include("visualize.jl")

g = [0, -9.81] # gravity
mu = .5     # friction constant
m = .2      # mass
j = 1.      # inertia

model = SquareObject(mu, m, j, g, [0], [0], [0])

n = 4
m = 4

N = 15
tf = 1.4
dt = tf/(N-1)

x0 = @SVector [1.,1.,0.,0.]
xf = @SVector zeros(4)

# Set up
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
tb = LinearConstraint(n, m, [0  1  0  -1.0], [0.], Equality(), 5:8)
add_constraint!(conSet, tb, 1:N-1)

# Max Grasp Force
A = [-1. 0 0 0;
      0 0 1 0]
f_bnd = 3.
b = [f_bnd, f_bnd]

tb = LinearConstraint(n, m, A, b, Inequality(), 5:8)
add_constraint!(conSet, tb, 1:N-1)

# Friction Cone
include("new_constraints.jl")
v1_0 = [-1., 0]
v2_0 = [1., 0]

A1 = (I - v1_0*v1_0')
c1 = model.mu*v1_0
nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), 5:6)
add_constraint!(conSet, nc1, 1:N-1)

A2 = (I - v2_0*v2_0')
c2 = model.mu*v2_0
nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), 7:8)
add_constraint!(conSet, nc2, 1:N-1)

prob = Problem(model, obj, xf, tf, x0=x0, constraints=conSet);

u0 = @SVector [0, model.m*9.81/2, 0, model.m*9.81/2]
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
solve!(altro);

# extract results
X = states(altro)
U = controls(altro)

x = [X[t][1] for t = 1:N]
y = [X[t][2] for t = 1:N]
xd = [X[t][3] for t = 1:N]
yd = [X[t][4] for t = 1:N]
F1 = [U[t][1:2] for t = 1:N-1]
F2 = [U[t][3:4] for t = 1:N-1]

# visualize
anim = Animation()
for t = 1:N-1
    global F1, F2, x, y, model
    p = [[1, 0], [-1, 0]]
    F = [F1[t], F2[t]]
    visualize_square([x[t],y[t]], 0., p, F, model.m*model.g, r=1)
    frame(anim)
end
gif(anim, string(@__DIR__,"/grasp_w_TO.gif"), fps=2)
