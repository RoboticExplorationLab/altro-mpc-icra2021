import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate();
using TrajectoryOptimization
using RobotDynamics
import RobotZoo.Cartpole
using StaticArrays, LinearAlgebra

using Plots
include("grasp_model.jl")
include("visualize.jl")

g = [0, -9.81] # gravity
mu = .5    # friction constant
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
Q = 1.0e-4*Diagonal(@SVector ones(n))
Qf = 100.0*Diagonal(@SVector ones(n))
R = 1.0e-1*Diagonal(@SVector ones(m))
obj = LQRObjective(Q,R,Qf,xf,N);

# Create Empty ConstraintList
conSet = ConstraintList(n,m,N)

# # Control Bounds
# u_bnd = 3.0
# bnd = BoundConstraint(n,m, u_min=-u_bnd, u_max=u_bnd)
# add_constraint!(conSet, bnd, 1:N-1)

# Goal Constraint
goal = GoalConstraint(xf)
add_constraint!(conSet, goal, N)

prob = Problem(model, obj, xf, tf, x0=x0, constraints=conSet);

u0 = @SVector fill(0.01,m)
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
gif(anim, "grasp_w_TO.gif", fps=2)
