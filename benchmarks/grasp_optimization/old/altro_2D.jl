import Pkg; Pkg.activate(string(@__DIR__,"/..")); Pkg.instantiate();
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays, LinearAlgebra
using RobotDynamics
using Plots
using BenchmarkTools

include("src/grasp_model.jl")
include("src/visualize.jl")

n = 4
m = 4
g = @SVector [0, -9.81] # gravity
mu = .5                 # friction constant
mass = .2               # mass
j = 1.                  # inertia
f = 3.                  # max grasp force
model = SquareObject(n, m, mu, mass, j, f, g, [0], [0], [0])

N = 15
tf = 1.4
dt = tf/(N-1)

x0 = @SVector [1.,1.,0.,0.]
xf = @SVector zeros(4)

# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n
F1_ind = 1+Int(n/2):n
F2_ind = 1+Int(n/2):n

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
tb = LinearConstraint(n, m, [0  1  0  -1.0], [0.], Equality(), 5:8)
add_constraint!(conSet, tb, 1:N-1)

# # Max Grasp Force
# A = [-1. 0 0 0;
#       0 0 1 0]
# f_bnd = 3.
# b = [f_bnd, f_bnd]
#
# tb = LinearConstraint(n, m, A, b, Inequality(), 5:8)
# add_constraint!(conSet, tb, 1:N-1)

# Friction Cone
include("src/new_constraints.jl")
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

u0 = @SVector [0.01, model.mass*9.81/2, -.01, model.mass*9.81/2]
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
    visualize_square([x[t],y[t]], 0., p, F, model.mass*model.g, r=1)
    frame(anim)
end
gif(anim, string(@__DIR__,"/altro_2D.gif"), fps=2)

"""
14.608 ms (133415 allocations: 11.23 MiB)

julia> U
14-element Array{SArray{Tuple{4},Float64,1,4},1}:
 [-0.0680814850750493, 0.28886181777182035, -0.06808148507502121, 0.2888618177718709]
 [-0.06616669330537372, 0.2902878765636018, -0.06616669330541214, 0.2902878765647139]
 [-0.06425190153596949, 0.2917139353558511, -0.06425190153600502, 0.2917139353565358]
 [-0.06233710976650353, 0.29313999414682645, -0.06233710976649909, 0.29313999414718384]
 [-0.06042231799698161, 0.2945660529400843, -0.06042231799696296, 0.29456605293950056]
 [-0.05850752622745947, 0.29599211173022244, -0.05850752622745592, 0.29599211172956574]
 [-0.05659273445800128, 0.2974181705240364, -0.05659273445796842, 0.2974181705229202]
 [-0.054677942688590164, 0.29884422931423194, -0.05467794268854398, 0.29884422931454413]
 [-0.05276315091910577, 0.3002702881057562, -0.05276315091909023, 0.300270288106951]
 [-0.05084835914958452, 0.3016963468989522, -0.050848359149579636, 0.30169634689878944]
 [-0.04893356738011079, 0.3031224056921389, -0.048933567380135656, 0.3031224056913555]
 [-0.0470187756106788, 0.30454846448147777, -0.047018775610697006, 0.3045484644817622]
 [-0.045103983841172646, 0.3059745232740514, -0.04510398384113046, 0.30597452327425034]
 [-0.04318919207166649, 0.3074005820651049, 2.391362161585882, 0.30740058206631793]
"""
