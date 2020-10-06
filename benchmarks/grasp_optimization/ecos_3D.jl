using LinearAlgebra
using Convex
using ECOS
using Plots
using RobotDynamics
using StaticArrays
using BenchmarkTools

include("src/visualize.jl")
include("src/grasp_model.jl")

g = @SVector [0, 0, -9.81] # gravity
mu = .5                 # friction constant
mass = .2               # mass
j = 1.                  # inertia
f = 3.                  # max grasp force

# dimensions
n = 6
m = 6
T = 15 # horizon
dt = .1

# generate p v B matrices
p1_0 = [0, -1, 0]
v1_0 = [0, 1, 0]
p2_0 = [0, 1, 0]
v2_0 = [0, -1, 0]

p1 = [p1_0 for t = 1:T]
v1 = [v1_0 for t = 1:T]
p2 = [p2_0 for t = 1:T]
v2 = [v2_0 for t = 1:T]
B1 = [skew(p1[t]) for t = 1:T]
B2 = [skew(p2[t]) for t = 1:T]

o = SquareObject(n, m, mu, mass, j, f, g, [p1, p2], [v1, v2], [B1, B2])

# variables
F1 = Variable(Int(m/2), T-1)
F2 = Variable(Int(m/2), T-1)
Z = Variable(n, T)

# objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
objective = Q*sumsquares(Z[:,1:T-1]) + Qf*sumsquares(Z[:,T]) + R*sumsquares([F1;F2])
prob = minimize(objective)

# initial and final positions
prob.constraints += Z[:,1] == [0.,1,1,0,0,0]
prob.constraints += Z[:,T] == zeros(n)

# indices for convenience
pos_ind = 1:Int(n/2)
vel_ind = 1+Int(n/2):n

# constraints
for t = 1:T-1
    global o, F1, F2, f, prob, Z, dt

    # friction cone
    prob.constraints += norm((I - o.v[1][t]*o.v[1][t]')*F1[:, t]) <= o.mu*o.v[1][t]'*F1[:, t] # friction cone
    prob.constraints += norm((I - o.v[2][t]*o.v[2][t]')*F2[:, t]) <= o.mu*o.v[2][t]'*F2[:, t] # friction cone

    # max grasp force
    prob.constraints += o.v[1][t]'*F1[:, t] <= o.f
    prob.constraints += o.v[2][t]'*F2[:, t] <= o.f

    # force balance
    u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
    prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
    prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2

    # torque balance
    prob.constraints += zeros(3) == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]
end

# solve
solve!(prob, ECOS.Optimizer)

# visualize
anim = Animation()
for t = 1:T-1
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1.value[2:3, t], F2.value[2:3, t]]
    visualize_square(Z.value[2:3,t], 0., p, F, o.mass*o.g[2:3], r=1)
    frame(anim)
end
gif(anim, string(@__DIR__, "/ecos_3D.gif"), fps=2)

plot(Z.value')
png(string(@__DIR__, "/ecos_3D.png"))

"""
OPTIMAL (within feastol=8.2e-09, reltol=8.9e-10, abstol=1.4e-07).
Runtime: 0.089287 seconds.
"""
