using LinearAlgebra
using Convex
using ECOS
using Plots
using RobotDynamics
using StaticArrays
using BenchmarkTools

include("src/visualize.jl")
include("src/grasp_model.jl")

g = @SVector [0, -9.81] # gravity
mu = .5                 # friction constant
mass = .2               # mass
j = 1.                  # inertia
f = 3.                  # max grasp force

# dimensions
n = 4
m = 4
T = 15 # horizon
dt = .1

# generate p v B matrices
p1_0 = [1, 0]
v1_0 = [-1, 0]
p2_0 = [-1, 0]
v2_0 = [1, 0]

θdd = .01*ones(T)
θd = [.01*t for t = 0:T-1]
θ = [.005*t^2 for t = 0:T-1]

p1 = [rot(θ[t])*p1_0 for t = 1:T]
v1 = [rot(θ[t])*v1_0 for t = 1:T]
p2 = [rot(θ[t])*p2_0 for t = 1:T]
v2 = [rot(θ[t])*v2_0 for t = 1:T]
B1 = [[-p1[t][2] p1[t][1]] for t = 1:T]
B2 = [[-p2[t][2] p2[t][1]] for t = 1:T]

o = SquareObject(n, m, mu, mass, j, f, g, [p1, p2], [v1, v2], [B1, B2])

# variables
F1 = Variable(2, T-1)
F2 = Variable(2, T-1)
Z = Variable(4, T)

# objective
Q = 1.0e-3
Qf = 10.0
R = 1.0
objective = Q*sumsquares(Z[:,1:T-1]) + Qf*sumsquares(Z[:,T]) + R*sumsquares([F1;F2])
prob = minimize(objective)

# initial and final positions
prob.constraints += Z[:,1] == [1.;1.;0.;0.]
prob.constraints += Z[:,T] == zeros(4)

# constraints
for t = 1:T-1
    global o, F1, F2, f, prob, Z, dt, qdd

    # friction cone
    prob.constraints += norm((I - o.v[1][t]*o.v[1][t]')*F1[:, t]) <= o.mu*o.v[1][t]'*F1[:, t] # friction cone
    prob.constraints += norm((I - o.v[2][t]*o.v[2][t]')*F2[:, t]) <= o.mu*o.v[2][t]'*F2[:, t] # friction cone

    # max grasp force
    prob.constraints += o.v[1][t]'*F1[:, t] <= o.f
    prob.constraints += o.v[2][t]'*F2[:, t] <= o.f

    # force balance
    u = 1/o.m * (F1[:, t] + F2[:, t]) + o.g
    prob.constraints += Z[3:4, t+1] == Z[3:4, t] + u*dt
    prob.constraints += Z[1:2, t+1] == Z[1:2, t] + Z[3:4, t]*dt + u*.5*dt^2

    # torque balance
    prob.constraints += θdd[t] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]
end

# solve
solve!(prob, ECOS.Optimizer)

# visualize
anim = Animation()
for t = 1:T-1
    p = [o.p[1][t], o.p[2][t]]
    F = [F1.value[:, t], F2.value[:, t]]
    visualize_square(Z.value[1:2,t], θ[t], p, F, o.m*o.g, r=1)
    frame(anim)
end
gif(anim, string(@__DIR__, "/ecos_2D_w_rotation.gif"), fps=2)

plot(Z.value')
