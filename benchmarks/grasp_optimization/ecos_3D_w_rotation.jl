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
# rotational trajectory
# θ0 = 0; θf= pi/4; θd0 = 0; θdf = .1
# t0 = 0; tf = (T-1)*dt
# c = compute_rot_traj_coeffs(t0, tf, [θ0; θf; θd0; θdf])
# θ = [dot(c, [t^3,t^2,t,1]) for t = 0:dt:tf]
# θdd = [dot(c, [6t,2,0,0]) for t = 0:dt:tf]
θdd = .01*ones(T)
θd = [.01*t for t = 0:T-1]
θ = [.005*t^2 for t = 0:T-1]

# generate p v B matrices
p1_0 = [.0,1, 0]; v1_0 = [.0,-1, 0]
p2_0 = [.0,-1, 0]; v2_0 = [.0,1, 0]
p1, v1, B1 = generate_pvB_3D(p1_0, v1_0, θ)
p2, v2, B2 = generate_pvB_3D(p2_0, v2_0, θ)

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
    global o, F1, F2, f, prob, Z, dt, θdd

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
    prob.constraints += [θdd[t],0,0] == o.B[1][t] * F1[:, t] + o.B[2][t] * F2[:, t]
end

# solve
Convex.solve!(prob, ECOS.Optimizer)

# visualize
anim = Animation()
for t = 1:T-1
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1.value[2:3, t], F2.value[2:3, t]]
    visualize_square(Z.value[2:3,t], θ[t], p, F, o.mass*o.g[2:3], r=1)
    frame(anim)
end
gif(anim, string(@__DIR__, "/ecos_3D_w_rotation.gif"), fps=2)

plot(Z.value')
png(string(@__DIR__, "/ecos_3D_w_rotation.png"))

"""
OPTIMAL (within feastol=8.2e-09, reltol=8.9e-10, abstol=1.4e-07).
Runtime: 0.089287 seconds.
"""
