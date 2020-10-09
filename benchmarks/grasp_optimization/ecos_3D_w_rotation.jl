using LinearAlgebra
using Convex
using ECOS
using Plots
using RobotDynamics
using StaticArrays
using BenchmarkTools

include("src/visualize.jl")
include("src/grasp_model.jl")

ex = 2
include("example$ex.jl")

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

# solve
# bm = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0))
Convex.solve!(prob, ECOS.Optimizer)

# visualize
anim = Animation()
for t = 1:N-1
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1.value[2:3, t], F2.value[2:3, t]]
    plot([])
    visualize_square(Z.value[2:3,t], θ[t], p, F, o.mass*o.g[2:3], r=1, xlims=(-3,7), ylims=(-4,4.5))
    frame(anim)
end
gif(anim, string(@__DIR__, "/ecos_3D_w_rotation$ex.gif"), fps=2)

plot(Z.value')
png(string(@__DIR__, "/ecos_3D_w_rotation$ex.png"))

plot([])
for t = 1:3:N-1
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1.value[2:3, t], F2.value[2:3, t]]
    visualize_square(Z.value[2:3,t], θ[t], p, F, o.mass*o.g[2:3], r=1, xlims=(-3,7), ylims=(-4,4.5), fa=t/(N-1))
end
plot!([])
