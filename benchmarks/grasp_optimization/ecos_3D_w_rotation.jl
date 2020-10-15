using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using StaticArrays, LinearAlgebra, Plots, BenchmarkTools
using RobotDynamics
using Convex, ECOS

include("src/visualize.jl")
include("src/grasp_model.jl")

ex = 2
include("settings$ex.jl")

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

# static version of gif
plot([])
for t = 1:3:N-1
    local p, F
    p = [o.p[1][t][2:3], o.p[2][t][2:3]]
    F = [F1.value[2:3, t], F2.value[2:3, t]]
    visualize_square(Z.value[2:3,t], θ[t], p, F, o.mass*o.g[2:3], r=1, xlims=(-3,7), ylims=(-4,4.5), fa=t/(N-1))
end
plot!([])

# plot of state vector trajectory
plot(Z.value',
    xlabel="Time Step",
    ylabel="Coordinate",
    label = ["x" "y" "z" "xd" "yd" "zd"])

# plot of normal forces
normal1 = [dot(o.v[1][i], F1.value[:,i]) for i = 1:N-1]
normal2 = [dot(o.v[2][i], F2.value[:,i]) for i = 1:N-1]
plot([normal1 normal2 o.f*ones(N-1)],
    xlabel="Time Step",
    ylabel="Normal Force",
    linestyle = [:solid :solid :dash],
    label = ["F_N1" "F_N1" "Max Grasping Force"])

# plot of tangential to normal force ratio
tangent1 = [norm((I - o.v[1][i]*o.v[1][i]')*F1.value[:,i]) for i = 1:N-1]
tangent2 = [norm((I - o.v[2][i]*o.v[2][i]')*F2.value[:,i]) for i = 1:N-1]
friction1 = tangent1 ./ normal1
friction2 = tangent2 ./ normal2
plot([friction1 friction2 o.mu*ones(N-1)],
    xlabel="Time Step",
    ylabel="Tangential Force/Normal Force",
    linestyle = [:solid :solid :dash],
    label = ["F_T1/F_N1" "F_T2/F_N2" "mu"])
