using LinearAlgebra
using StaticArrays
using Rotations

include("Woofer/QuadrupedDynamics.jl")
include("Woofer/MPCControl/MPCControl.jl")
include("Woofer/Utilities.jl")
include("Woofer/Config.jl")

using .QuadrupedDynamics
import .MPCControl

# IMPORTANT: to change which solver is used modify solver parameter in Woofer/MPCControl/MPC.yaml
param = MPCControl.ControllerParams(Float64, Int64)
N = param.N

x = [0.0, 0.0, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
MPCControl.reference_trajectory!(x, param)
t = 0.0
MPCControl.foot_history!(t, param)
(X, U) = MPCControl.foot_forces!(x, t, param)
# MPCControl.foot_forces!(x, t, param)

using YAML
data = YAML.load(open("Woofer/MPCControl/MPC.yaml"))

Q = Diagonal(data["q"])
R = Diagonal(data["r"])

select(i, n) = (n*(i-1)+1):(n*(i-1)+n)
ϵ = 1e-6
μ = data["mu"]
min_vert_force = data["min_vert_force"]
max_vert_force = data["max_vert_force"]

# check dynamics
for i=1:N-1
    @assert maximum(abs.(X[i+1] - (param.optimizer.A_vec[i]*X[i] + param.optimizer.B_vec[i]*U[i] + param.optimizer.d_vec[i]))) <= 1e-5

    for j=1:4
        u_j = U[i][select(j, 3)]

        @assert abs(u_j[1]) <= μ*u_j[3] + ϵ
        @assert abs(u_j[2]) <= μ*u_j[3] + ϵ

        @assert u_j[3] <= max_vert_force + ϵ
        @assert u_j[3] >= min_vert_force - ϵ
    end
end



objective_value_(X, U) = sum([0.5*(X[i] - param.x_ref[i])'*Q*(X[i] - param.x_ref[i]) for i=1:N]) + sum([0.5*U[i]'*R*U[i] for i=1:N-1])
@show objective_value_(X, U)

using Plots

plot(1:(N-1), [X[i][10] for i=1:(N-1)])

plot(1:(N-1), [U[i][1] for i=1:(N-1)])

# using MuJoCo # MuJoCo.jl is in the Lyceum Registry

# m = jlModel("woofer.xml")
# d = jlData(m)

# function get_state(d)
#     q = d.qpos
#     q̇ = d.qvel
#     rot = UnitQuaternion(q[4], q[5], q[6], q[7])
#     mrp = MRP(rot)
#     ω = rot \ q̇[SUnitRange(4, 6)]

#     x = [   q[SUnitRange(1, 3)]; 
#             Rotations.params(mrp); 
#             q̇[SUnitRange(1, 3)]; 
#             ω   ]

#     return x
# end

# # annoying way to get rid of knee joint measurements
# get_joint_pos(d) = d.qpos[@SVector [8,9,11,13,14,16,18,19,21,23,24,26]]
# get_joint_vel(d) = d.qvel[@SVector [7,8,10,12,13,15,17,18,20,22,23,25]]

# for i=1:100


#     mj_step(m, d);
#     println(d.qpos)
# end

