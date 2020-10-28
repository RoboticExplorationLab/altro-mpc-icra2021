using LinearAlgebra
using StaticArrays
import StaticArrays: SUnitRange
using Rotations
using BenchmarkTools
using MuJoCo # MuJoCo.jl is in the Lyceum Registry

include("Woofer/QuadrupedDynamics.jl")
include("Woofer/Utilities.jl")
include("Woofer/Config.jl")

using .QuadrupedDynamics

include("Woofer/MPCControl/MPCControl.jl")
import .MPCControl

#################################################################################################################
# Closed Loop Simulation in MuJoCo:
# Trotting in Place
function mujoco_controller_test(controller, tf, mpc_dt, direct_time = false)
    m = jlModel("Woofer/woofer.xml")
    d = jlData(m)

    low_level_control_dt = 0.001
    last_control_update = -0.1
    last_mpc_update = -0.1

    steps = round(tf/m.opt.timestep)

    τ = zeros(12)

    num_samples = Integer(floor(tf/mpc_dt + 1))
    solve_times = zeros(num_samples)
    state_history = [zeros(12) for i=1:num_samples]
    time_history = zeros(num_samples)
    j = 1

    for i=1:steps
        x = get_state(d)
        q = get_joint_pos(d)
        q̇ = get_joint_vel(d)
        t = d.time

        if (t - last_control_update) >= low_level_control_dt
            # pull benchmark out of control function
            if (t-last_mpc_update) >= mpc_dt
                MPCControl.reference_trajectory!(x, controller)
                MPCControl.foot_history!(t, controller)
                b = MPCControl.foot_forces!(x, t, controller)

                solve_times[j] = direct_time ? b : b.times[1]
                state_history[j] = x
                time_history[j] = t

                j += 1

                last_mpc_update = t
            end

            τ = MPCControl.control!(τ, x, q, q̇, t, controller)

            d.ctrl .= τ

            last_control_update = t
        end

        mj_step(m, d);
    end

    return solve_times, state_history, time_history
end

function get_state(d)
        q = d.qpos
        q̇ = d.qvel
        rot = UnitQuaternion(q[4], q[5], q[6], q[7])
        mrp = MRP(rot)
        ω = rot \ q̇[SUnitRange(4, 6)]

        x = [   q[SUnitRange(1, 3)]; 
                Rotations.params(mrp); 
                q̇[SUnitRange(1, 3)]; 
                ω   ]

        return x
    end

# annoying way to get rid of knee joint measurements
get_joint_pos(d) = d.qpos[@SVector [8,9,11,13,14,16,18,19,21,23,24,26]]
get_joint_vel(d) = d.qvel[@SVector [7,8,10,12,13,15,17,18,20,22,23,25]]