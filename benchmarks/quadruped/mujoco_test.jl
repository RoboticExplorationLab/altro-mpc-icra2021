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
function mujoco_simulate(controller, tf, mpc_dt)
    m = jlModel(joinpath(@__DIR__, "Woofer/woofer.xml"))
    d = jlData(m)

    low_level_control_dt = 0.001
    last_control_update = -0.1
    last_mpc_update = -0.1

    steps = round(tf/m.opt.timestep)

    τ = zeros(12)

    num_samples = Integer(floor(tf/mpc_dt + 1))
    solve_times = zeros(num_samples)
    solve_iterations = zeros(Int64, num_samples)
    state_history = [zeros(12) for i=1:num_samples]
    force_history = [zeros(12) for i=1:num_samples]
    time_history = zeros(num_samples)
    j = 1

    for i=1:steps
        x = get_state(d)
        q = get_joint_pos(d)
        q̇ = get_joint_vel(d)
        t = d.time

        if (t - last_control_update) >= low_level_control_dt
            τ = MPCControl.control!(τ, x, q, q̇, t, controller)

            if controller.new_info 
                solve_times[j] = controller.last_solve_time
                state_history[j] = x
                time_history[j] = t
                force_history[j] = controller.forces
                solve_iterations[j] = controller.last_solve_iterations

                j += 1 

                controller.new_info = false
            end

            d.ctrl .= τ

            last_control_update = t
        end

        mj_step(m, d);
    end

    return solve_times, state_history, force_history, time_history, solve_iterations
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


### Test controllers simulatenously and ensure equal foot forces:
# c1 - altro controller
# c2 - general solver controller
function test_same_solution(c1, c2, mpc_dt; tf=0.5, direct_time = false, linearized_friction=true, test_tol=1e-5)
    m = jlModel(joinpath(@__DIR__, "Woofer/woofer.xml"))
    d = jlData(m)

    low_level_control_dt = 0.001
    last_control_update = -0.1
    last_mpc_update = -0.1

    steps = round(tf/m.opt.timestep)

    τ = zeros(12)

    num_samples = Integer(floor(tf/mpc_dt + 1))
    is_same = zeros(Bool, num_samples)
    max_diff = zeros(num_samples)
    j = 1

    for i=1:steps
        x = get_state(d)
        q = get_joint_pos(d)
        q̇ = get_joint_vel(d)
        t = d.time

        if (t - last_control_update) >= low_level_control_dt
            τ = MPCControl.control!(τ, x, q, q̇, t, c1)

            # calculate but don't actually use outpt from c2
            MPCControl.control!(τ, x, q, q̇, t, c2)
            if c1.new_info 
                @assert c2.new_info

                # max_x_deviation = 0.0
                # max_u_deviation = 0.0
                # for i=1:c1.N
                #     max_x_deviation = max(max_x_deviation, maximum(abs.(c1.optimizer.X0[i] - c2.optimizer.X0[i])))
                #     i == c1.N && continue
                #     max_u_deviation = max(max_u_deviation, maximum(abs.(c1.optimizer.U0[i] - c2.optimizer.U0[i])))
                # end
                # @show max_x_deviation
                # @show max_u_deviation

                # both solvers have same dynamics and meet friction constraints!
                for i=1:c1.N-1
                    @assert c1.optimizer.model.A[i] ≈ c2.optimizer.A_vec[i]
                    @assert c1.optimizer.model.B[i] ≈ c2.optimizer.B_vec[i]
                    @assert c1.optimizer.model.d[i] ≈ c2.optimizer.d_vec[i]

                    @assert c2.optimizer.X0[i+1] ≈ (c2.optimizer.A_vec[i]*c2.optimizer.X0[i] + c2.optimizer.B_vec[i]*c2.optimizer.U0[i] + c2.optimizer.d_vec[i]) atol=dynamics_tol
                    @assert c1.optimizer.X0[i+1] ≈ (c1.optimizer.model.A[i]*c1.optimizer.X0[i] + c1.optimizer.model.B[i]*c1.optimizer.U0[i] + c1.optimizer.model.d[i]) atol=dynamics_tol

                    for j=1:4
                        if linearized_friction
                            check_linearized_friction(c1.optimizer.U0[i], j)
                            check_linearized_friction(c2.optimizer.U0[i], j)
                        else
                            check_friction(c1.optimizer.U0[i], j)
                            check_friction(c2.optimizer.U0[i], j)
                        end
                    end
                end

                # check costs: 
                @assert c1.x_des ≈ c2.x_des

                Q = c1.optimizer.objective[1].Q
                R = c1.optimizer.objective[1].R
                cost_altro = get_cost(c1.optimizer.X0, c1.optimizer.U0, Q, R, c1.x_des)
                cost_other = get_cost(c2.optimizer.X0, c2.optimizer.U0, Q, R, c2.x_des)

                @show cost_altro
                @show cost_other

                is_same[j] = c1.forces ≈ c2.forces
                j += 1 

                c1.new_info = false
                c2.new_info = false
            end

            d.ctrl .= τ

            last_control_update = t
        end

        mj_step(m, d);
    end

    return count(isone, is_same)/j
end

function check_linearized_friction(u, j, μ=0.5, ϵ=1e-4)
    xind = 3*(j-1)+1
    yind = 3*(j-1)+2
    zind = 3*(j-1)+3

    @assert u[xind] <= μ*u[zind] + ϵ
    @assert u[xind] >= -μ*u[zind] - ϵ

    @assert u[yind] <= μ*u[zind] + ϵ
    @assert u[yind] >= -μ*u[zind] - ϵ

    @assert u[zind] >= -ϵ
end

function check_friction(u, j, μ=0.5, ϵ=1e-4)
    xind = 3*(j-1)+1
    yind = 3*(j-1)+2
    zind = 3*(j-1)+3

    @assert norm(u[xind:yind]) <= μ*u[zind] + ϵ
    @assert u[zind] >= -ϵ
end

function get_cost(X, U, Q, R, x̄)
    return sum([(X[i] - x̄)'*Q*(X[i] - x̄) for i=1:length(X)]) + sum([(U[i])'*R*(U[i]) for i=1:length(U)])
end

## SOLUTION TESTING CODE:
# altro_socp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = false, tol=1e-8)
# ecos_controller = MPCControl.ControllerParams(solver = "ECOS", linearized_friction = false, tol=1e-8)
# test_same_solution(altro_socp_controller, ecos_controller, altro_socp_controller.mpc_update, linearized_friction=false, direct_time = false, test_tol=1e-4)
# altro_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = true, tol=1e-6)
# osqp_controller = MPCControl.ControllerParams(solver = "OSQP", linearized_friction = true, tol=1e-10)
# test_same_solution(altro_controller, osqp_controller, osqp_controller.mpc_update, linearized_friction=true, direct_time = false)