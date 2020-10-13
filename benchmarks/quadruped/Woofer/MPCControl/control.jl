function control!(
    torques::AbstractVector,
    x_est::AbstractVector,
    joint_pos::AbstractVector,
    joint_vel::AbstractVector,
    t::AbstractFloat,
    param::ControllerParams,
) 
    # get current leg positions
    param.cur_foot_loc = ForwardKinematicsAll(joint_pos)

    # prev phase -> cur_phase check contacts to regenerate swing
    param.cur_phase = get_phase(t, param)
    param.cur_phase_time = get_phase_time(t, param.cur_phase, param)
    
    param.active_feet = param.gait.contact_phases[param.cur_phase]
    coordinate_expander!(param.active_feet_12, param.active_feet)

    rot = MRP(x_est[4], x_est[5], x_est[6])
    regen_footstep = false

    # swing leg
    for i = 1:4
        if (t - param.last_replan_t) > param.replan_update
            regen_footstep = true
            param.last_replan_t = t
        end


        # calculate footstep and generate trajectory (stored in swing params) if needed
        if param.gait.contact_phases[param.prev_phase][i] == 1
            if param.gait.contact_phases[param.cur_phase][i] == 0
                param.next_foot_loc[i] =
                    footstep_location(x_est, rot, param.cur_phase, i, param)

                # make sure MPC accounts for this next foot location
                param.planner_foot_loc[i] = param.next_foot_loc[i]

                J = LegJacobian(joint_pos[SLegIndexToRange(i)], i)
                cur_foot_vel_i = J * joint_vel[SLegIndexToRange(i)]

                foot_trajectory(
                    x_est,
                    rot,
                    cur_foot_vel_i,
                    t,
                    t + param.gait.phase_times[param.cur_phase],
                    i,
                    param,
                    regen_z = true,
                )
                param.last_replan_t = t
            end
        end

        # actually calculate swing torques
        if param.gait.contact_phases[param.cur_phase][i] == 0
            # calculate current foot tip velocity
            J = LegJacobian(joint_pos[SLegIndexToRange(i)], i)
            cur_foot_vel_i = J * joint_vel[SLegIndexToRange(i)]

            if regen_footstep
                param.next_foot_loc[i] =
                    footstep_location(x_est, rot, param.cur_phase, i, param)

                # make sure MPC accounts for this next foot location
                param.planner_foot_loc[i] = param.next_foot_loc[i]

                foot_trajectory(
                    x_est,
                    rot,
                    cur_foot_vel_i,
                    t,
                    (t - param.cur_phase_time) +
                    param.gait.phase_times[param.cur_phase],
                    i,
                    param,
                    regen_z = false,
                )
            end

            swing_torque_i = swing_torques(
                x_est, 
                rot,
                cur_foot_vel_i,
                joint_pos[SLegIndexToRange(i)],
                t,
                i,
                param,
            )
            param.swing_torques[LegIndexToRange(i)] .= swing_torque_i
        end
    end
    param.prev_phase = param.cur_phase

    if (t - param.last_t) >= param.mpc_update
        # update MPC forces
        reference_trajectory!(x_est, param)
        foot_history!(t, param)
        # @time foot_forces!(x_est, param)
        foot_forces!(x_est, t, param)

        param.last_t = t
    end

    # needs to be negative so force is exerted by body on world
    param.mpc_torques = Force2Torque(-param.forces, joint_pos)

    torques .=
        param.active_feet_12 .* param.mpc_torques +
        (ones(12) - param.active_feet_12) .* param.swing_torques
end
