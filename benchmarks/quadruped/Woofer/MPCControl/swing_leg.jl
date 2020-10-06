function foot_trajectory(
    x_est::AbstractVector{T},
    rot::Rotation,
    foot_vel_b::AbstractVector{T},
    t0::T,
    tf::T,
    i::Int,
    param::ControllerParams;
    regen_z = true,
) where {T<:Number}
    #=
    	Generate a body relative trajectory for the ith foot via spline interpolation
    	This function is called when the phase is switched to a no contact phase

    	Updated foot location taken from param.next_foot_loc
    =#

    foot_loc_cur_b = param.cur_foot_loc[i]
    foot_loc_cur_n = x_est[SUnitRange(1, 3)] + rot*foot_loc_cur_b

    foot_vel_n = x_est[SUnitRange(7, 9)] + rot*foot_vel_b #+ rot*Rotations.skew(x_est[SUnitRange(10, 12)])*foot_loc_cur_b

    # generate cubic spline in x,y to get body relative foot trajectory
    A = @SMatrix [
        t0^3 t0^2 t0 1.0;
        tf^3 tf^2 tf 1.0;
        3 * t0^2 2 * t0 1.0 0.0;
        3 * tf^2 2 * tf 1.0 0.0
    ]

    b_x = @SVector [foot_loc_cur_n[1], param.next_foot_loc[i][1], foot_vel_n[1], 0.0]
    b_y = @SVector [foot_loc_cur_n[2], param.next_foot_loc[i][2], foot_vel_n[2], 0.0]

    if regen_z
        # generate cubic spline in z to enforce height constraint and terminal velocity constraint
        A_z = @SMatrix [
            t0^3 t0^2 t0 1.0;
            tf^3 tf^2 tf 1.0;
            (0.5 * (tf + t0))^3 (0.5 * (tf + t0))^2 (0.5 * (tf + t0)) 1;
            3 * tf^2 2 * tf 1.0 0.0
        ]

        b_z = @SVector [
            foot_loc_cur_n[3],
            woofer.geometry.foot_radius,
            param.swing.step_height,
            0.0,
        ]

        # FIXME: singular exception here?
        foot_trajectory = [A \ b_x; A \ b_y; A_z \ b_z]
    else
        # FIXME: singular exception here?
        foot_trajectory = [
            A \ b_x
            A \ b_y
            param.swing.foot_trajectories[i][SUnitRange(9,12)]
        ]
    end

    param.swing.foot_trajectories[i] = foot_trajectory
end

function swing_torques(
    x_est::AbstractVector{T},
    rot::Rotation,
    cur_vel_b::AbstractVector{T},
    α::AbstractVector{T},
    t::T,
    i::Integer,
    param::ControllerParams,
) where {T<:Number}
    #=
    	PD cartesian controller around swing leg trajectory

    	puts swing_torques into param
    =#

    t_p = @SVector [t^3, t^2, t, 1]
    t_v = @SVector [3 * t^2, 2 * t, 1, 0]

    r_des_n = @SVector [
        dot(param.swing.foot_trajectories[i][SUnitRange(1,4)], t_p),
        dot(param.swing.foot_trajectories[i][SUnitRange(5,8)], t_p),
        dot(param.swing.foot_trajectories[i][SUnitRange(9,12)], t_p),
    ]

    v_des_n = @SVector [
        dot(param.swing.foot_trajectories[i][SUnitRange(1,4)], t_v),
        dot(param.swing.foot_trajectories[i][SUnitRange(5,8)], t_v),
        dot(param.swing.foot_trajectories[i][SUnitRange(9,12)], t_v),
    ]

    r_des_b = rot \ (r_des_n - x_est[SUnitRange(1,3)])
    v_des_b = rot \ (v_des_n - x_est[SUnitRange(7,9)])

    J = LegJacobian(α, i)

    return J' * (
        param.swing.kp_cart * (r_des_b - param.cur_foot_loc[i]) +
        param.swing.kd_cart * (v_des_b - cur_vel_b)
    )
end
