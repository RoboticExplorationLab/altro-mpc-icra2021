function footstep_location(x_est::AbstractVector{T}, rot::Rotation, cur_phase::Integer, i::Integer, param::ControllerParams) where {T<:Number}
	# implement body velocity heuristic to get next body relative foot location
	ω_n = rot*x_est[SUnitRange(10,12)]
	v_n = x_est[SUnitRange(7,9)]
	p = x_est[SUnitRange(1,3)]

	next_phase = get_next_phase(cur_phase, param)
	t_next = param.gait.phase_times[next_phase]

	# k = sqrt(param.x_des[3]/9.81) # capture point heuristic
	k = T(0.0)

	v_des = param.x_des[SUnitRange(7,9)]

	nom_foot_loc_n = p + rot*param.nom_foot_loc[i]

	# projection matrix removes z component 
	projection_matrix = @SMatrix [	1 0 0;
									0 1 0;
									0 0 0	]

	next_foot_loc = nom_foot_loc_n +
						param.gait.alpha*t_next*v_n + k*(v_n - v_des)
						param.gait.beta*RotZ(t_next*ω_n[3])*param.cur_foot_loc[i]

	return projection_matrix * next_foot_loc + (@SVector [0., 0., woofer.geometry.foot_radius])
end

function foot_history!(t::Number, param::ControllerParams)
	# construct the contact and foot location history for MPC solver
	# inputs:
	# state reference trajectory

	# outputs (in param):
	# contacts
	# foot_locs
	# temporarily: integrated velocity -> position reference

	t_i = t + param.optimizer.dt

	prev_phase = get_phase(t, param)

	prev_foot_locs = zero(FootstepLocation)
	rot = MRP(param.x_ref[1][4], param.x_ref[1][5], param.x_ref[1][6])
	prev_foot_locs = param.x_ref[1][SUnitRange(1,3)] + rot*param.cur_foot_loc

	# current contact is first entry
	param.contacts[1] = param.gait.contact_phases[prev_phase]

	# cur_foot_loc is first entry
	param.foot_locs[1] .= prev_foot_locs

	for i in 2:(param.N+1)
		next_phase = get_phase(t_i, param)

		param.contacts[i] = param.gait.contact_phases[next_phase]

		x_ref_i = param.x_ref[i]
		v_i = x_ref_i[SUnitRange(7,9)]
		ω_b_i = x_ref_i[SUnitRange(10,12)]

		# rotation matrix from body to inertial
		rot = MRP(x_ref_i[4], x_ref_i[5], x_ref_i[6])

		for j in 1:4
			if param.gait.contact_phases[prev_phase][j] == 1
				if param.gait.contact_phases[next_phase][j] == 0
					# next foot placement must be planned prior to foot being released
					param.planner_foot_loc[j] = footstep_location(x_ref_i, rot, next_phase, j, param)
					# prev_foot_locs[j] = param.planner_foot_loc[j]
				end
			else
				if param.gait.contact_phases[next_phase][j] == 1
					prev_foot_locs[j] = param.planner_foot_loc[j]
				end
			end
		end

		param.foot_locs[i] .= prev_foot_locs

		t_i += param.optimizer.dt
		prev_phase = next_phase
	end
end
