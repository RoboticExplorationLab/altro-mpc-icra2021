function get_phase(t::AbstractFloat, param::ControllerParams)
	phase_time = t % param.gait.phase_length

	for i in 1:param.gait.num_phases
		if phase_time < sum(param.gait.phase_times[get_phase_index(i)])
			return i
		end
	end
end

function get_phase_time(t::AbstractFloat, phase::Integer, param::ControllerParams)
	phase_time = t % param.gait.phase_length

	return phase_time - sum(param.gait.phase_times[get_phase_index(phase-1)])
end

function get_next_phase(phase::Integer, param::ControllerParams)
	if (phase == param.gait.num_phases)
		return 1
	else
		return phase+1
	end
end

function get_phase_index(i::Integer)
	return SVector{i}(1:i)
end

# function nextFootTimeOnGround(i, phase::Integer, gait_params::GaitParams)
# """
# Finds the next length of time for which foot i will be on the ground
# """
# 	# TODO: makes planner more general
# end

function coordinate_expander!(expanded::Vector, compact::AbstractVector)
	expanded[1:3] .= compact[1]
	expanded[4:6] .= compact[2]
	expanded[7:9] .= compact[3]
	expanded[10:12] .= compact[4]
end
