struct SwingLegParams{T}
	foot_trajectories::Vector{SVector{12, T}}

	# z distance below robot of foot at middle of swing leg trajectory
	step_height::T

	kp_cart::T
	kd_cart::T

	function SwingLegParams(step_height::T, wn_cart::T, zeta_cart::T) where {T<:Number}
		foot_trajectories = [@SVector zeros(T, 12) for i=1:4]

		kp_cart = wn_cart^2
		kd_cart = 2*wn_cart*zeta_cart

		new{T}(foot_trajectories, step_height, kp_cart, kd_cart)
	end
end
