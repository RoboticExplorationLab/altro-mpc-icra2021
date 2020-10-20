struct FrictionConstraint{T, S} <: TO.ControlConstraint
	m::S
	μ::T
	i::S
	function FrictionConstraint(m::S, μ::T, i::S) where {T <: Number, S <: Real}
		new{T, S}(m, μ, i)
	end
end

@inline TrajectoryOptimization.control_dim(con::FrictionConstraint) = con.m
@inline TrajectoryOptimization.sense(::FrictionConstraint) = Inequality()
Base.length(::FrictionConstraint) = 4

function TrajectoryOptimization.evaluate(con::FrictionConstraint, u::SVector)
	x_ind = 3*(con.i-1)+1
	y_ind = 3*(con.i-1)+2
	z_ind = 3*(con.i-1)+3

	return @SVector [
						u[x_ind] - con.μ*u[z_ind],
						-con.μ*u[z_ind] - u[x_ind],
						u[y_ind] - con.μ*u[z_ind],
						-con.μ*u[z_ind] - u[y_ind]
					]
end
