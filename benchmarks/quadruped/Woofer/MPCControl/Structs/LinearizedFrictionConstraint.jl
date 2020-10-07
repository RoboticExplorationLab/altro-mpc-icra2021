struct FrictionConstraint{T} <: TO.ControlConstraint
	m::Int
	μ::T
	function FrictionConstraint(m::Int, μ::T) where T
		new{T}(m, μ)
	end
end

@inline TrajectoryOptimization.control_dim(con::FrictionConstraint) = con.m
@inline TrajectoryOptimization.sense(con::FrictionConstraint) = Inequality()
Base.length(::FrictionConstraint) = 16

function TrajectoryOptimization.evaluate(con::FrictionConstraint, u::SVector)
	return @SVector [
						u[1] - con.μ*u[3],
						-con.μ*u[3] - u[1],
						u[2] - con.μ*u[3],
						-con.μ*u[3] - u[2],
						u[4] - con.μ*u[6],
						-con.μ*u[6] - u[4],
						u[5] - con.μ*u[6],
						-con.μ*u[6] - u[5],
						u[7] - con.μ*u[9],
						-con.μ*u[9] - u[7],
						u[8] - con.μ*u[9],
						-con.μ*u[9] - u[8],
						u[10] - con.μ*u[12],
						-con.μ*u[12] - u[10],
						u[11] - con.μ*u[12],
						-con.μ*u[12] - u[11],
					]
end
