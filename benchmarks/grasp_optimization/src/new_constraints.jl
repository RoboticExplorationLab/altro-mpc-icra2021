import RobotDynamics: state_dim, control_dim
import TrajectoryOptimization: StageConstraint, ConstraintSense, SecondOrderCone

"""
norm(Ay) <= c'y
(Ay, c'y) lies in soc
"""
struct NormConstraint2{S,D} <: StageConstraint
	n::Int
	m::Int
	A::SizedMatrix
	c::SVector
	sense::S
	inds::SVector{D,Int}
	function NormConstraint2(n::Int, m::Int, A::SizedMatrix, c::SVector, sense::ConstraintSense,
			inds=SVector{n+m}(1:n+m))
		if inds == :state
			inds = SVector{n}(1:n)
		elseif inds == :control
			inds = SVector{m}(n .+ (1:m))
		end
		@assert size(A,2) == length(c)
		@assert eltype(A) == eltype(c)
		new{typeof(sense),length(inds)}(n,m,A,c,sense,inds)
	end
end

function FrictionConstraint(n::Int, m::Int, A::AbstractMatrix, c::AbstractVector,
		sense::S, inds=1:n+m) where {S<:ConstraintSense}
	@assert size(A,1) == length(c)
	@assert eltype(A) == eltype(c)
	p,q = size(A)
	A = SizedMatrix{p,q}(A)
	c = SVector{q}(c)
	NormConstraint2(n, m, A, c, sense, inds)
end

@inline state_dim(con::NormConstraint2) = con.n
@inline control_dim(con::NormConstraint2) = con.m
@inline TO.sense(con::NormConstraint2) = con.sense
# @inline Base.length(::NormConstraint2) = 1
@inline Base.length(::NormConstraint2{SecondOrderCone,D}) where D = D + 1

# function TO.evaluate(con::NormConstraint2, z::AbstractKnotPoint)
# 	x = z.z[con.inds]
# 	return @SVector [x'x - con.val*con.val]
# end

function TO.evaluate(con::NormConstraint2{SecondOrderCone}, z::AbstractKnotPoint)
	v = con.A*z.z[con.inds]
	return push(v, con.c'*z.z[con.inds])
end

# function TO.jacobian!(∇c, con::NormConstraint2, z::AbstractKnotPoint)
# 	x = z.z[con.inds]
# 	∇c[1,con.inds] .= 2*x
# 	return false
# end

function TO.jacobian!(∇c, con::NormConstraint2{SecondOrderCone}, z::AbstractKnotPoint)
	∇c[1:length(con.inds),con.inds] = con.A
	∇c[1+length(con.inds),con.inds] = con.c
	return true
end

function TO.change_dimension(con::NormConstraint2, n::Int, m::Int, ix=1:n, iu=1:m)
	NormConstraint2(n, m, con.val, con.sense, ix[con.inds])
end
