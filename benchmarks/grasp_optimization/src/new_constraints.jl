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
	c::SizedVector
	sense::S
	inds::SVector{D,Int}
	function NormConstraint2(n::Int, m::Int, A::SizedMatrix, c::SizedVector, sense::ConstraintSense,
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
	c = SizedVector{q}(c)
	NormConstraint2(n, m, A, c, sense, inds)
end

@inline state_dim(con::NormConstraint2) = con.n
@inline control_dim(con::NormConstraint2) = con.m
@inline TO.sense(con::NormConstraint2) = con.sense
@inline Base.length(::NormConstraint2{TO.SecondOrderCone,D}) where D = D + 1

function TO.evaluate(con::NormConstraint2{TO.SecondOrderCone}, z::AbstractKnotPoint)
	v = con.A*z.z[con.inds]
	return push(v, con.c'*z.z[con.inds])
end

function TO.jacobian!(∇c, con::NormConstraint2{TO.SecondOrderCone}, z::AbstractKnotPoint)
	∇c[1:length(con.inds),con.inds] = con.A
	∇c[1+length(con.inds),con.inds] = con.c
	return true
end

function TO.change_dimension(con::NormConstraint2, n::Int, m::Int, ix=1:n, iu=1:m)
	NormConstraint2(n, m, con.val, con.sense, ix[con.inds])
end

"""
LinearConstraint with mutable b vector
"""
struct LinearConstraint2{S,P,W,T} <: StageConstraint
	n::Int
	m::Int
	A::SizedMatrix{P,W,T,2}
	b::SizedVector{P,T}
	sense::S
	inds::SVector{W,Int}
	function LinearConstraint2(n::Int, m::Int, A::StaticMatrix{P,W,T}, b::SizedVector{P,T},
			sense::ConstraintSense, inds=1:n+m) where {P,W,T}
		@assert length(inds) == W
		inds = SVector{W}(inds)
		new{typeof(sense),P,W,T}(n,m,A,b,sense,inds)
	end
end

function LinearConstraint2(n::Int, m::Int, A::AbstractMatrix, b::AbstractVector,
		sense::S, inds=1:n+m) where {S<:ConstraintSense}
	@assert size(A,1) == length(b)
	p,q = size(A)
	A = SizedMatrix{p,q}(A)
	b = SizedVector{p}(b)
	LinearConstraint2(n,m, A, b, sense, inds)
end

Base.copy(con::LinearConstraint2{S}) where S =
	LinearConstraint2(con.n, con.m, copy(con.A), copy(con.b), S(), con.inds)

@inline TO.sense(con::LinearConstraint2) = con.sense
@inline Base.length(con::LinearConstraint2{<:Any,P}) where P = P
@inline TO.state_dim(con::LinearConstraint2) = con.n
@inline TO.control_dim(con::LinearConstraint2) = con.m
TO.evaluate(con::LinearConstraint2, z::AbstractKnotPoint) = con.A*z.z[con.inds] .- con.b
function TO.jacobian!(∇c, con::LinearConstraint2, z::AbstractKnotPoint)
	∇c[:,con.inds] .= con.A
	return true
end

function TO.change_dimension(con::LinearConstraint2, n::Int, m::Int, ix=1:n, iu=1:m)
	inds0 = [ix; n .+ iu]  # indices of original z in new z
	inds = inds0[con.inds] # indices of elements in new z
	LinearConstraint2(n, m, con.A, con.b, con.sense, inds)
end
