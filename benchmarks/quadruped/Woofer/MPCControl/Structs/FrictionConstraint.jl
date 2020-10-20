function FrictionConstraint(m::Integer, μ::T, i::Integer) where {T<:Real}
	inds = SLegIndexToRange(i)

	A = SizedMatrix{3,3, T}([1 0 0; 0 1 0; 0 0 0])
	c = μ*SVector{3, T}([0, 0, 1])

	return NormConstraint2(m, A, c, TO.SecondOrderCone(), inds)
end

struct NormConstraint2{S,D,p,q} <: TrajectoryOptimization.ControlConstraint 
    m::Int
    A::SizedMatrix{p, q, Float64, 2}
    c::SVector{p, Float64}
    sense::S
    inds::SVector{D,Int}
    function NormConstraint2(m::Int, A::SizedMatrix, c::SVector, sense::TrajectoryOptimization.ConstraintSense,
            inds = SVector{m}(n .+ (1:m)))
        @assert size(A,2) == length(c)
        @assert eltype(A) == eltype(c)

        new{typeof(sense),length(inds), size(A, 1), size(A, 2)}(m,A,c,sense,inds)
    end
end

@inline TrajectoryOptimization.control_dim(con::NormConstraint2) = con.m
@inline TrajectoryOptimization.sense(con::NormConstraint2) = con.sense
@inline Base.length(::NormConstraint2{TO.SecondOrderCone,D}) where D = D + 1

function TrajectoryOptimization.evaluate(con::NormConstraint2{TrajectoryOptimization.SecondOrderCone}, u::SVector)
    v = con.A*u[con.inds]
    return push(v, con.c'*u[con.inds])
end

function TrajectoryOptimization.jacobian!(∇c, con::NormConstraint2{TrajectoryOptimization.SecondOrderCone}, 
                                                    u::SVector)
    ∇c[1:length(con.inds),con.inds] = con.A
    ∇c[1+length(con.inds),con.inds] = con.c
    return true
end