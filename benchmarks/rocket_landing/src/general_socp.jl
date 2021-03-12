
using Altro, TrajectoryOptimization

"""
    NormConstraint2{S,D,p,q}

Generalized Norm Constraints for cases of the form
||Ax|| ≤ c'x
"""
struct NormConstraint2{S,D,p,q} <: TrajectoryOptimization.StageConstraint
    n::Int
    m::Int
    A::SizedMatrix{p, q, Float64, 2}
    c::SVector{p, Float64}
    sense::S
    inds::SVector{D,Int}
    function NormConstraint2(n::Int, m::Int, A::StaticMatrix, c::StaticVector,
                                sense::TrajectoryOptimization.ConstraintSense,
                                inds = SVector{n+m}(1:n+m))
        if inds == :state
			inds = SVector{n}(1:n)
		elseif inds == :control
			inds = SVector{m}(n .+ (1:m))
		end
        @assert size(A,2) == length(c)
        @assert eltype(A) == eltype(c)

        new{typeof(sense),length(inds), size(A, 1), size(A, 2)}(n,m,A,c,sense,inds)
    end
end

@inline TrajectoryOptimization.state_dim(con::NormConstraint2) = con.n
@inline TrajectoryOptimization.control_dim(con::NormConstraint2) = con.m
@inline TrajectoryOptimization.sense(con::NormConstraint2) = con.sense
# @inline Base.length(::NormConstraint2) = 1
@inline Base.length(::NormConstraint2{TrajectoryOptimization.SecondOrderCone,D}) where D = D + 1
@inline Base.length(::NormConstraint2) = 1 

function TrajectoryOptimization.evaluate(
                con::NormConstraint2{TrajectoryOptimization.SecondOrderCone},
                z::AbstractKnotPoint)
    v = con.A*z.z[con.inds]
    return push(v, con.c'*z.z[con.inds])
end

function TrajectoryOptimization.jacobian!(∇c,
                con::NormConstraint2{TrajectoryOptimization.SecondOrderCone},
                z::AbstractKnotPoint)
    ∇c[1:length(con.inds),con.inds] = con.A
    ∇c[1+length(con.inds),con.inds] = con.c
    return true
end

function TrajectoryOptimization.evaluate(
                con::NormConstraint2,
                z::AbstractKnotPoint)
    x = z.z[con.inds]
    v = con.A*x
    s = con.c'x
    return SA[v'v - s*s]
end

function TrajectoryOptimization.jacobian!(∇c,
                con::NormConstraint2,
                z::AbstractKnotPoint)
    A,c = SMatrix(con.A), con.c
    jac = A'A - c*c' 
    x = z.z[con.inds]
    jac = jac * x
    for (i,j) in enumerate(con.inds)
        ∇c[j] = 2*jac[i]
    end
    return false
end

function TrajectoryOptimization.change_dimension(con::NormConstraint2,
                n::Int, m::Int, ix=1:n, iu=1:m)
    NormConstraint2(n, m, con.val, con.sense, ix[con.inds])
end
