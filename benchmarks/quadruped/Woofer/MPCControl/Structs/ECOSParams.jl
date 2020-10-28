# ECOS Dependencies
using JuMP
import ECOS

struct ECOSParams{T,S}
    dt::T
    N::S

    Q::Diagonal{T}
    R::Diagonal{T}
    μ::T
    max_vert_force::T

    # Store dynamics matrices:
    A_vec::Vector{SMatrix{12,12,T,144}}
    B_vec::Vector{SMatrix{12,12,T,144}}
    d_vec::Vector{SVector{12,T}}

    X0::Vector{SVector{12,T}}
    U0::Vector{SVector{12,T}}

    tol::T

    u_ref::Vector{SVector{12, T}}
    J::SMatrix{3,3,T,9}
    sprung_mass::T
end

function ECOSParams(
    dt::T,
    N::S,
    q::AbstractVector{T},
    r::AbstractVector{T},
    μ,
    min_vert_force,
    max_vert_force;
    tol = 1e-4
) where {T<:Number,S<:Integer}
    J = woofer.inertial.body_inertia
    sprung_mass = woofer.inertial.sprung_mass

    A_vec = [@SMatrix zeros(n,n) for i=1:(N-1)]
    B_vec = [@SMatrix zeros(n,m) for i=1:(N-1)]
    d_vec = [@SVector zeros(n) for i=1:(N-1)]

    X0 = [@SVector zeros(n) for i=1:N]
    U0 = [@SVector zeros(m) for i=1:N-1]

    u_ref = [@SVector zeros(m) for i=1:(N-1)]

    return ECOSParams{T,S}(
        dt,
        N,
        Diagonal(q),
        Diagonal(r),
        μ,
        max_vert_force,
        A_vec,
        B_vec,
        d_vec,
        X0, 
        U0,
        tol,
        u_ref,
        J,
        sprung_mass
    )
end