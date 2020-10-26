# OSQP solver deps:
using OSQP
using SparseArrays

struct OSQPParams{T,S,M}
    dt::T
    N::S

    model::M
    x_ref_osqp::Vector{T}
    P_osqp::SparseMatrixCSC{T,S}
    q_osqp::Vector{T}
    A_osqp::SparseMatrixCSC{T,S}
    l::Vector{T}
    u::Vector{T}

    # Store dynamics matrices:
    A_vec::Vector{SMatrix{12,12,T,144}}
    B_vec::Vector{SMatrix{12,12,T,144}}
    d_vec::Vector{SVector{12,T}}

    u_ref::Vector{SVector{12, T}}
    J::SMatrix{3,3,T,9}
    sprung_mass::T
end

select(i, n) = (n*(i-1)+1):(n*(i-1)+n)
selectn_l(i, j, k, n=12, l=3) = n*(i-1) + l*(j-1) + k

function OSQPParams(
    dt::T,
    N::S,
    q::AbstractVector{T},
    r::AbstractVector{T},
    μ,
    min_vert_force,
    max_vert_force;
    n::Integer = 12,
    m::Integer = 12,
    tol = 1e-4
) where {T<:Number,S<:Integer}
    # initialize model and variables
    model = OSQP.Model()    
    # initialize quadratic cost parameters
    p = [repeat(q, N-1); repeat(r, N-1)]
    P_osqp = Diagonal(p)

    # x = [x(2), ..., x(N), u(1), ..., u(N-1)]
    osqp_dim = (m+n)*(N-1)
    x_ref_osqp = zeros(osqp_dim)
    u_ref = [zeros(m) for i=1:(N-1)]

    q_osqp = P_osqp*x_ref_osqp

    # (N-1)*12 Dynamics Constraints
    # (N-1)*20 Friction/Control Constraints
    constraint_dim = n*(N-1) + (N-1)*20
    l = zeros(constraint_dim)
    u = zeros(constraint_dim)
    C = zeros(constraint_dim, osqp_dim)

    A_ones = ones(n,n)
    B_ones = ones(n,m)
    d_ones = ones(n)

    x_end = n*(N-1)

    # Dynamics Constraints
    C[select(1,n), x_end .+ select(1,m)] = B_ones 
    C[select(1,n), select(1,n)] = -I(n)
    l[select(1,n)] = -d_ones - A_ones*ones(n)
    u[select(1,n)] = -d_ones - A_ones*ones(n)
    for i=2:N-1
        C[select(i,n), select(i-1,n)] = A_ones
        C[select(i,n), x_end .+ select(i,m)] = B_ones 
        C[select(i,n), select(i,n)] = -I(n)
        l[select(i,n)] = -d_ones
        u[select(i,n)] = -d_ones
    end

    ∞ = Inf
  
    # Control Constraints
    dyn_end = n*(N-1)
    for i=1:N-1
        # loop over each foot
        for j=1:4
            # |f_x| ≤ μf_z
            # -μf_z ≤ f_x ≤ μf_z
            # -∞ ≤ f_x - μf_z ≤ 0
            C[dyn_end .+ selectn_l(i, j, 1, 20, 5), x_end .+ selectn_l(i, j, 1)] = 1
            C[dyn_end .+ selectn_l(i, j, 1, 20, 5), x_end .+ selectn_l(i, j, 3)] = -μ
            l[dyn_end .+ selectn_l(i, j, 1, 20, 5)] = -∞
            u[dyn_end .+ selectn_l(i, j, 1, 20, 5)] = 0
            # 0 ≤ f_x + μf_z ≤ ∞
            C[dyn_end .+ selectn_l(i, j, 2, 20, 5), x_end .+ selectn_l(i, j, 1)] = 1
            C[dyn_end .+ selectn_l(i, j, 2, 20, 5), x_end .+ selectn_l(i, j, 3)] = μ
            l[dyn_end .+ selectn_l(i, j, 2, 20, 5)] = 0
            u[dyn_end .+ selectn_l(i, j, 2, 20, 5)] = ∞
            
            # |f_y| ≤ μf_z
            # -μf_z ≤ f_y ≤ μf_z
            # -∞ ≤ f_y - μf_z ≤ 0
            C[dyn_end .+ selectn_l(i, j, 3, 20, 5), x_end .+ selectn_l(i, j, 2)] = 1
            C[dyn_end .+ selectn_l(i, j, 3, 20, 5), x_end .+ selectn_l(i, j, 3)] = -μ
            l[dyn_end .+ selectn_l(i, j, 3, 20, 5)] = -∞
            u[dyn_end .+ selectn_l(i, j, 3, 20, 5)] = 0
            # 0 ≤ f_y + μf_z ≤ ∞
            C[dyn_end .+ selectn_l(i, j, 4, 20, 5), x_end .+ selectn_l(i, j, 2)] = 1
            C[dyn_end .+ selectn_l(i, j, 4, 20, 5), x_end .+ selectn_l(i, j, 3)] = μ
            l[dyn_end .+ selectn_l(i, j, 4, 20, 5)] = 0
            u[dyn_end .+ selectn_l(i, j, 4, 20, 5)] = ∞
            
            # min_vert_force ≤ f_z ≤ max_vert_force
            C[dyn_end .+ selectn_l(i, j, 5, 20, 5), x_end .+ selectn_l(i, j, 3)] = 1
            l[dyn_end .+ selectn_l(i, j, 5, 20, 5)] = min_vert_force
            u[dyn_end .+ selectn_l(i, j, 5, 20, 5)] = max_vert_force
        end
    end

    A_osqp = sparse(C)
    P_osqp = sparse(P_osqp)

    OSQP.setup!(model, P=P_osqp, q=q_osqp, A=A_osqp, l=l, u=u, verbose=false)
    OSQP.update_settings!(model, eps_abs = tol, eps_rel = tol, eps_prim_inf = tol, eps_dual_inf = tol)

    J = woofer.inertial.body_inertia
    sprung_mass = woofer.inertial.sprung_mass

    M = typeof(model)
    R = typeof(results)

    A_vec = [@SMatrix zeros(n,n) for i=1:(N-1)]
    B_vec = [@SMatrix zeros(n,m) for i=1:(N-1)]
    d_vec = [@SVector zeros(n) for i=1:(N-1)]

    return OSQPParams{T,S,M}(
        dt,
        N,
        model,
        x_ref_osqp,
        P_osqp,
        q_osqp,
        A_osqp,
        l,
        u,
        A_vec,
        B_vec,
        d_vec,
        u_ref,
        J,
        sprung_mass
    )
end