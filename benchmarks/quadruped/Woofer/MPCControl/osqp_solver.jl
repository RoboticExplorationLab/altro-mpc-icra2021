"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

using ForwardDiff

function NonLinearContinuousDynamics(
    x::SVector,
    u::SVector,
    r::FootstepLocation,
    contacts::SVector,
    J::SMatrix,
    sprung_mass::AbstractFloat,
)
    rot = MRP(x[4], x[5], x[6])

    p = x[SUnitRange(1,3)]
    v = x[SUnitRange(7,9)]
    ω = x[SUnitRange(10,12)]

    x_d_1_3 = v
    x_d_4_6 = Rotations.kinematics(rot, ω)

    torque_sum = @SVector zeros(3)
    force_sum = @SVector [0, 0, -9.81]
    for i = 1:4
        force_sum += 1 / sprung_mass * contacts[i] * u[SLegIndexToRange(i)]

        # foot position in body frame:
        r_b = rot'*(r[i] - p)

        torque_sum +=
            contacts[i] *
            Rotations.skew(r_b) *
            rot' *
            u[SLegIndexToRange(i)]
    end
    x_d_7_9 = force_sum
    x_d_10_12 = inv(J) * (-Rotations.skew(ω) * J * ω + torque_sum)

    return [x_d_1_3; x_d_4_6; x_d_7_9; x_d_10_12]
end

function LinearizedContinuousDynamicsA(
    x::SVector{n,T},
    u::SVector{m,T},
    r,
    contacts,
    J,
    sprung_mass
)::SMatrix{n,n,T,n * n} where {T,n,m}
    return ForwardDiff.jacobian(
        (x_var) ->
            NonLinearContinuousDynamics(x_var, u, r, contacts, J, sprung_mass),
        x,
    )
end

function LinearizedContinuousDynamicsB(
    x::SVector{n,T},
    u::SVector{m,T},
    r,
    contacts,
    J,
    sprung_mass
)::SMatrix{n,m,T,n * m} where {T,n,m}
    return ForwardDiff.jacobian(
        (u_var) ->
            NonLinearContinuousDynamics(x, u_var, r, contacts, J, sprung_mass),
        u,
    )
end

function reference_trajectory!(
    x_curr::AbstractVector{T},
    param::ControllerParams,
) where {T<:Number}
    α = collect(range(0, 1, length = param.N + 1))

    for i = 1:param.N+1
        param.x_ref[i] = (1 - α[i]) * x_curr + α[i] * param.x_des
    end
end

function foot_forces!(
    x_curr::AbstractVector{T},
    t::T,
    param::ControllerParams,
) where {T<:Number}
    # x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
    # contacts: 4xN+1 matrix of foot contacts over the planning horizon
    # foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

    opt = param.optimizer
    N = opt.N
    dt = opt.dt

    for i = 1:N
        opt.x_ref[i][] = param.x_ref[i]

        A_c_i = LinearizedContinuousDynamicsA(opt.x_ref[i][], opt.u_ref[i][], param.foot_locs[i], param.contacts[i], opt.J, opt.sprung_mass)
		B_c_i = LinearizedContinuousDynamicsB(opt.x_ref[i][], opt.u_ref[i][], param.foot_locs[i], param.contacts[i], opt.J, opt.sprung_mass)
		d_c_i = NonLinearContinuousDynamics(opt.x_ref[i][], opt.u_ref[i][], param.foot_locs[i], param.contacts[i], opt.J, opt.sprung_mass)

		# Euler Discretization
		opt.A_d[i][] = oneunit(SMatrix{12,12,T}) + A_c_i*dt
		opt.B_d[i][] = B_c_i*dt
		opt.d_d[i][] = d_c_i*dt + opt.x_ref[i][]

		opt.q[i][] = 2*opt.Q_i*opt.x_ref[i][]
		opt.r[i][] = 2*opt.R_i*opt.u_ref[i][]
    end
    opt.x_ref[N+1][] = param.x_ref[N+1]
	opt.q[N+1][] = 2*opt.Q_f*opt.x_ref[N+1][]

    allocs = @allocated(solve!(opt.model))

    # println("Allocations: ", allocs)
    
    rot = MRP(x_curr[4], x_curr[5], x_curr[6])

    forces_inertial = value.(opt.model, opt.u)[select12(1)]

    param.forces = [rot \ forces_inertial[SLegIndexToRange(1)];
                    rot \ forces_inertial[SLegIndexToRange(2)];
                    rot \ forces_inertial[SLegIndexToRange(3)];
                    rot \ forces_inertial[SLegIndexToRange(4)]]
end

function select12(i)
    return SVector{12}((12*(i-1)+1):(12*(i-1)+12))
end

function select12_3(i, j, k)
    return 12*(i-1) + 3*(j-1) + k
end
