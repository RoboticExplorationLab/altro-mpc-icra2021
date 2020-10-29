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
        r_b = rot' * (r[i] - p)

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
    sprung_mass,
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
    sprung_mass,
)::SMatrix{n,m,T,n * m} where {T,n,m}
    return ForwardDiff.jacobian(
        (u_var) ->
            NonLinearContinuousDynamics(x, u_var, r, contacts, J, sprung_mass),
        u,
    )
end

# override dynamics evaluation to improve speed
function RD.discrete_dynamics(::Type{PassThrough}, model::LinearModel, x, u, t, dt)
    k = RD.get_k(model, t)
    dt_model = RD.is_timevarying(model) ? model.times[k+1] - model.times[k] : model.dt
    
    A = SMatrix(model.A[k])
    B = SMatrix(model.B[k])
    d = SVector(model.d[k])

    p_ind = SLegIndexToRange(1)
    ϕ_ind = SLegIndexToRange(2)
    v_ind = SLegIndexToRange(3)
    ω_ind = SLegIndexToRange(4)

    p₋ = x[p_ind]
    ϕ₋ = x[ϕ_ind]
    v₋ = x[v_ind]
    ω₋ = x[ω_ind]

    p₊ = p₋ + A[p_ind, v_ind]*v₋ + d[p_ind]
    ϕ₊ = ϕ₋ + A[ϕ_ind, ω_ind]*ω₋ + d[ϕ_ind]
    v₊ = v₋ + B[v_ind, SLegIndexToRange(1)]*u[SLegIndexToRange(1)] + B[v_ind, SLegIndexToRange(2)]*u[SLegIndexToRange(2)] +
            B[v_ind, SLegIndexToRange(3)]*u[SLegIndexToRange(3)] + B[v_ind, SLegIndexToRange(4)]*u[SLegIndexToRange(4)] + d[v_ind]
    ω₊ = ω₋ + B[ω_ind, SLegIndexToRange(1)]*u[SLegIndexToRange(1)] + B[ω_ind, SLegIndexToRange(2)]*u[SLegIndexToRange(2)] +
            B[ω_ind, SLegIndexToRange(3)]*u[SLegIndexToRange(3)] + B[ω_ind, SLegIndexToRange(4)]*u[SLegIndexToRange(4)] + d[ω_ind]

    xdot = [p₊; ϕ₊; v₊; ω₊]
    xdot
end

function reference_trajectory!(
    x_curr::AbstractVector{T},
    param::ControllerParams,
) where {T<:Number}
    # integrate the x,y position from the reference

    # interp_length = Int(round(param.N/3))
    interp_length = 2
    α = [collect(range(0, 1, length = interp_length)); ones(param.N - interp_length)]

    # for i = 1:param.N+1
    #     param.x_ref[i] = (1 - α[i]) * x_curr + α[i] * param.x_des
    # end

    p_integrate = x_curr[SUnitRange(1,2)]
    v_i = param.x_des[SUnitRange(7,8)]

    for i=1:param.N
        # if param.vel_ctrl
		# 	p_integrate += v_i*param.optimizer.dt

		# 	param.x_ref[i] = [p_integrate; ((1 - α[i]) * x_curr + α[i] * param.x_des)[SUnitRange(3, 12)]]
        # else
        #     param.x_ref[i] = ((1 - α[i]) * x_curr + α[i] * param.x_des)
        # end

        param.x_ref[i] = param.x_des
    end
end