"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

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
        if param.vel_ctrl
			p_integrate += v_i*param.optimizer.dt

			param.x_ref[i] = [p_integrate; ((1 - α[i]) * x_curr + α[i] * param.x_des)[SUnitRange(3, 12)]]
        else
            param.x_ref[i] = ((1 - α[i]) * x_curr + α[i] * param.x_des)
        end
    end
end

function update_dynamics_matrices!(param::ControllerParams)
    opt = param.optimizer

    for i = 1:(param.N)
        A_c_i = LinearizedContinuousDynamicsA(
            param.x_ref[i],
            opt.u_ref[i],
            param.foot_locs[i],
            param.contacts[i],
            opt.J,
            opt.sprung_mass,
        )
        B_c_i = LinearizedContinuousDynamicsB(
            param.x_ref[i],
            opt.u_ref[i],
            param.foot_locs[i],
            param.contacts[i],
            opt.J,
            opt.sprung_mass,
        )
        d_c_i = - A_c_i*param.x_ref[i] - B_c_i*opt.u_ref[i] + NonLinearContinuousDynamics(
            param.x_ref[i],
            opt.u_ref[i],
            param.foot_locs[i],
            param.contacts[i],
            opt.J,
            opt.sprung_mass,
        ) 

        # Midpoint Discretization
        opt.model.A[i] = oneunit(SMatrix{12,12}) + A_c_i * opt.dt + A_c_i^2 * opt.dt^2/2
        opt.model.B[i] = B_c_i * opt.dt + A_c_i*B_c_i*opt.dt^2/2
        opt.model.d[i] = d_c_i * opt.dt + A_c_i*d_c_i*opt.dt^2/2
    end

    # Z = Traj(param.x_ref, opt.u_ref, opt.dt)

    # update_trajectory!(opt.objective, Z)
end

function foot_forces!(
    x_curr::AbstractVector{T},
    t::T,
    param::ControllerParams{O},
) where {T<:Number, O<:AltroParams}
    # x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
    # contacts: 4xN+1 matrix of foot contacts over the planning horizon
    # foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon
    opt = param.optimizer

    tf = (param.N - 1) * opt.dt
    n = 12
    m = 12

    # TrajectoryOptimization.reset!(opt.constraints) 
    update_dynamics_matrices!(param)
    opt.solver.solver_al.solver_uncon.x0 .= x_curr

    initial_states!(opt.problem, opt.X0)
    initial_controls!(opt.problem, opt.U0)
    b = @benchmark Altro.solve!($(opt.solver)) samples=1 evals=1

    opt.X0 .= states(opt.solver)
    opt.U0 .= controls(opt.solver)

    if status(opt.solver) != Altro.SOLVE_SUCCEEDED
        @warn "Solver status: $(opt.solver.stats.status)"
    end

    param.forces = opt.U0[1]

    return b
end
