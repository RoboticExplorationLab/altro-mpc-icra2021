"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""
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

    param.new_info = true
    param.last_solve_time = b.times[1]

    opt.X0 .= states(opt.solver)
    opt.U0 .= controls(opt.solver)

    if status(opt.solver) != Altro.SOLVE_SUCCEEDED
        @warn "Solver status: $(opt.solver.stats.status)"

        @show opt.solver.stats
    end

    param.forces = opt.U0[1]
end
