"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

const n = 12
const m = 12

function update_osqp_model!(param::ControllerParams, x_curr)
    opt = param.optimizer
    N = param.N

    x_end = n*(N-1)

    for i = 1:(N-1)
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

        opt.A_vec[i] = oneunit(SMatrix{12,12}) + A_c_i * opt.dt #+ A_c_i^2 * opt.dt^2/2
        opt.B_vec[i] = B_c_i * opt.dt #+ A_c_i*B_c_i*opt.dt^2/2
        opt.d_vec[i] = d_c_i * opt.dt #+ A_c_i*d_c_i*opt.dt^2/2

        if i==1
            opt.A_osqp[select(1,n), x_end .+ select(1,m)] = opt.B_vec[i]
            opt.A_osqp[select(1,n), select(1,n)] = -I(n)
            opt.l[select(1,n)] = -opt.d_vec[i] - opt.A_vec[i]*x_curr
            opt.u[select(1,n)] = -opt.d_vec[i] - opt.A_vec[i]*x_curr
        else
            opt.A_osqp[select(i,n), select(i-1,n)] = opt.A_vec[i]
            opt.A_osqp[select(i,n), x_end .+ select(i,m)] = opt.B_vec[i]
            opt.A_osqp[select(i,n), select(i,n)] = -I(n)
            opt.l[select(i,n)] = -opt.d_vec[i]
            opt.u[select(i,n)] = -opt.d_vec[i]

            # x_ref:
            opt.x_ref_osqp[select(i-1,n)] = param.x_ref[N]

            # u_ref:
            opt.x_ref_osqp[x_end .+ select(i,m)] = opt.u_ref[i]
        end
    end

    opt.q_osqp .= -opt.P_osqp*opt.x_ref_osqp

    OSQP.update!(opt.model, q=opt.q_osqp, Ax=opt.A_osqp.nzval, l=opt.l, u=opt.u)
end

function foot_forces!(
    x_curr::AbstractVector{T},
    t::T,
    param::ControllerParams{O},
) where {T<:Number, O<:OSQPParams}
    # x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
    # contacts: 4xN+1 matrix of foot contacts over the planning horizon
    # foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

    N = param.N

    update_osqp_model!(param, x_curr)

    # gets benchmark to return before populating results
    results = OSQP.solve!(param.optimizer.model)

    param.new_info = true
    param.last_solve_time = results.info.solve_time
    param.last_solve_iterations = results.info.iter

    for i=1:N
        param.optimizer.X0[i] = SVector{12}(i == 1 ? x_curr : results.x[select12(i-1)])
        i==N && continue
        param.optimizer.U0[i] = SVector{12}(results.x[n*(N-1) .+ select12(i)])
    end

    param.forces = results.x[n*(N-1) .+ select12(1)]
end

function select12(i)
    return SVector{12}((12*(i-1)+1):(12*(i-1)+12))
end

function select12_3(i, j, k)
    return 12*(i-1) + 3*(j-1) + k
end
