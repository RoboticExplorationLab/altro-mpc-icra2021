"""
This file turns the discrete MPC problem into a quadratic problem via a sparse
formulation.
"""

function update_dynamics_matrices_ecos!(param::ControllerParams)
    opt = param.optimizer
    N = param.N

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

        opt.A_vec[i] = oneunit(SMatrix{12,12}) + A_c_i * opt.dt + A_c_i^2 * opt.dt^2/2
        opt.B_vec[i] = B_c_i * opt.dt + A_c_i*B_c_i*opt.dt^2/2
        opt.d_vec[i] = d_c_i * opt.dt + A_c_i*d_c_i*opt.dt^2/2
    end
end

function foot_forces!(
    x_curr::AbstractVector{T},
    t::T,
    param::ControllerParams{O},
) where {T<:Number, O<:ECOSParams}
    # x_ref: 12xN+1 matrix of state reference trajectory (where first column is x0)
    # contacts: 4xN+1 matrix of foot contacts over the planning horizon
    # foot_locs: 12xN+1 matrix of foot location in body frame over planning horizon

    N = param.N
    opt = param.optimizer

    update_dynamics_matrices_ecos!(param)

    jump_model = Model(
        optimizer_with_attributes(
            ECOS.Optimizer, "feastol"=>1e-4, "abstol"=>1e-4, "reltol"=>1e-4
        )
    )
    set_silent(jump_model)

    X = [@variable(jump_model, [i=1:n]) for i=1:N]
    U = [@variable(jump_model, [i=1:m]) for i=1:N-1]

    @objective(jump_model, Min, sum([(X[i] - param.x_ref[i])'*opt.Q*(X[i] - param.x_ref[i]) for i=1:N]) + 
                    sum([(U[i] - opt.u_ref[i])'*opt.R*(U[i] - opt.u_ref[i]) for i=1:N-1]))

    @constraint(jump_model, X[1] .== x_curr)
    for i=1:N-1
        @constraint(jump_model, X[i+1] .== opt.A_vec[i]*X[i] + opt.B_vec[i]*U[i] + opt.d_vec[i])
        @constraint(jump_model, [opt.μ * U[i][3]; U[i][1:2]] in SecondOrderCone())
        @constraint(jump_model, 0 <= U[i][3])
        @constraint(jump_model, opt.max_vert_force >= U[i][3])
    end

    # gets benchmark to return before populating results
    optimize!(jump_model)
    b = solve_time(jump_model)
    param.forces = value.(U[1])

    return b
end

function select12(i)
    return SVector{12}((12*(i-1)+1):(12*(i-1)+12))
end

function select12_3(i, j, k)
    return 12*(i-1) + 3*(j-1) + k
end
