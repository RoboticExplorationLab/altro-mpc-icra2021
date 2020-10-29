function mpc_update!(prob_mpc, o::SquareObject, k_mpc, Z_track)
    n, m, N_mpc = size(prob_mpc)

    # Update Initial Time
    dt = prob_mpc.Z[1].dt
    TO.set_initial_time!(prob_mpc, k_mpc*dt)

    # Update Initial State (use 1st control and add some noise)
    x0 = discrete_dynamics(TO.integration(prob_mpc), prob_mpc.model, prob_mpc.Z[1])
    x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
    TO.set_initial_state!(prob_mpc, x0)

    # Update the Reference Trajectory
    TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc+1)

    # Warm start
    RD.shift_fill!(prob_mpc.Z)

    # Setup ECOS
    prob_mpc_ecos, X_ecos, U_ecos = gen_ECOS(prob_mpc, k_mpc+1, Z_track)
    con = @constraint(prob_mpc_ecos, X_ecos[:,1] .== x0) # Initial condition
    set_name(con, "initial_condition")

    # Update ALTRO and ECOS Stage Constraints
    cons = prob_mpc.constraints
    traj_inds = k_mpc .+ (1:N_mpc)

    # Torque Balance
    Atorque = [[o.B[1][i] o.B[2][i]] for i in traj_inds]
    btorque = [[o.θdd[i],0,0] for i in traj_inds]

    # Friction 1
    Afriction1 = map(traj_inds) do i
        v = o.v[1][i]
        (I - v*v')
    end
    cfriction1 = [o.mu * o.v[1][i] for i in traj_inds]

    # Friction 2
    Afriction2 = map(traj_inds) do i
        v = o.v[2][i]
        (I - v*v')
    end
    cfriction2 = [o.mu * o.v[2][i] for i in traj_inds]

    # Update the constraints
    for i = 1:N_mpc-1
        cons[1].A[i] = Atorque[i]
        cons[1].b[i] = btorque[i]
        cons[2].A[i][1,1:m÷2]     .= o.v[1][i+k_mpc]
        cons[2].A[i][2,1+m÷2:end] .= o.v[2][i+k_mpc]
        cons[3].A[i] = Afriction1[i]
        cons[3].c[i] = cfriction1[i]
        cons[4].A[i] = Afriction2[i]
        cons[4].c[i] = cfriction2[i]
    end


    i_c = 1 # constraint index
    for i = 1:N_mpc-1
        i_s = i + k_mpc # shift index for v and B matrices

        # for convenience with JuMP
        F1 = @expression(prob_mpc_ecos, U_ecos[1:3, i]) # contact force 1
        F2 = @expression(prob_mpc_ecos, U_ecos[4:6, i]) # contact force 2

        # Torque Balance
        A = [o.B[1][i_s] o.B[2][i_s]]
        b = [o.θdd[i_s], 0, 0]
        # prob_mpc.constraints[i_c].A .= A
        # prob_mpc.constraints[i_c].b .= b
        # i_c += 1

        @constraint(prob_mpc_ecos, Matrix(Atorque[i])*U_ecos[:, i] .== Vector(btorque[i]))

        # Max Grasp Force
        # prob_mpc.constraints[i_c].A[1,1:Int(m/2)] .= o.v[1][i_s]
        # prob_mpc.constraints[i_c].A[2,1+Int(m/2):end] .= o.v[2][i_s]
        # A = copy(prob_mpc.constraints[i_c].A)
        # i_c += 1
        A = Matrix((cons[2].A[i]))
        @constraint(prob_mpc_ecos, A*U_ecos[:, i] .<= o.f*ones(2))


        # SOCP Friction Cone 1
        # v1_i = o.v[1][i_s]
        # prob_mpc.constraints[i_c].A .= (I - v1_i*v1_i')
        # prob_mpc.constraints[i_c].c .= o.mu*v1_i

        # A = copy(prob_mpc.constraints[i_c].A)
        # c = copy(prob_mpc.constraints[i_c].c)
        # i_c += 1
        A = Matrix(cons[3].A[i])
        c = Vector(cons[3].c[i])
        @constraint(prob_mpc_ecos, [c'*F1; A[1:end,1:end]*F1] in JuMP.SecondOrderCone())


        # SOCP Friction Cone 2
        # v2_i = o.v[2][i_s]
        # prob_mpc.constraints[i_c].A .= (I - v2_i*v2_i')
        # prob_mpc.constraints[i_c].c .= o.mu*v2_i
        # i_c += 1

        # A = copy(prob_mpc.constraints[i_c].A)
        # c = copy(prob_mpc.constraints[i_c].c)
        A = Matrix(cons[4].A[i])
        c = Vector(cons[4].c[i])
        @constraint(prob_mpc_ecos, [c'*F2; A[1:end,1:end]*F2] in JuMP.SecondOrderCone())


        # Dynamics Constraints for ECOS
        pos_ind = 1:3; vel_ind = 4:6 # indices
        a = 1/o.mass * (F1 + F2) + o.g # xyz acceleration
        @constraint(prob_mpc_ecos, X_ecos[vel_ind, i+1] .== X_ecos[vel_ind, i] + a*dt)
        @constraint(prob_mpc_ecos, X_ecos[pos_ind, i+1] .== X_ecos[pos_ind, i] + X_ecos[vel_ind, i]*dt + .5*a*dt^2)
    end

    return prob_mpc_ecos, X_ecos, U_ecos
end

function gen_ECOS(prob_mpc, k_mpc, Z_track)
    n, m, N_mpc = size(prob_mpc) # ALTRO problem dimensions
    prob_mpc_ecos = Model() # JuMP Model for ECOS

    # Variables
    @variable(prob_mpc_ecos, X[1:n, 1:N_mpc])
    @variable(prob_mpc_ecos, U[1:m, 1:N_mpc-1])

    # Get Tracking Trajectories
    X_ref = states(Z_track)[(k_mpc-1) .+ (1:N_mpc)]
    X_ref = hcat(Vector.(X_ref)...)

    U_ref = controls(Z_track)[(k_mpc-1) .+ (1:N_mpc-1)]
    U_ref = hcat(Vector.(U_ref)...)

    # Set Tracking Objective
    Q = prob_mpc.obj[1].Q[1]
    R = prob_mpc.obj[1].R[1]
    Qf = prob_mpc.obj[end].Q[1]

    state = @expression(prob_mpc_ecos, Q*sum((X[:,i]-X_ref[:,i])'*(X[:,i]-X_ref[:,i]) for i in 1:N_mpc-1))
    control = @expression(prob_mpc_ecos, R*sum((U[:,i]-U_ref[:,i])'*(U[:,i]-U_ref[:,i]) for i in 1:N_mpc-1))
    terminal = @expression(prob_mpc_ecos, Qf*(X[:,end]-X_ref[:,end])'*(X[:,end]-X_ref[:,end]))

    @objective(prob_mpc_ecos, Min, state+control+terminal)

    return prob_mpc_ecos, X, U
end
