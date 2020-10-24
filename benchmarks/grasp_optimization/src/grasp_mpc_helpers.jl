function mpc_update!(prob_mpc, o::SquareObject, k_mpc, Z_track)
    n, m, N_mpc = size(prob_mpc)

    # Update Initial Time
    dt = prob_mpc.Z[1].dt
    TO.set_initial_time!(prob_mpc, k_mpc*dt)

    # Update initial state by using 1st control, and adding some noise
    x0 = discrete_dynamics(TO.integration(prob_mpc), prob_mpc.model, prob_mpc.Z[1])
    x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
    TO.set_initial_state!(prob_mpc, x0)

    # Update the Reference Trajectory
    TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

    # Warm start
    RD.shift_fill!(prob_mpc.Z)

    # Setup ECOS
    prob_mpc_ecos, X_ecos, U_ecos = gen_ECOS(prob_mpc, k_mpc, Z_track)
    e_cons = prob_mpc_ecos.constraints
    e_cons += X_ecos[:,1] == x0 # Initial condition

    # Update ALTRO and ECOS Stage Constraints
    i_c = 1 # constraint index
    for i = 1:N_mpc-1
        # for convenience
        i_s = i + k_mpc # shift index for v and B matrices
        F1 = U_ecos[1:3, i] # contact force 1
        F2 = U_ecos[4:6, i] # contact force 2

        # Torque Balance
        A = [o.B[1][i_s] o.B[2][i_s]]
        b = [o.θdd[i_s], 0, 0]
        prob_mpc.constraints[i_c].A .= A
        prob_mpc.constraints[i_c].b .= b

        e_cons += A*U_ecos[:, i] == b

        i_c += 1

        # Max Grasp Force
        prob_mpc.constraints[i_c].A[1,1:Int(m/2)] .= o.v[1][i_s]
        prob_mpc.constraints[i_c].A[2,1+Int(m/2):end] .= o.v[2][i_s]

        A = copy(prob_mpc.constraints[i_c].A)
        e_cons += A*U_ecos[:, i] <= o.f*ones(2)

        i_c += 1

        # SOCP friction cone 1
        v1_i = o.v[1][i_s]
        prob_mpc.constraints[i_c].A .= (I - v1_i*v1_i')
        prob_mpc.constraints[i_c].c .= o.mu*v1_i

        A = copy(prob_mpc.constraints[i_c].A)
        c = copy(prob_mpc.constraints[i_c].c)
        e_cons += norm(A*F1) <= c'*F1

        i_c += 1

        # SOCP friction cone 2
        v2_i = o.v[2][i_s]
        prob_mpc.constraints[i_c].A .= (I - v2_i*v2_i')
        prob_mpc.constraints[i_c].c .= o.mu*v2_i

        A = copy(prob_mpc.constraints[i_c].A)
        c = copy(prob_mpc.constraints[i_c].c)
        e_cons += norm(A*F2) <= c'*F2

        i_c += 1

        # Dynamics Constraints for ECOS
        u = 1/o.mass * (F1 + F2) + o.g
        pos_ind = 1:3; vel_ind = 4:6
        e_cons += X_ecos[vel_ind, i+1] == X_ecos[vel_ind, i] + u*dt
        e_cons += X_ecos[pos_ind, i+1] == X_ecos[pos_ind, i] + X_ecos[vel_ind, i]*dt + u*.5*dt^2
    end

    return prob_mpc_ecos, U_ecos
end

function gen_ECOS(prob_mpc, k_mpc, Z_track)
    n, m, N_mpc = size(prob_mpc)

    # Variables
    X = Variable(n, N_mpc)
    U = Variable(m, N_mpc-1)

    # Warm Starting
    X_warm = states(prob_mpc.Z)[1:N_mpc]
    set_value!(X, hcat(Vector.(X_warm)...))

    U_warm = controls(prob_mpc.Z)[1:N_mpc-1]
    set_value!(U, hcat(Vector.(U_warm)...))

    # Get Tracking Trajectory
    X_ref = states(Z_track)[k_mpc .+ (1:N_mpc)]
    X_ref = hcat(Vector.(X_ref)...)

    U_ref = states(Z_track)[k_mpc .+ (1:N_mpc-1)]
    U_ref = hcat(Vector.(U_ref)...)

    # Set Tracking Objective
    Q = prob_mpc.obj[1].Q[1]
    Qf = prob_mpc.obj[1].R[1]
    R = prob_mpc.obj[end].Q[1]
    objective = Q*sumsquares(X[:,1:N_mpc-1]-X_ref[:, 1:N_mpc-1]) + Qf*sumsquares(X[:,N_mpc]-X_ref[:,N_mpc]) + R*sumsquares(U-U_ref)
    prob_mpc_ecos = minimize(objective)

    return prob_mpc_ecos, X, U
end
