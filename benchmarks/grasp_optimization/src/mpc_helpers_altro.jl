function altro_mpc_setup(o::SquareObject, X_warm, U_warm, N, shift=0)
    n, m = size(o)
    xf = X_warm[end]

    # objective
    Q = 1.0e-3*Diagonal(@SVector ones(n))
    Qf = 10.0*Diagonal(@SVector ones(n))
    R = 1.0*Diagonal(@SVector ones(m))
    obj = LQRObjective(Q,R,Qf,xf,N);

    # Update the Reference Trajectory
    for k in 1:N-1
        TO.set_LQR_goal!(obj.cost[k], X_warm[k], U_warm[k])
    end
    TO.set_LQR_goal!(obj.cost[end], X_warm[end])

    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,N)

    # Goal Constraint
    goal = GoalConstraint(SVector{n}(xf))
    add_constraint!(conSet, goal, N)

    for i = 1:N-1
        i_s = i + shift

        # Torque Balance
        B = [o.B[1][i_s] o.B[2][i_s]]
        t_bal = LinearConstraint2(n, m, B, [θdd[i_s],0,0], Equality(), u_ind)
        add_constraint!(conSet, t_bal, i:i)

        # Max Grasp Force
        A = zeros(2, m)
        A[1,1:Int(m/2)] = o.v[1][i_s]
        A[2,1+Int(m/2):end] = o.v[2][i_s]
        max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
        add_constraint!(conSet, max_f, i:i)

        # SOCP friction cone
        v1_i = o.v[1][i_s]
        A1 = (I - v1_i*v1_i')
        c1 = o.mu*v1_i
        nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
        add_constraint!(conSet, nc1, i:i)

        v2_i = o.v[2][i_s]
        A2 = (I - v2_i*v2_i')
        c2 = o.mu*v2_i
        nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
        add_constraint!(conSet, nc2, i:i)
    end

    return obj, conSet
end

function altro_mpc_update!(prob, o::SquareObject, X_warm, U_warm, N, shift)
    n, m = size(o)
    xf = X_warm[end]

    # Update the Reference Trajectory
    TO.set_initial_state!(prob, X_warm[1])
    for k in 1:N-1
        TO.set_LQR_goal!(prob.obj.cost[k], X_warm[k], U_warm[k])
    end
    TO.set_LQR_goal!(prob.obj.cost[end], X_warm[end])

    # Goal Constraint
    prob.constraints[1].xf .= xf

    i_c = 2 # constraint index
    for i = 1:N-1
        i_s = i + shift # index for v an B matrices

        # Torque Balance
        prob.constraints[i_c].A .= [o.B[1][i_s] o.B[2][i_s]]
        prob.constraints[i_c].b .= [θdd[i_s], 0, 0]
        i_c += 1

        # Max Grasp Force
        prob.constraints[i_c].A[1,1:Int(m/2)] .= o.v[1][i_s]
        prob.constraints[i_c].A[2,1+Int(m/2):end] .= o.v[2][i_s]
        i_c += 1

        # SOCP friction cone
        v1_i = o.v[1][i_s]
        prob.constraints[i_c].A .= (I - v1_i*v1_i')
        prob.constraints[i_c].c .= o.mu*v1_i
        i_c += 1

        v2_i = o.v[2][i_s]
        prob.constraints[i_c].A .= (I - v2_i*v2_i')
        prob.constraints[i_c].c .= o.mu*v2_i
        i_c += 1
    end

    initial_controls!(prob, U_warm)
    initial_states!(prob, X_warm)
    rollout!(prob)
    
    return
end
