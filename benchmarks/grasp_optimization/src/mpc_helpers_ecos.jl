function ecos_mpc_setup(o::SquareObject, X_ref, U_ref, N, shift, X_warm, U_warm)
    x0 = X_ref[:, 1]
    xf = X_ref[:, end]

    # Variables
    Z = Variable(n, N)
    F1 = Variable(Int(m/2), N-1)
    F2 = Variable(Int(m/2), N-1)

    # Warm start
    set_value!(Z, X_warm)
    set_value!(F1, U_warm[1:Int(o.n/2),:])
    set_value!(F2, U_warm[1+Int(o.n/2):end,:])

    # Objective
    Q = 1.0e-3
    Qf = 10.0
    R = 1.0
    Z_cold = X_cold[:, shift .+ (1:N)]
    objective = Q*sumsquares(Z[:,1:N-1]-X_ref[:, 1:N-1]) + Qf*sumsquares(Z[:,N]-xf) + R*sumsquares([F1;F2]-U_ref)
    prob = minimize(objective)

    # Start and Goal Constraint
    prob.constraints += Z[:,1] == x0
    prob.constraints += Z[:,N] == xf

    # Stage Constraints
    for t = 1:N-1
        t_s = t + shift

        # Torque Balance
        prob.constraints += [θdd[t_s],0,0] == o.B[1][t_s] * F1[:, t] + o.B[2][t_s] * F2[:, t]

        # Max Grasp Force
        prob.constraints += o.v[1][t_s]'*F1[:, t] <= o.f
        prob.constraints += o.v[2][t_s]'*F2[:, t] <= o.f

        # SOCP friction cone
        prob.constraints += norm((I - o.v[1][t_s]*o.v[1][t_s]')*F1[:, t]) <= o.mu*o.v[1][t_s]'*F1[:, t] # friction cone
        prob.constraints += norm((I - o.v[2][t_s]*o.v[2][t_s]')*F2[:, t]) <= o.mu*o.v[2][t_s]'*F2[:, t] # friction cone

        # Dynamics
        u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
        prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
        prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2
    end

    return prob, Z, F1, F2
end
