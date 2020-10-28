function GraspProblem(o::SquareObject, N = 61, tf = 6.0,
                    x0 = [0.,3.,3.,0.,0.,0.], # initial position
                    xf = zeros(6) # final position
                    )

    n, m = size(o)  # state and control size
    dt = tf/(N-1)   # time step

    # set SquareObject orientation trajectories if blank
    o.p == [] && set_orientation_traj!(o, dt, tf)

    # Indices for convenience
    pos_ind = 1:Int(n/2)
    vel_ind = 1+Int(n/2):n
    u_ind = n .+ (1:m)
    F1_ind = n .+ (1:Int(m/2))
    F2_ind = n .+ (1+Int(m/2):m)

    # Objective
    Q = 1.0e-3*Diagonal(@SVector ones(n))
    Qf = 10.0*Diagonal(@SVector ones(n))
    R = 1.0*Diagonal(@SVector ones(m))
    obj = LQRObjective(Q,R,Qf,xf,N);

    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,N)

    # Goal Constraint
    goal = GoalConstraint(SVector{n}(xf))
    TO.add_constraint!(conSet, goal, N)

    # Stage Constraints
    Atorque = [[o.B[1][i] o.B[2][i]] for i = 1:N-1]
    btorque = [[o.θdd[i],0,0] for i = 1:N-1]
    torque_balance = LinearConstraintTraj(n, m, Atorque, btorque, Equality(), u_ind)
    add_constraint!(conSet, torque_balance, 1:N-1)

    Agrasp = map(1:N-1) do i
        A = zeros(2,m)
        A[1,1:m÷2] = o.v[1][i]
        A[2,1+m÷2:end] = o.v[2][i]
        A
    end
    bgrasp = [o.f * ones(2) for i = 1:N-1]
    max_force = LinearConstraintTraj(n, m, Agrasp, bgrasp, Inequality(), u_ind)
    add_constraint!(conSet, max_force, 1:N-1)
    
    for i = 1:N-1
        # Torque Balance
        B = [o.B[1][i] o.B[2][i]]
        t_bal = LinearConstraint2(n, m, B, [o.θdd[i],0,0], Equality(), u_ind)
        # TO.add_constraint!(conSet, t_bal, i:i)

        # Max Grasp Force
        A = zeros(2, m)
        A[1,1:Int(m/2)] = o.v[1][i]
        A[2,1+Int(m/2):end] = o.v[2][i]
        max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
        # TO.add_constraint!(conSet, max_f, i:i)

        # SOCP Friction Cone 1
        v1_i = o.v[1][i]
        A1 = (I - v1_i*v1_i')
        c1 = o.mu*v1_i
        nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
        TO.add_constraint!(conSet, nc1, i:i)

        # SOCP Friction Cone 2
        v2_i = o.v[2][i]
        A2 = (I - v2_i*v2_i')
        c2 = o.mu*v2_i
        nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
        TO.add_constraint!(conSet, nc2, i:i)
    end

    # Problem
    prob = TO.Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet, integration=RD.PassThrough);

    # Intialize Trajectory
    u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
    U0 = [u0 for k = 1:N-1]
    initial_controls!(prob, U0)
    rollout!(prob);

    return prob
end
