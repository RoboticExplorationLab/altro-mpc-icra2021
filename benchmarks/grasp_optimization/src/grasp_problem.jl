function GraspProblem(o::SquareObject, N = 51, tf = 5.0)
    n, m = size(o)
    dt = tf/(N-1)   # time step

    # set SquareObject orientation trajectories if blank
    o.p == [] && set_orientation_traj!(o, dt, tf)

    # initial and final positions
    x0 = [0.,3.,3.,0.,0.,0.]
    xf = zeros(n)

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
    # TODO set flag for tracking vs not obj = TO.TrackingObjective(Q, R, Z, Qf=Qf)
    obj = LQRObjective(Q,R,Qf,xf,N);

    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,N)

    # Goal Constraint
    goal = GoalConstraint(SVector{n}(xf))
    TO.add_constraint!(conSet, goal, N)

    # Stage Constraints
    for i = 1:N-1
        # Torque Balance
        B = [o.B[1][i] o.B[2][i]]
        t_bal = LinearConstraint2(n, m, B, [o.θdd[i],0,0], Equality(), u_ind)
        TO.add_constraint!(conSet, t_bal, i:i)

        # Max Grasp Force
        A = zeros(2, m)
        A[1,1:Int(m/2)] = o.v[1][i]
        A[2,1+Int(m/2):end] = o.v[2][i]
        max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
        TO.add_constraint!(conSet, max_f, i:i)

        # SOCP friction cone
        v1_i = o.v[1][i]
        A1 = (I - v1_i*v1_i')
        c1 = o.mu*v1_i
        nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
        TO.add_constraint!(conSet, nc1, i:i)

        v2_i = o.v[2][i]
        A2 = (I - v2_i*v2_i')
        c2 = o.mu*v2_i
        nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
        TO.add_constraint!(conSet, nc2, i:i)
    end

    # Problem
    prob = TO.Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet);

    u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
    U0 = [u0 for k = 1:N-1]
    initial_controls!(prob, U0)
    rollout!(prob);

    return prob
end
