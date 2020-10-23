function GraspProblem(N = 31, tf = 3.0;
        mu = 0.5,     # coefficient of friction
        mass = 0.2,   # kg
        inertia = 1.,
        f_max = 3.0   # max grasp force
    )


    N = 31          # horizon
    tf = 3.         # final time
    dt = tf/(N-1)   # time step
    n = 6           # state size
    m = 6           # control size

    g = @SVector [0, 0, -9.81]  # gravity
    mu = .5                     # friction constant
    mass = .2                   # mass

    # rotational trajectory
    θ0 = 0; θf= pi/4; θd0 = 0; θdf = .15; t0 = 0;
    c = compute_rot_traj_coeffs(t0, tf, [θ0; θf; θd0; θdf])
    θ = [dot(c, [t^3,t^2,t,1]) for t = 0:dt:tf]
    θdd = [dot(c, [6t,2,0,0]) for t = 0:dt:tf]

    # generate p v B matrices
    p1_0 = [.0,1, 0]; v1_0 = [.0,-1, 0]
    p2_0 = [.0,-1, 0]; v2_0 = [.0,1, 0]
    p2_0 = [.0,1, 0]; v2_0 = [.0,-1, 0]
    p1_0 = [.0,-1, 0]; v1_0 = [.0,1, 0]
    p1, v1, B1 = generate_pvB_3D(p1_0, v1_0, θ)
    p2, v2, B2 = generate_pvB_3D(p2_0, v2_0, θ)

    # model
    o = SquareObject(n, m, mu, mass, inertia, f_max, g, [p1, p2], [v1, v2], [B1, B2])

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
    obj = LQRObjective(Q,R,Qf,xf,N);
    
    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,N)
    
    # Goal Constraint
    goal = GoalConstraint(SVector{n}(xf))
    add_constraint!(conSet, goal, N)
    
    # Stage Constraints
    for i = 1:N-1
        # global model, u_ind, F1_ind, F2_ind, n, m
        # local B, t_bal, A, max_f, A1, c1, nc1, A2, c2, nc2
    
        # Torque Balance
        B = [o.B[1][i] o.B[2][i]]
        t_bal = LinearConstraint2(n, m, B, [θdd[i],0,0], Equality(), u_ind)
        add_constraint!(conSet, t_bal, i:i)
    
        # Max Grasp Force
        A = zeros(2, m)
        A[1,1:Int(m/2)] = o.v[1][i]
        A[2,1+Int(m/2):end] = o.v[2][i]
        max_f = LinearConstraint(n, m, A, o.f*ones(2), Inequality(), u_ind)
        add_constraint!(conSet, max_f, i:i)
    
        # SOCP friction cone
        v1_i = o.v[1][i]
        A1 = (I - v1_i*v1_i')
        c1 = o.mu*v1_i
        nc1 = FrictionConstraint(n, m, A1, c1, TO.SecondOrderCone(), F1_ind)
        add_constraint!(conSet, nc1, i:i)
    
        v2_i = o.v[2][i]
        A2 = (I - v2_i*v2_i')
        c2 = o.mu*v2_i
        nc2 = FrictionConstraint(n, m, A2, c2, TO.SecondOrderCone(), F2_ind)
        add_constraint!(conSet, nc2, i:i)
    end
    
    # Problem
    prob = Problem(o, obj, xf, tf, x0=SVector{n}(x0), constraints=conSet);
    
    u0 = @SVector [0, -1.5, o.mass*9.81/2, 0, 1.5, o.mass*9.81/2]
    U0 = [u0 for k = 1:N-1]
    initial_controls!(prob, U0)
    rollout!(prob);

    return prob
end



