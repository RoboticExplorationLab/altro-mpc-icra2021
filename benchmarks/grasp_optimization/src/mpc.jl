using TrajectoryOptimization
const TO = TrajectoryOptimization

function MPC_SetUp(o::SquareObject, x0, xf, N, shift)
    n, m = size(o)

    # objective
    Q = 1.0e-3*Diagonal(@SVector ones(n))
    Qf = 10.0*Diagonal(@SVector ones(n))
    R = 1.0*Diagonal(@SVector ones(m))
    obj = LQRObjective(Q,R,Qf,xf,N);

    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,N)

    # Goal Constraint
    goal = GoalConstraint(SVector{n}(xf))
    add_constraint!(conSet, goal, N)

    include("src/new_constraints.jl")
    for i = 1:N-1
        i_s = i + shift

        # Torque Balance
        B = [o.B[1][i_s] o.B[2][i_s]]
        t_bal = LinearConstraint(n, m, B, [θdd[i_s],0,0], Equality(), u_ind)
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

function changeGoal!(c::ConstraintList, xf_new::SArray{Tuple{6},Float64,1,6})
    @assert typeof(c.constraints[end]) == GoalConstraint{6, Float64}
    c.constraints[end] = GoalConstraint(xf_new)

    return nothing
end


# Array{TrajectoryOptimization.AbstractConstraint,1}

function changeGoal!(p::TO.Problem, xf_new::SArray{Tuple{6},Float64,1,6})
    return changeGoal!(p.constraints, xf_new)
end

function warmstart_problem(model, obj::Objective, cset::ConstraintList, warm_states, warm_controls, tf)

    prob = TO.Problem(model, obj, warm_states[end], tf, x0=warm_states[1], constraints=cset)

    initial_controls!(prob, warm_controls)
    initial_states!(prob, warm_states)
    rollout!(prob)

    return prob
end

function MPC_Solve(prob::TO.Problem, opts; summary::Bool = false)
    altro = ALTROSolver(prob, opts)

    set_options!(altro, show_summary=false)
    b = benchmark_solve!(altro)

    set_options!(altro, show_summary=summary)
    solve!(altro);

    return b, altro
end

function ecos_mpc_step(o::SquareObject, x0, xf, N, shift)
    # variables
    F1 = Variable(Int(m/2), N-1)
    F2 = Variable(Int(m/2), N-1)
    Z = Variable(n, N)

    # objective
    Q = 1.0e-3
    Qf = 10.0
    R = 1.0
    Z_cold = X_cold[:, shift .+ (1:N)]
    objective = Q*sumsquares(Z[:,1:N-1]-Z_cold[:, 1:N-1]) + Qf*sumsquares(Z[:,N]-Z_cold[:, N]) + R*sumsquares([F1;F2])
    prob = minimize(objective)

    # start and goal constraint
    prob.constraints += Z[:,1] == x0
    prob.constraints += Z[:,N] == xf

    # stage constraints
    for t = 1:N-1
        t_s = t + shift

        # torque balance
        prob.constraints += [θdd[t_s],0,0] == o.B[1][t_s] * F1[:, t] + o.B[2][t_s] * F2[:, t]

        # max grasp force
        prob.constraints += o.v[1][t_s]'*F1[:, t] <= o.f
        prob.constraints += o.v[2][t_s]'*F2[:, t] <= o.f

        # friction cone
        prob.constraints += norm((I - o.v[1][t_s]*o.v[1][t_s]')*F1[:, t]) <= o.mu*o.v[1][t_s]'*F1[:, t] # friction cone
        prob.constraints += norm((I - o.v[2][t_s]*o.v[2][t_s]')*F2[:, t]) <= o.mu*o.v[2][t_s]'*F2[:, t] # friction cone

        # dynamics
        u = 1/o.mass * (F1[:, t] + F2[:, t]) + o.g
        prob.constraints += Z[vel_ind, t+1] == Z[vel_ind, t] + u*dt
        prob.constraints += Z[pos_ind, t+1] == Z[pos_ind, t] + Z[vel_ind, t]*dt + u*.5*dt^2
    end

    b = @benchmark Convex.solve!($prob, ECOS.Optimizer(verbose=0))
    Convex.solve!(prob, ECOS.Optimizer(verbose=0))
    println(prob.status)
    u_curr = [F1.value[:, 1]; F2.value[:, 1]]

    return u_curr, b
end
