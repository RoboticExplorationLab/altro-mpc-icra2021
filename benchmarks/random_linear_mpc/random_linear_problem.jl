"""
Generate the random linear problem with `n` states, `m` controls, and `N` knot points
"""
function gen_random_linear(n,m,N, dt=0.1)
    # Create model
    A,B = gendiscrete(n,m)
    model = RD.LinearModel(A, B; dt=dt)

    # Objective
    Q = Diagonal(10*rand(n))
    R = Diagonal(0.1*ones(m))
    Qf = Q * (N-1)

    u_bnd = 2  # 2 sigma

    x0 = @SVector zeros(n)
    xf = copy(x0) 
    obj = LQRObjective(Q, R, Qf, xf, N)

    # Constraints
    constraints = ConstraintList(n, m, N)
    bound = BoundConstraint(n, m, u_min=-u_bnd, u_max=u_bnd) 
    add_constraint!(constraints, bound, 1:N-1)
    # add_constraint!(constraints, GoalConstraint(xf), N)

    # Problem
    tf = (N-1)*dt
    prob = Problem(model, obj, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)
    rollout!(prob)
    return prob
end

"""
Generate an OSQP problem from a TrajectoryOptimization Problem
"""
function gen_OSQP(prob0::Problem, opts::SolverOptions)
    prob = copy(prob0)
    TO.add_dynamics_constraints!(prob)
    n,m,N = size(prob)
    nlp = TrajOptNLP(prob, remove_bounds=true)
    NN = N*n + (N-1)*m

    # Cost function
    TO.hess_f!(nlp)
    P = nlp.data.G
    q = zeros(n+m, N)
    for k = 1:N
        q[1:n,k] .= prob.obj[k].q
        q[n+1:end,k] .= prob.obj[k].r
    end
    dt = prob.Z[1].dt
    q[:,1:N-1] .*= dt
    q = q[1:end-m]

    # Constraints
    TO.jac_c!(nlp)
    A = nlp.data.D
    gL, gU = TO.constraint_bounds(nlp)
    zL, zU = TO.primal_bounds!(nlp)
    active = isfinite.(zL) .| isfinite.(zU)  # filter out the bounds at infinity

    # Put dynamics on top of bound constraints
    A = [A; I(NN)[active,:]]
    u = [gU; zU[active]]
    l = [gL; zL[active]]

    model = OSQP.Model()
    OSQP.setup!(model, P=P, q=q, A=A, l=l, u=u;
        verbose=opts.verbose>0,
        eps_abs=opts.cost_tolerance,
        eps_rel=opts.cost_tolerance,
        eps_prim_inf=opts.constraint_tolerance,
        eps_dual_inf=opts.constraint_tolerance,
    )
    return model, l,u
end

"""
Create a Trajectory Optimization problem that tracks the trajectory in `prob`,
using the same constraints, minus the goal constraint. Tracks the first `N`
time steps.
"""
function gen_tracking_problem(prob::Problem, N)
    n,m = size(prob)
    dt = prob.Z[1].dt
    tf = (N-1)*dt

    # Get sub-trajectory
    Z = Traj(prob.Z[1:N])
    x0 = state(Z[1])
    xf = state(Z[N])  # this actually doesn't effect anything

    # Generate a cost that tracks the trajectory
    Q = Diagonal(@SVector fill(10.0, n))
    R = Diagonal(@SVector fill(0.1, m))
    obj = TO.TrackingObjective(Q, R, Z) 

    # Use the same constraints, except the Goal constraint
    cons = ConstraintList(n,m,N)
    for (inds, con) in zip(prob.constraints)
        if !(con isa GoalConstraint)
            if inds.stop > N
                inds = inds.start:N
            end
            add_constraint!(cons, con, inds)
        end
    end

    prob = Problem(prob.model, obj, xf, tf, x0=x0, constraints=cons, 
        integration=TO.integration(prob)
    )
    initial_trajectory!(prob, Z)
    return prob
end


"""
Run the MPC problem to track `Z_track`

Compares the OSQP and ALTRO solutions and timing results
"""
function run_MPC(prob_mpc, opts_mpc, Z_track)
    # Generate ALTRO solver
    altro = ALTROSolver(prob_mpc, opts_mpc)

    # Generate OSQP solver
    osqp,l,u = gen_OSQP(prob_mpc, opts_mpc)

    # Some initialization
    N_mpc = prob_mpc.N
    xinds = [(k-1)*(n+m) .+ (1:n) for k = 1:N_mpc]
    uinds = [(k-1)*(n+m) + n .+ (1:m) for k = 1:N_mpc-1]
    xi = vcat(xinds...)
    ui = vcat(uinds...)
    x0_l = view(l, (N_mpc-1)*n .+ (1:n))  # views into OSQP constraints
    x0_u = view(u, (N_mpc-1)*n .+ (1:n))
    q = zeros(RD.num_vars(prob_mpc.Z))

    dt = prob_mpc.Z[1].dt
    t0 = prob_mpc.t0
    k_mpc = 1
    obj = prob_mpc.obj

    num_iters = length(Z_track) - prob_mpc.N
    err_traj = zeros(num_iters,2)
    err_x0 = zeros(num_iters,2)
    iters = zeros(Int, num_iters,2)
    times = zeros(num_iters,2)

    # Solve initial iteration
    solve!(altro)
    res = OSQP.solve!(osqp)

    # Create views into OSQP dual variables
    y = res.y
    λ = view(y, 1:(N_mpc-1)*n)        # dynamics multipliers
    μ = view(y, N_mpc*n+1:length(y))  # bounds multipliers

    for i = 1:num_iters
        # Update initial time
        t0 += dt
        k_mpc += 1
        TO.set_initial_time!(prob_mpc, t0)

        # Update initial state by using 1st control, and adding some noise 
        x0 = discrete_dynamics(TO.integration(prob), prob_mpc.model, prob_mpc.Z[1])
        x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
        TO.set_initial_state!(prob_mpc, x0)

        # Update tracking cost
        TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

        # Shift the initial trajectory
        RD.shift_fill!(prob_mpc.Z)

        # Shift the multipliers and penalties
        Altro.shift_fill!(TO.get_constraints(altro))

        # Update OSQP
        x0_l .= x0  # initial condition
        x0_u .= x0
        for k = 1:length(obj) - 1  # objective
            q[xinds[k]] .= obj[k].q * dt
            q[uinds[k]] .= obj[k].r * dt
        end
        q[xinds[end]] .= obj[end].q
        OSQP.update!(osqp, q=q, l=l, u=u)
        xf_tmp = res.x[end-n+1:end]
        x = circshift(res.x, -(n+m))  # shift the primals
        x[end-n+1:end] .= xf_tmp
        λ_new = circshift(λ, -n)      # shift the dynamics multipliers
        μ_new = circshift(μ, -m)      # shift the bounds multipliers
        λ .= λ_new
        μ .= μ_new
        OSQP.warm_start!(osqp, x=x, y=y)

        # Solve the updated problem
        solve!(altro)
        OSQP.solve!(osqp, res)

        iters[i,1] = iterations(altro)
        iters[i,2] = res.info.iter
        times[i,1] = altro.stats.tsolve
        times[i,2] = res.info.solve_time * 1000 

        # Compare the solutions
        X_altro = vcat(Vector.(states(altro))...)
        X_osqp = res.x[xi] 
        err_traj[i,1] = norm(X_altro - X_osqp, Inf)

        U_altro = vcat(Vector.(controls(altro))...)
        U_osqp = res.x[ui] 
        err_traj[i,2] = norm(U_altro - U_osqp, Inf)

        err_x0[i,1] = norm(X_altro[1:n] - x0)
        err_x0[i,2] = norm(X_osqp[1:n] - x0)
    end
    return Dict(:time=>times, :iter=>iters, :err_traj=>err_traj, :err_x0=>err_x0) 
end
