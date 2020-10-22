
function gen_random_linear(n,m,N, dt=0.1)
    # Create model
    A,B = gendiscrete(n,m)
    model = RD.LinearModel(A, B; dt=dt)

    # Objective
    Q = Diagonal(10*rand(n))
    R = Diagonal(0.1*ones(m))
    Qf = Q * (N-1)

    x̄ = rand(n) .+ 1
    ū = rand(m) * 10 / (N-1)
    # ū = rand(m) * 0.1 
    x0 = (rand(n) .- 1) .* x̄ * 0.5
    xf = zeros(n)
    obj = LQRObjective(Q, R, Qf, xf, N)

    # Constraints
    constraints = ConstraintList(n, m, N)
    bound = BoundConstraint(n, m, x_min=-x̄, x_max=x̄, u_min=-ū, u_max=ū)
    add_constraint!(constraints, bound, 1:N-1)
    add_constraint!(constraints, GoalConstraint(xf), N)

    # Problem
    tf = (N-1)*dt
    prob = Problem(model, obj, xf, tf, x0=x0, constraints=constraints, integration=RD.PassThrough)
    rollout!(prob)
    return prob
end

function gen_OSQP_JuMP(prob::Problem)
    n,m,N = size(prob)
    A,B = Matrix(RD.get_A(prob.model)), Matrix(RD.get_B(prob.model))
    x0 = prob.x0
    xf = prob.xf
    x̄ = prob.constraints[1].z_max[1:n]
    ū = prob.constraints[1].z_max[n+1:n+m]
    dt = prob.Z[1].dt
    Q = prob.obj[1].Q * dt
    R = prob.obj[1].R * dt
    Qf = prob.obj[end].Q

    select(i, n) = (n*(i-1)+1):(n*(i-1)+n)

    jump_model = Model(OSQP.Optimizer)
    set_silent(jump_model)

    @variable(jump_model, x[1:N*n])
    @variable(jump_model, u[1:(N-1)*m])

    objective_exp = @expression(jump_model, 0.5*transpose(x[select(N, n)]) * Qf * x[select(N, n)])

    @constraint(jump_model, x[select(1, n)] .== x0)

    for i=1:N-1
        # dynamics constraints
        @constraint(jump_model, A*x[select(i, n)] + B*u[select(i, m)] .== x[select(i+1, n)])

        # stagewise state cost
        add_to_expression!(objective_exp, 0.5*transpose(x[select(i, n)]) * Q * x[select(i, n)])

        # stagewise control cost
        add_to_expression!(objective_exp, 0.5*transpose(u[select(i, m)]) * R * u[select(i, m)])

        # control/state bound constraints
        @constraint(jump_model, x[select(i, n)] .<= x̄)
        @constraint(jump_model, x[select(i, n)] .>= -x̄)
        @constraint(jump_model, u[select(i, m)] .<= ū)
        @constraint(jump_model, u[select(i, m)] .>= -ū)
    end

    @objective(jump_model, Min, objective_exp)
    return jump_model
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

    # Put dynamics on top of bound constraints
    A = [A; I(NN)]
    u = [gU; zU]
    l = [gL; zL]

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
Generate a set of random initial conditions
"""
function gen_ICs(prob, iters=10)
    n,m,N = size(prob)
    x̄ = prob.constraints[1].z_max[1:n]
    [(rand(n) .- 1) .* x̄ * 0.5 for i = 1:iters]
end

function run_MPC(prob_mpc, opts_mpc, Z_track)
    # Generate ALTRO solver
    altro = ALTROSolver(prob_mpc, opts_mpc)

    # Generate OSQP solver
    osqp,l,u = gen_OSQP(prob_mpc, opts_mpc)

    # Some initialization
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
    num_iters = length(Z_track) - prob_mpc.N
    err_traj = zeros(num_iters,2)
    err_x0 = zeros(num_iters,2)
    obj = prob_mpc.obj

    # Solve initial iteration
    solve!(altro)
    res = OSQP.solve!(osqp)

    X_altro = vcat(Vector.(states(altro))...)
    X_osqp = res.x[xi] 
    err_traj[1,1] = norm(X_altro - X_osqp, Inf)

    U_altro = vcat(Vector.(controls(altro))...)
    U_osqp = res.x[ui] 
    err_traj[1,2] = norm(U_altro - U_osqp, Inf)

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
        

        # Solve the updated problem
        solve!(altro)
        res = OSQP.solve!(osqp)

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
    return err_traj, err_x0 
end

function MPC_OSQP(model, l, u, x0_l, x0_u, ICs)
    times = Float64[]
    for ic in ICs
        t = @elapsed begin
            x0_l .= ic 
            x0_u .= ic 
            OSQP.update!(model, u=u, l=l)
            results = OSQP.solve!(model)
        end
        push!(times, t)
    end
    return times
end

"""
Compare ALTRO and OSQP on a randomly-generated linear problem of the given size.

Run an MPC controller that generates `steps` initial conditions and solves the problem,
warm-starting from the previous initial condition. Reports the median times of the MPC 
iterations.
"""
function run_comparison(n,m,N,steps=100; opts=SolverOptions())
    # Generate the problem
    prob = gen_random_linear(n,m,N)

    # Convert to OSQP
    osqp,l,u = gen_OSQP(prob, opts) 
    NN = N*n + (N-1)*m
    xi = vcat([(k-1)*(n+m) .+ (1:n) for k = 1:N]...)
    ui = vcat([(k-1)*(n+m) + n .+ (1:m) for k = 1:N-1]...)
    x0_l = view(l, (N-1)*n .+ (1:n))
    x0_u = view(u, (N-1)*n .+ (1:n))

    # Solve the first time
    altro = ALTROSolver(prob, opts)
    solve!(altro)
    res = OSQP.solve!(osqp)

    # Compare the results
    # println("Difference in Cost: ", abs(res.info.obj_val - cost(altro)))

    X_altro = vcat(Vector.(states(altro))...)
    X_osqp = res.x[xi] 
    # println("Difference in states: ", norm(X_altro - X_osqp, Inf))

    U_altro = vcat(Vector.(controls(altro))...)
    U_osqp = res.x[ui] 
    # println("Difference in controls: ", norm(U_altro - U_osqp, Inf))

    # Generate the initial conditions 
    ICs = gen_ICs(prob, steps)

    # Change the ALTRO solver options for MPC
    set_options!(altro, reset_duals=false, penalty_initial=1e-1, penalty_scaling=1000., 
        show_summary=false, verbose=0)

    # Run MPC
    t_altro = MPC_Altro(altro, ICs)
    t_osqp = MPC_OSQP(osqp, l, u, x0_l, x0_u, ICs)

    println("Median ALTRO time: ", median(t_altro))
    println("Median OSQP time:  ", median(t_osqp))
    return median(t_altro), median(t_osqp)
end

function comp_plot(xs, times_altro, times_osqp; kwargs...)
    times_altro *= 1000
    times_osqp *= 1000
    avg_altro = mean.(eachrow(times_altro))
    std_altro = std.(eachrow(times_altro))
    avg_osqp = mean.(eachrow(times_osqp))
    std_osqp = std.(eachrow(times_osqp))
    p = plot(ylabel="time (ms)"; kwargs...) 
    plot!(xs, avg_altro, yerr=std_altro, markerstrokecolor=:auto, label="ALTRO")
    plot!(xs, avg_osqp, yerr=std_osqp, markerstrokecolor=:auto, label="OSQP")
    return p
end