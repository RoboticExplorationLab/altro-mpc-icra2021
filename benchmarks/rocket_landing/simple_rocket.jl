function gen_JuMP_rocket(prob::Problem, opts::SolverOptions, optimizer;
        goal_constraint = false,
        warmstart = true
    )
    n,m,N = size(prob)
    NN = N*n + (N-1)*m
    dt = prob.Z[1].dt

    inds = reshape(1:(n+m)*N,n+m,N)
    xinds = [z[1:n] for z in eachcol(inds)]
    uinds = [z[n+1:end] for z in eachcol(inds)][1:N-1]

    # Create the model 
    model = Model(optimizer)
    @variable(model, z[1:NN]) 
    if warmstart 
        z0 = vcat([Vector(RD.get_z(z)) for z in get_trajectory(prob)]...)
        JuMP.set_start_value.(z, z0)    
    end

    # Cost function
    P,q,c = get_batch_objective(prob)
    @objective(model, Min, 0.5*dot(z,P,z) + dot(q,z) + c)

    # Dynamics Constraints
    A = Array(RD.get_A(prob.model.linmodel,1))
    B = Array(RD.get_B(prob.model.linmodel,1))
    d = Array(RD.get_d(prob.model.linmodel,1))
    for k = 1:N-1
        @constraint(model, A*z[xinds[k]] .+ B*z[uinds[k]] .+ d .== z[xinds[k+1]])
    end

    # Initial condition 
    @constraint(model, z[xinds[1]] .== prob.x0)

    # Goal constraint
    if goal_constraint 
        @constraint(model, z[xinds[N]] .== prob.xf)
    end

    # Thrust cone constraint
    con = prob_mpc.constraints[findfirst(x-> x isa NormConstraint, prob_mpc.constraints)]
    maxThrust = con.val
    for k = 1:N-1
        @constraint(model, [maxThrust, z[uinds[k]]...] in JuMP.SecondOrderCone()) 
    end

    # Thrust angle constraint
    cones = prob_mpc.constraints[findall(x-> x isa NormConstraint2, prob_mpc.constraints)]
    α_max = cones[1].c[3]
    for k = 1:N-1
        u1,u2,u3 = z[uinds[k]]
        @constraint(model, [α_max * u3, u1, u2] in JuMP.SecondOrderCone())
    end

    return model, z, (P,q,c)
end

function mpc_update(altro, prob_mpc, Z_track, t0, k_mpc)
    TO.set_initial_time!(prob_mpc, t0)

    # Propagate the system forward w/ noise
    x0 = discrete_dynamics(TO.integration(prob_mpc),
                                prob_mpc.model, prob_mpc.Z[1])
    pos_mag = norm(x0[1:3])
    vel_mag = norm(x0[4:6])
    noise = [
        (@SVector randn(3)) * pos_mag / 1000; 
        (@SVector randn(3)) * vel_mag / 100
    ]
    x0 += noise
    TO.set_initial_state!(prob_mpc, x0)

    # Update tracking cost
    TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

    # Shift the initial trajectory
    RD.shift_fill!(prob_mpc.Z)

    # Shift the multipliers and penalties
    Altro.shift_fill!(TO.get_constraints(altro.solver_al))
end

function get_batch_objective(prob)
    prob_ = copy(prob)
    TO.add_dynamics_constraints!(prob_)
    n,m,N = size(prob_)
    nlp = TrajOptNLP(prob_, remove_bounds=true)

    TO.hess_f!(nlp)
    P = nlp.data.G
    q = zeros(n+m, N)
    for k = 1:N
        q[1:n,k] .= prob_.obj[k].q
        q[n+1:end,k] .= prob_.obj[k].r
    end
    dt = prob_.Z[1].dt
    q[:,1:N-1] .*= dt
    q = q[1:end-m]
    dts = fill(prob_.Z[1].dt,N)
    dts[end] = 1
    c = [c.c for c in prob_.obj]'dts
    return P,q,c
end

function run_Rocket_MPC(prob_mpc, opts_mpc, Z_track; 
        num_iters = length(Z_track) - prob_mpc.N,
        ecos_tol = opts_mpc.constraint_tolerance,
        optimizer = JuMP.optimizer_with_attributes(ECOS.Optimizer, 
            "verbose"=>false,
            "feastol"=>ecos_tol,
            "abstol"=>opts_mpc.cost_tolerance,
            "reltol"=>opts_mpc.cost_tolerance
        ),
        benchmark=false
    )
    # Setup and solve the first iteration
    altro = ALTROSolver(prob_mpc, opts_mpc)
    solve!(altro)

    # Initialize results
    X_traj = [zero(state(Z_track[1])) for i = 1:num_iters+1]
    X_traj[1] = state(Z_track[1])
    err_traj = zeros(num_iters,3)  # columns: [state error, control error, ECOS dynamics error]
    err_x0 = zeros(num_iters,2)
    iters = zeros(Int, num_iters,2)
    times = zeros(num_iters,2)
    costs = zeros(num_iters,3)
    X = zero.(states(prob_mpc))
    U = zero.(controls(prob_mpc))

    n,m,N = size(prob_mpc)
    dt = Z_track[1].dt
    inds = reshape(1:(n+m)*N,n+m,N)
    xinds = [z[1:n] for z in eachcol(inds)]
    uinds = [z[n+1:end] for z in eachcol(inds)][1:N-1]
    dts = fill(dt,N)
    dts[end] = 0

    t0 = 0
    k_mpc = 1
    for i = 1:num_iters
        t0 += dt
        k_mpc += 1

        # Update the ALTRO solution, advancing forward by 1 time step
        mpc_update(altro, prob_mpc, Z_track, t0, k_mpc)
        X_traj[i+1] = prob_mpc.x0

        # Generate a JuMP Model for solving the problem with ECOS
        model,z,obj0 = gen_JuMP_rocket(prob_mpc, altro.opts, optimizer, warmstart=false)

        # Solve both problems
        if benchmark
            b = benchmark_solve!(altro, samples=2, evals=2)
            altro.stats.tsolve = median(b).time / 1e6
        else
           solve!(altro)
        end
        optimize!(model)
        
        times[i,1] = altro.stats.tsolve
        times[i,2] = JuMP.solve_time(model) * 1000
        iters[i,1] = iterations(altro)

        # Compare trajectories
        x_ecos = value.(z)
        for k = 1:N
            X[k] = x_ecos[xinds[k]]
            k < N && (U[k] = x_ecos[uinds[k]])
        end
        err_X = maximum(norm.(X - states(altro), Inf))
        err_U = maximum(norm.(U - controls(altro), Inf))
        err_dyn = dynamics_violation(prob_mpc, X, U)
        err_traj[i,:] = [err_X, err_U, err_dyn]
    
        # Compare cost
        P,q,c = get_batch_objective(prob_mpc)
        Z = Traj(X, U, dts)
        x_altro = vcat([Vector(RD.get_z(z)) for z in get_trajectory(altro)]...)
        J_altro = cost(prob_mpc)
        J_ecos = objective_value(model)
        J_ecos2 = cost(prob_mpc.obj, Z) 
        costs[i,1] = J_altro
        costs[i,2] = J_ecos2  # equivalent cost
        costs[i,3] = J_ecos   # includes slack variables?
    end
    X_traj, Dict(:time=>times, :iter=>iters, :cost=>costs, :err_traj=>err_traj), X, U
end

function dynamics_violation(prob,X,U)
    err = [zero(X[1]) for u in U]
    for k = 1:length(U)
        t = prob.Z[k].t
        dt = prob.Z[k].dt
        err = discrete_dynamics(TO.integration(prob), prob.model, X[k], U[k], t, dt) - X[k+1]
    end
    return maximum(norm.(err,Inf))
end


