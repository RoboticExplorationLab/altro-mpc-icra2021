# MPC in ALTRO

using Altro, TrajectoryOptimization

"""
    update_MPC_controller!(prob::TrajectoryOptimization.Problem,
                k_new, x_curr, obj_opts::ObjectiveOptions,
                t_opts::TrajectoryOptionsWARM, out_opts::OutputOptions)

Updates ALTRO MPC in place. (Specifically, specify k_new and x_curr)
"""
function update_MPC_controller!(prob::TrajectoryOptimization.Problem,
                k_new, x_curr, obj_opts::ObjectiveOptions,
                t_opts::TrajectoryOptionsWARM, out_opts::OutputOptions)

    hor = t_opts.Horizon

    TrajectoryOptimization.set_initial_state!(prob, x_curr)

    X_warm = [[x_curr]; t_opts.ref_traj_x[k_new + 1:k_new + 1 + hor]]

    if out_opts.verbose
        println("Updated Initial State")
    end

    # Update the Reference Trajectory
    for k in 1:(hor - 1)
        TrajectoryOptimization.set_LQR_goal!(prob.obj.cost[k], X_warm[k],
                                            t_opts.ref_traj_u[k + k_new - 1])
    end
    TrajectoryOptimization.set_LQR_goal!(prob.obj.cost[end], X_warm[end])

    if out_opts.verbose
        println("Updated LQR")
    end

    # Update the Goal Constraint
    TrajectoryOptimization.set_goal_state!(prob.constraints.constraints[end],
                                                X_warm[end])

    if out_opts.verbose
        println("Updated Goal")
    end

    initial_controls!(prob, t_opts.ref_traj_u[k_new:k_new + hor])
    initial_states!(prob, X_warm)

    if out_opts.verbose
        println("Warm Start Set")
    end

    rollout!(prob)

    return prob
end


"""
    run_MPC_controller!(prob::TrajectoryOptimization.Problem, opts,
                                out_opts::OutputOptions)

Solves the MPC problem and returns the states and controls.
"""
function run_MPC_controller!(prob::TrajectoryOptimization.Problem, opts,
                                out_opts::OutputOptions)

    altro_mpc = ALTROSolver(prob, opts)
    set_options!(altro_mpc, show_summary=false)
    solve!(altro_mpc)

    vio = max_violation(altro_mpc)

    if out_opts.verbose
        println(vio)
    end

    # Grab the new optimized trajectory
    X_new = states(altro_mpc)
    U_new = controls(altro_mpc)

    return X_new, U_new, vio

end


"""
    loop_MPC_controller!(prob::TrajectoryOptimization.Problem, opts,
                obj_opts::ObjectiveOptions, t_opts::TrajectoryOptionsWARM,
                out_opts::OutputOptions; x_threshold::Float64 = 0.002,
                disturbance = zeros(6, 1))

Runs the full loop for the MPC Controller.
"""
function loop_MPC_controller!(r::Rocket, prob::TrajectoryOptimization.Problem,
                opts, obj_opts::ObjectiveOptions, t_opts::TrajectoryOptionsWARM,
                out_opts::OutputOptions; x_threshold::Float64 = 0.002,
                disturbance = zeros(6, 1), max_vio::Float64 = 1.0)

    hor = t_opts.Horizon
    t_iter = t_opts.dt * hor
    MPC_Iteration_States = []
    MPC_Iteration_Controls = []

    x_curr = t_opts.ref_traj_x[1]

    for iter in 1:(t_opts.N - hor - 1)
        # Step 1: Run the MPC
        if norm(x_curr - t_opts.ref_traj_x[iter]) > x_threshold
            # Solve MPC Trajectory
            X_mpc_k, U_mpc_k, vio = run_MPC_controller!(prob, opts, out_opts)

            if vio > max_vio
                # Constraint is too high
                X_mpc_k = t_opts.ref_traj_x[iter:iter + hor]
                U_mpc_k = t_opts.ref_traj_u[iter:iter + hor - 1]

                if true
                    println("Solve Failed")
                end
            end
        else
            # Use the reference trajectory
            X_mpc_k = t_opts.ref_traj_x[iter:iter + hor]
            U_mpc_k = t_opts.ref_traj_u[iter:iter + hor - 1]

            if out_opts.verbose
                println("Reference Trajectory is Applicable")
            end
        end

        # Step 2: Save the Current Controls
        u_curr = U_mpc_k[1]
        push!(MPC_Iteration_Controls, u_curr)

        # Step 3: Propogate the dynamics
        # updated_model = full_model(rocket, x_curr, dt, [u_curr], wind,
        #                             num_timesteps = num_tsteps, t0 = t_curr)
        x_curr = vec(propogate_dynamics(r, x_curr, u_curr, disturbance))

        # Step 4: Save the Current State
        push!(MPC_Iteration_States, x_curr)

        # Step 5: Update the MPC for the next iteration
        update_MPC_controller!(prob, iter, x_curr, obj_opts, t_opts, out_opts)

        if true
            print("At Iteration $iter @ ")
            println(norm(x_curr - t_opts.ref_traj_x[iter + 1]))
        end
    end

    return MPC_Iteration_States, MPC_Iteration_Controls
end
