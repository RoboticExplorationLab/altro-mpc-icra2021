# Create objective for ALTRO

using LinearAlgebra, SparseArrays, StaticArrays
using Convex, ECOS

include("utils.jl")
include("general_socp.jl")
include("struct_setup.jl")


function make_problem_CONVEX_COLD(r::Rocket, obj_opts::ObjectiveOptions,
                    t_opts::TrajectoryOptions, out_opts::OutputOptions;
                    includeGoal::Bool = true, ground_level::Float64 = 0.0)
    println("Chose Convex with Cold Start")

    n, m = r.num_states, r.num_controls
    X = Variable(n, t_opts.N) # State Trajectory
    U = Variable(m, t_opts.N - 1) # State Trajectory

    # First, we build the cost function.
    objective = obj_opts.Qk * sumsquares(X[:,1:t_opts.N-1]) +
                        obj_opts.Qfk * sumsquares(X[:,t_opts.N]) +
                        obj_opts.Rk * sumsquares(U)
    prob = minimize(objective)

    if out_opts.verbose
        println("Objective Set")
    end

    # Next, we put the constraints

    # First up are the dynamics constraints
    constraints = Constraint[ X[:,k+1] == r.Adis * X[:,k] + r.Bdis * U[:,k] +
                                                r.Gdis for k in 1:t_opts.N-1]

    # Including the initial conditions
    push!(constraints, X[:,1] == t_opts.x0)

    if out_opts.verbose
        println("Dynamics Constraint Set")
    end

    # Second up is the goal constraint
    if includeGoal
        push!(constraints, X[:, end] == t_opts.xf)
        if out_opts.verbose
            println("Goal Constraint Set")
        end
    end

    # Third up is the ground constraint
    push!(constraints, X[3, :] >= ground_level)

    if out_opts.verbose
        println("Ground Constraint Set")
    end


    # Fourth up is the max thrust constraint
    [push!(constraints, norm(U[:,i]) <= r.u_max) for i in t_opts.N - 1]

    if out_opts.verbose
        println("Max Thrust Constraint Set")
    end

    # Fifth up is the max thrust angle constraint
    maxTAalpha = getAlpha(r.theta)
    [push!(constraints, norm(U[1:2, i]) <= maxTAalpha * U[3, i])
                                                        for i in t_opts.N - 1]

    if out_opts.verbose
        println("Max Thrust Angle Constraint Set")
    end

    # Now we are done with the constraints
    # Lastly, we set the initial controls to a hover
    set_value!(U, hcat([r1.grav for k = 1:t1.N - 1]...))

    return prob
end



function make_problem_CONVEX_WARM(r::Rocket, obj_opts::ObjectiveOptions,
                    t_opts::TrajectoryOptionsWARM, out_opts::OutputOptions;
                    includeGoal::Bool = true, ground_level::Float64 = 0.0)
    println("Chose Convex with Warm Start")

    # model = LinearModel(r.Acont, r.Bcont, r.Gcont)
    # n, m = size(model)
    #
    # xTrack = t_opts.ref_traj_x[1:t_opts.Horizon]
    # uTrack = t_opts.ref_traj_u[1:(t_opts.Horizon - 1)]
    #
    # n, m = r.num_states, r.num_controls
    # println("n = $n, m = $m")
    #
    # # Set the Objective
    # cost_list = TrajectoryOptimization.CostFunction[]
    #
    # # Set-UP the LQR Tracking Objective Function
    # Q = obj_opts.Qk * Diagonal(@SVector ones(n))
    # Qf = obj_opts.Qfk * Diagonal(@SVector ones(n))
    # R = obj_opts.Rk * Diagonal(@SVector ones(m))
    #
    # for k in 1:(t_opts.Horizon-1)
    #     xf = xTrack[k]
    #     uf = uTrack[k]
    #     push!(cost_list, LQRCost(Q, R, xf, uf))
    # end
    #
    # xf = xTrack[end]
    # qf = -Qf*xf
    # cf = 0.5*xf'*Qf*xf
    # push!(cost_list, QuadraticCost(Qf, zeros(0,0),zeros(0,size(Qf,1)), qf,
    #                     zeros(0), cf, terminal = true))
    #
    # obj = Objective(cost_list)
    #
    # if out_opts.verbose
    #     println("Objective Set")
    # end
    #
    # # Create Empty ConstraintList
    # conSet = ConstraintList(n, m, t_opts.Horizon)
    #
    # # Bound Constraint that the rocket doesn't crash into the ground
    # # This constraint can be made more complicated for difficult terrain
    # # This constraint can also be a glidescope constraint.
    # bnd = BoundConstraint(n,m, x_min=[-Inf, -Inf, ground_level,
    #                                   -Inf, -Inf, -Inf])
    # if out_opts.verbose
    #     println("Ground Constraint Set")
    # end
    #
    # # Norm Constrant that reflects the max thrust the rocket can provide
    # maxT = NormConstraint(n, m, r.u_max,
    #                     TrajectoryOptimization.SecondOrderCone(), :control)
    # add_constraint!(conSet, maxT, 1:(t_opts.Horizon-1))
    # if out_opts.verbose
    #     println("Max Thrust Constraint Set")
    # end
    #
    # # Generalized Norm Constraint that reflects the max thrust angle constraint
    # # Based on the rocket gimbal
    # maxTAalpha = getAlpha(r.theta)
    # thetaARocket = SizedMatrix{3,3}([1.0 0 0; 0 1.0 0; 0 0 0])
    # thetacRocket = SVector{3}([0; 0; maxTAalpha])
    # maxTA = NormConstraint2(n, m, thetaARocket, thetacRocket,
    #                     TrajectoryOptimization.SecondOrderCone(), :control)
    # add_constraint!(conSet, maxTA, 1:(t_opts.Horizon-1))
    # if out_opts.verbose
    #     println("Max Thrust Angle Constraint Set")
    # end
    #
    # # Goal Constraint that the rocket must reach the landing site.
    # goal = GoalConstraint(xTrack[end])
    # add_constraint!(conSet, goal, t_opts.Horizon)
    # if out_opts.verbose
    #     println("Goal Constraint Set at $xf")
    # end
    #
    # # Package the objective and constraints into a "problem" type
    # tf = t_opts.Horizon * t_opts.dt
    # prob = Problem(model, obj, t_opts.xf, tf, x0=xTrack[1], constraints=conSet)
    # if out_opts.verbose
    #     println("Problem Packaged: \nxf = $xf, \nx0 = $(prob.x0), \ntf = $tf")
    # end
    #
    # initial_controls!(prob, uTrack)
    # initial_states!(prob, xTrack)
    # rollout!(prob)
    #
    # return prob
end
