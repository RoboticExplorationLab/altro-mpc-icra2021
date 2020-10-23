# Create objective for ALTRO

using LinearAlgebra, SparseArrays, StaticArrays
using RobotDynamics, TrajectoryOptimization, Altro

include("utils.jl")
include("general_socp.jl")

function make_problem_ALTRO_COLD(r::Rocket, obj_opts::ObjectiveOptions,
                    t_opts::TrajectoryOptions, out_opts::OutputOptions;
                    includeGoal::Bool = true, ground_level::Float64 = 0.0)
    println("Chose ALTRO with Cold Start")

    model = LinearModel(r.Acont, r.Bcont, r.Gcont)
    n, m = size(model)

    # Set-UP the LQR Objective Function
    Q = obj_opts.Qk * Diagonal(@SVector ones(n))
    Qf = obj_opts.Qfk * Diagonal(@SVector ones(n))
    R = obj_opts.Rk * Diagonal(@SVector ones(m))
    obj = LQRObjective(Q,R,Qf,t_opts.xf,t_opts.N)
    if out_opts.verbose
        println("Objective Set")
    end

    # Create Empty ConstraintList
    conSet = ConstraintList(n,m,t_opts.N)

    if includeGoal
        # Goal Constraint that the rocket must reach the landing site.
        goal = GoalConstraint(t_opts.xf)
        TrajectoryOptimization.add_constraint!(conSet, goal, t_opts.N)
        if out_opts.verbose
            println("Goal Constraint Set")
        end
    end

    # Bound Constraint that the rocket doesn't crash into the ground
    # This constraint can be made more complicated for difficult terrain
    # This constraint can also be a glidescope constraint.
    bnd = BoundConstraint(n,m, x_min=[-Inf, -Inf, ground_level,
                                      -Inf, -Inf, -Inf])
    TrajectoryOptimization.add_constraint!(conSet, bnd, 1:t_opts.N-1)
    if out_opts.verbose
        println("Ground Constraint Set")
    end

    # Norm Constrant that reflects the max thrust the rocket can provide
    maxT = NormConstraint(n, m, r.u_max,
                        TrajectoryOptimization.SecondOrderCone(), :control)
    TrajectoryOptimization.add_constraint!(conSet, maxT, 1:t_opts.N-1)
    if out_opts.verbose
        println("Max Thrust Constraint Set")
    end

    # Generalized Norm Constraint that reflects the max thrust angle constraint
    # Based on the rocket gimbal
    maxTAalpha = getAlpha(r.theta)
    thetaARocket = SizedMatrix{3,3}([1.0 0 0; 0 1.0 0; 0 0 0])
    thetacRocket = SVector{3}([0; 0; maxTAalpha])
    maxTA = NormConstraint2(n, m, thetaARocket, thetacRocket,
                        TrajectoryOptimization.SecondOrderCone(), :control)
    TrajectoryOptimization.add_constraint!(conSet, maxTA, 1:t_opts.N-1)
    if out_opts.verbose
        println("Max Thrust Angle Constraint Set")
    end

    # Package the objective and constraints into a "problem" type
    tf = (t_opts.N-1) * t_opts.dt
    prob = Problem(model, obj, t_opts.xf, tf, x0=t_opts.x0, constraints=conSet)
    if out_opts.verbose
        println("Problem Made with tf = $tf")
    end

    # Set the initial controls to a hover
    u0 = r.grav # controls that would nominally hover
    U0 = [u0 for k = 1:t_opts.N-1] # vector of the small controls
    initial_controls!(prob, U0)
    rollout!(prob)
    if out_opts.verbose
        println("Set Initial Controls")
    end

    return prob
end



function make_problem_ALTRO_WARM(r::Rocket, obj_opts::ObjectiveOptions,
                    t_opts::TrajectoryOptionsWARM, out_opts::OutputOptions;
                    includeGoal::Bool = true, ground_level::Float64 = 0.0)

    model = LinearModel(r.Acont, r.Bcont, r.Gcont)
    n, m = size(model)

    xTrack = t_opts.ref_traj_x[1:t_opts.Horizon]
    uTrack = t_opts.ref_traj_u[1:(t_opts.Horizon - 1)]

    n, m = r.num_states, r.num_controls
    println("n = $n, m = $m")

    # Set the Objective
    cost_list = TrajectoryOptimization.CostFunction[]

    # Set-UP the LQR Tracking Objective Function
    Q = obj_opts.Qk * Diagonal(@SVector ones(n))
    Qf = obj_opts.Qfk * Diagonal(@SVector ones(n))
    R = obj_opts.Rk * Diagonal(@SVector ones(m))

    for k in 1:(t_opts.Horizon-1)
        xf = xTrack[k]
        uf = uTrack[k]
        push!(cost_list, LQRCost(Q, R, xf, uf))
    end

    xf = xTrack[end]
    qf = -Qf*xf
    cf = 0.5*xf'*Qf*xf
    push!(cost_list, QuadraticCost(Qf, zeros(0,0),zeros(0,size(Qf,1)), qf,
                        zeros(0), cf, terminal = true))

    obj = Objective(cost_list)

    if out_opts.verbose
        println("Objective Set")
    end

    # Create Empty ConstraintList
    conSet = ConstraintList(n, m, t_opts.Horizon)

    # Bound Constraint that the rocket doesn't crash into the ground
    # This constraint can be made more complicated for difficult terrain
    # This constraint can also be a glidescope constraint.
    bnd = BoundConstraint(n,m, x_min=[-Inf, -Inf, ground_level,
                                      -Inf, -Inf, -Inf])
    TrajectoryOptimization.add_constraint!(conSet, bnd, 1:(t_opts.Horizon-1))
    if out_opts.verbose
        println("Ground Constraint Set")
    end

    # Norm Constrant that reflects the max thrust the rocket can provide
    maxT = NormConstraint(n, m, r.u_max,
                        TrajectoryOptimization.SecondOrderCone(), :control)
    TrajectoryOptimization.add_constraint!(conSet, maxT, 1:(t_opts.Horizon-1))
    if out_opts.verbose
        println("Max Thrust Constraint Set")
    end

    # Generalized Norm Constraint that reflects the max thrust angle constraint
    # Based on the rocket gimbal
    maxTAalpha = getAlpha(r.theta)
    thetaARocket = SizedMatrix{3,3}([1.0 0 0; 0 1.0 0; 0 0 0])
    thetacRocket = SVector{3}([0; 0; maxTAalpha])
    maxTA = NormConstraint2(n, m, thetaARocket, thetacRocket,
                        TrajectoryOptimization.SecondOrderCone(), :control)
    TrajectoryOptimization.add_constraint!(conSet, maxTA, 1:(t_opts.Horizon-1))
    if out_opts.verbose
        println("Max Thrust Angle Constraint Set")
    end

    # Goal Constraint that the rocket must reach the landing site.
    goal = GoalConstraint(xTrack[end])
    TrajectoryOptimization.add_constraint!(conSet, goal, t_opts.Horizon)
    if out_opts.verbose
        println("Goal Constraint Set at $xf")
    end

    # Package the objective and constraints into a "problem" type
    tf = t_opts.Horizon * t_opts.dt
    prob = Problem(model, obj, t_opts.xf, tf, x0=xTrack[1], constraints=conSet)
    if out_opts.verbose
        println("Problem Packaged: \nxf = $xf, \nx0 = $(prob.x0), \ntf = $tf")
    end

    initial_controls!(prob, uTrack)
    initial_states!(prob, xTrack)
    rollout!(prob)

    return prob
end
