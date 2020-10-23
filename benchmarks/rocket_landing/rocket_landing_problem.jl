include("src/general_socp.jl")
include(joinpath(@__DIR__,"..","mpc.jl"))

"""
Build the continuous time linear rocket model, given the mass.
"""
function RocketModel(mass, grav, dt; integration=RD.Exponential)
    A = [
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3);
    ]
    B = [
        zeros(3,3);
        -1/mass * I(3)
    ]
    d = [
        zeros(3);
        grav
    ]

    # Continuous Model
    cmodel =  LinearModel(A,B,d)
    if dt == zero(dt)
        return cmodel
    end

    # Discrete Model
    model = LinearizedModel(cmodel, dt=dt, is_affine=true, integration=integration)
end

getalpha(θ, deg=true) = deg ? tan(θ) : tand(θ)

function RocketProblem(N=101, tf=10.0;
        x0 = SA_F64[4, 2, 20, -3, 2, -5],
        Qk = 1e-2, 
        Qfk = 100, 
        Rk = 1e-1,
        ground_level = 0.0,
        mass = 10.0,
        perWeightMax = 2.0,
        θ_max = 7.0,  # deg
        include_goal = true,
        integration=RD.Exponential
    )
    """
    Model
    """
    n,m = 6,3
    dt = tf / (N-1)
    grav = SA[0,0,9.81]
    model = RocketModel(mass, grav, dt, integration=integration)

    # Initial and Final Condition
    xf = @SVector zeros(6)

    """
    Objective
    """
    Q = Diagonal(@SVector fill(Float64(Qk), n))
    R = Diagonal(@SVector fill(Float64(Rk), m))
    Qf = Diagonal(@SVector fill(Float64(Qfk), n))
    obj = LQRObjective(Q,R,Qf,xf,N)

    """
    Constraints
    """
    cons = ConstraintList(n,m,N)

    # Goal
    include_goal && add_constraint!(cons, GoalConstraint(xf), N)  # goal constraint

    # Floor
    x_min = SA[-Inf, -Inf, ground_level, -Inf, -Inf, -Inf]
    add_constraint!(cons, BoundConstraint(n,m, x_min = x_min), 2:N-1)  # floor constraint

    # Max Thrust
    u_bnd = mass * abs(grav[3]) * perWeightMax
    maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
    add_constraint!(cons, maxThrust, 1:N-1)

    # Max Thrust Angle
    α_max = getalpha(θ_max)
    ARocket = SA_F64[
        1 0 0; 
        0 1 0; 
        0 0 0
    ]
    cRocket = SA_F64[0, 0, α_max]
    maxAngle = NormConstraint2(n, m, ARocket, cRocket, TO.SecondOrderCone(), :control)
    add_constraint!(cons, maxAngle, 1:N-1)

    """
    Problem
    """
    prob = Problem(model, obj, xf, tf, x0=x0, constraints=cons, integration=RD.PassThrough)

    # Initialization
    U0 = [grav for k = 1:N-1]
    initial_controls!(prob, U0)
    rollout!(prob)

    return prob
end

function RocketProblemMPC(prob::Problem, N; kwargs...)
    prob = gen_tracking_problem(prob, N; kwargs...) 
end

function run_Rocket_MPC(prob_mpc, opts_mpc, Z_track, num_iters = length(Z_track) - prob_mpc.N)
    altro = ALTROSolver(prob_mpc, opts_mpc)
    
    # ecos = gen_ECOS(prob_mpc, opts_mpc)

    # Solve initial iteration
    solve!(altro)

    err_traj = zeros(num_iters,2)
    err_x0 = zeros(num_iters,2)
    iters = zeros(Int, num_iters,2)
    times = zeros(num_iters,2)

    t0 = 0
    k_mpc = 1
    x0 = SVector(prob_mpc.x0)
    X_traj = [copy(x0) for k = 1:num_iters+1]
    for i = 1:num_iters
        # Update initial time
        t0 += dt
        k_mpc += 1
        TO.set_initial_time!(prob_mpc, t0)

        # Update initial state by using 1st control, and adding some noise 
        x0 = discrete_dynamics(TO.integration(prob_mpc), prob_mpc.model, prob_mpc.Z[1])
        x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
        TO.set_initial_state!(prob_mpc, x0)
        X_traj[i+1] = x0

        # Update tracking cost
        TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

        # Shift the initial trajectory
        RD.shift_fill!(prob_mpc.Z)

        # Shift the multipliers and penalties
        Altro.shift_fill!(TO.get_constraints(altro))

        # Solve the updated problem
        solve!(altro)
        # solve ECOS...

        iters[i,1] = iterations(altro)
        # iters[i,2] = res.info.iter
        # times[i,1] = median(b).time * 1e-6 
        times[i,1] = altro.stats.tsolve
        # times[i,2] = res.info.solve_time * 1000 

        # compare the solutions....
    end
    return X_traj, Dict(:time=>times, :iter=>iters, :err_traj=>err_traj, :err_x0=>err_x0) 
end