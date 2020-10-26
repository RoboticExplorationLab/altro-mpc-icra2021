include("src/general_socp.jl")
include(joinpath(@__DIR__,"..","mpc.jl"))
include("src/utils.jl")

using Altro
using TrajectoryOptimization
const TO = TrajectoryOptimization
using StaticArrays

"""
Build the continuous time linear rocket model, given the mass and angular
velocity of the planet.

For example, ωMars = [0; 6.62e-5; 2.53e-5]
"""
function RocketModel(mass, grav, dt, ωPlanet = [0.0; 0.0; 0.0];
                                                integration=RD.Exponential)
    A = [
        zeros(3,3)      I(3);
        -S(ωPlanet)^2   -2*S(ωPlanet);
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

getalpha(θ, deg=true) = deg ? tand(θ) : tan(θ)

function RocketProblem(N=101, tf=10.0;
        x0 = SA_F64[4, 2, 20, -3, 2, -5],
        Qk = 1e-2,
        Qfk = 100,
        Rk = 1e-1,
        ground_level = 0.0,
        gravity = [0; 0; 9.81],
        mass = 10.0,
        ωPlanet = [0.0; 0.0; 0.0],
        perWeightMax = 2.0,
        θ_thrust_max = 7.0,  # deg
        θ_glideslope = 60.0, # deg
        glide_recover_k = 8,
        include_goal = true,
        integration=RD.Exponential
    )
    """
    Model
    """
    # n = 6 for 3 position entries and 3 velocity entries
    # m = 3 for a 3 thrust entries
    n,m = 6,3
    dt = tf / (N-1)
    grav = gravity # @SVector gravity
    model = RocketModel(mass, grav, dt, ωPlanet, integration=integration)

    # Initial and Final Condition
    xf = @SVector zeros(6)

    """
    Objective
    """
    # Fundamentally assuming that xf is at the origin.
    Q = Diagonal(@SVector fill(Float64(Qk), n))
    R = Diagonal(@SVector fill(Float64(Rk), m))
    Qf = Diagonal(@SVector fill(Float64(Qfk), n))
    obj = LQRObjective(Q,R,Qf,xf,N)

    """
    Constraints
    """
    cons = ConstraintList(n,m,N)

    # Goal Constraint
    #=
    Generally we want the goal as a hard constraint for a precision landing
    (such as on a barge or launch pad). In the case that the rocket is landing
    on a more forgiving area, like a desert, we can remove the goal constraint.
    =#
    include_goal && TO.add_constraint!(cons, GoalConstraint(xf), N)

    # Floor Constraint
    #=
    This prevents trajectories that go in the ground. This is sufficient for
    flat locals, such as oceans, but not more complex locations.
    =#
    # x_min = SA[-Inf, -Inf, ground_level, -Inf, -Inf, -Inf]
    # TO.add_constraint!(cons, BoundConstraint(n,m, x_min = x_min), 2:N-1)

    # Max Thrust
    #=
    This accounts for the limits of the thruster magnitude. As a maximum thrust,
    it is convex. A naive minimum thrust would be non-convex.

    For the sake of inuition, it is easier to understand the max thrust in
    terms of acceleration. Hence we can write it as (mg) * k, where k is a
    per weight factor.

    To have a safe, soft landing, we need k > 1.

    Lastly, the constraint is written as ||u|| ≤ u_max
    =#
    if true
        @assert perWeightMax > 1
        u_bnd = mass * abs(grav[3]) * perWeightMax
        maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
        TO.add_constraint!(cons, maxThrust, 1:N-1)
    end

    # Max Thrust Angle
    #=
    This accounts for the gimbal limits of the thruster.

    We write this as || [u_x; u_y] || ≤ tan(θ) u_z
    To write this in the SOCP form of ||Ax|| ≤ c'x, we take the A and c below.
    =#
    if true
        α_max = getalpha(θ_thrust_max)
        ARocket = SA_F64[
            1 0 0;
            0 1 0;
            0 0 0
        ]
        cRocket = SA_F64[0, 0, α_max]
        maxAngle = NormConstraint2(n, m, ARocket, cRocket,
                                        TO.SecondOrderCone(), :control)
        TO.add_constraint!(cons, maxAngle, 1:N-1)
    end

    # Glideslope Constraint
    #=
    This prevents trajectories that drifts near the surface (which is a safety
    risk). As we make the that angle larger, the constraint is more lenient.
    The angle cannot be larger than 90 deg.
    =#

    if true
        α_glide = getalpha(θ_glideslope)
        AGlide = SA_F64[
            1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0
        ]
        cGlide = SA_F64[0, 0, α_glide, 0, 0, 0]
        glideslope = NormConstraint2(n, m, AGlide, cGlide,
                                        TO.SecondOrderCone(), :state)
        TO.add_constraint!(cons, glideslope, glide_recover_k:N-1)
    end
    # α_glide = cosd(θ_glideslope)
    # AGlide = SA_F64[
    #     1 0 0;
    #     0 1 0;
    #     0 0 1
    # ]
    # cGlide = SA_F64[0, 0, 1/α_glide]
    # glideslope = NormConstraint2(n, m, AGlide, cGlide,
    #                                 TO.SecondOrderCone(), :control)
    # TO.add_constraint!(cons, glideslope, 1:N-1)

    """
    Problem
    """
    prob = TO.Problem(model, obj, xf, tf, x0=x0, constraints=cons,
                                            integration=RD.PassThrough)

    # Initialization
    #=
    We can immediately warm start the solution using a hover, which is the
    exact force to balance gravity (mg).
    =#
    U0 = [mass * grav for k = 1:N-1]
    initial_controls!(prob, U0)
    rollout!(prob)

    return prob
end

function RocketProblemMPC(prob::TO.Problem, N; kwargs...)
    prob = gen_tracking_problem(prob, N; kwargs...)
end

function sanity_check_cost(X, U, X_track, U_track,
                                Qk = 1,
                                Qfk = 1,
                                Rk = 1;
                                verbose = false)
    # (X-X_t)'Q(X-X_t) + (Xf-Xf_t)'Qf(Xf-Xf_t) + (U-U_t)'R(U-U_t)
    N = size(X, 2)

    xCost = Qk * sum([(X[:,i] - X_track[:,i])' * (X[:,i] - X_track[:,i])
                                                    for i in 1:N - 1])
    xfCost = Qfk * (X[:,end] - X_track[:,N])' * (X[:,end] - X_track[:,N])
    uCost = Rk * sum([(U[:,i] - U_track[:,i])' * (U[:,i] - U_track[:,i])
                                                    for i in 1:N - 1])

    verbose && println("xCost = $xCost, xfCost = $xfCost, uCost = $uCost")
    return xCost + xfCost + uCost
end

function get_circle_proj(x, y, rad)
    @assert rad > 0

    mag = norm([x; y])
    if mag < rad
        return x, y
    end

    return (rad / mag) .* (x, y)
end


function run_Rocket_MPC(prob_mpc, opts_mpc, Z_track,
                            num_iters = length(Z_track) - prob_mpc.N)
    # First, Let's generate the ALTRO problem
    altro = ALTROSolver(prob_mpc, opts_mpc)

    # To match the altro problem, we can derive it from the altro
    ecos, X_ecos, U_ecos = gen_ECOS_Rocket(prob_mpc, Z_track, 1, verbose = true)

    reformat_X_track = getX_toECOS(Z_track)
    reformat_U_track = getU_toECOS(Z_track)

    Qk = prob_mpc.obj[1].Q[1]
    Rk = prob_mpc.obj[1].R[1]
    Qfk = prob_mpc.obj[end].Q[1]

    # Solve initial iteration
    Altro.solve!(altro)

    err_traj = zeros(num_iters,2)
    err_x0 = zeros(num_iters,2)
    iters = zeros(Int, num_iters,2)
    times = zeros(num_iters,2)

    # Setup ECOS optimizer
    ecos_optimizer = ECOS.Optimizer(
        verbose=opts_mpc.verbose > 0,
        feastol=opts_mpc.constraint_tolerance,
        abstol=opts_mpc.cost_tolerance,
        reltol=opts_mpc.cost_tolerance,
        max_iter=1e6
    )

    n,m = size(prob_mpc)

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
        x0 = discrete_dynamics(TO.integration(prob_mpc),
                                    prob_mpc.model, prob_mpc.Z[1])
        # x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
        # noise = @SVector [d for d in [randn(3); zeros(3)]]
        # x0 += noise * norm(x0,Inf) / 100  # 1% noise
        noise_pos = @SVector randn(3)
        noise_vel = @SVector randn(3)

        pos_norm = norm(x0[1:3], Inf) / 100 # 1% noise
        vel_norm = norm(x0[4:6], Inf) / 1e6 # 1ppm noise

        x0 += [noise_pos * pos_norm; noise_vel * vel_norm]

        # rad_lim = tand(45) * x0[3]
        # xy_new = 0.9999 .* get_circle_proj(x0[1], x0[2], rad_lim)
        #
        # x0 = @SVector [xy_new[1], xy_new[2], x0[3], x0[4], x0[5], x0[6]]

        lateral = norm(x0[1:2])
        angle = atand(lateral, x0[3])
        if angle > 45
            println("Iter $i : Angle = $angle")
        end
        # To assert that the error does not drive the rocket into the ground.
        # x0_NoGround = [x for x in [x0[1:2]..., max(x0[3], 0.0), x0[4:6]...]]
        #
        # if x0 != x0_NoGround
        #     println("Iter $i : x0 = $x0 && x0_NoGround = $x0_NoGround")
        # end

        TO.set_initial_state!(prob_mpc, x0)
        X_traj[i+1] = x0

        # Update tracking cost
        TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

        # Shift the initial trajectory
        RD.shift_fill!(prob_mpc.Z)

        # Shift the multipliers and penalties
        Altro.shift_fill!(TO.get_constraints(altro))

        # Again, to ensure they are solving the same problem, we get regenerate
        # the ECOS problem
        ecos, X_ecos, U_ecos = gen_ECOS_Rocket(prob_mpc, Z_track, k_mpc)

        # Solve the updated problem
        Altro.solve!(altro)
        Convex.solve!(ecos, ecos_optimizer)

        iters[i,1] = iterations(altro)
        # iters[i,2] = res.info.iter
        # times[i,1] = median(b).time * 1e-6
        times[i,1] = altro.stats.tsolve
        times[i,2] = ecos_optimizer.sol.solve_time * 1000

        # compare the solutions....
        # Grab the X and U trajectories from ALTRO
        X_altro = getX_toECOS(get_trajectory(altro))
        U_altro = getU_toECOS(get_trajectory(altro))
        # Grab the X and U trajectories from ECOS
        X_ecos_eval = evaluate(X_ecos)
        U_ecos_eval = evaluate(U_ecos)

        # Use the infinity norm to determine the maximum violation between the
        # two trajectories
        err_traj[i,1] = norm(X_altro - X_ecos_eval, Inf)
        err_traj[i,2] = norm(U_altro - U_ecos_eval, Inf)

        ecos_not_optimal = ecos.status != Convex.MathOptInterface.OPTIMAL

        # Get the Euclidean norm beteen the initial states.
        err_x0[i,1] = norm(X_altro[:,1] - x0)
        err_x0[i,2] = norm(X_ecos_eval[:,1] - x0)

        if err_traj[i,1] > 1e-2 || ecos_not_optimal
            println("Iteration $i: State Error = $(err_traj[i,1])")
            plot_3setRef(states(altro), X_ecos,
                    states(Z_track)[k_mpc:k_mpc + prob_mpc.N - 1],
                    title = "Position between ALTRO and ECOS at iter $i")

            # println("ALTRO cost = $(altro.stats.cost[end])")
            # println("ECOS cost  = $(ecos_optimizer.sol.objective_value)")
            # plot_3setRef(controls(altro), U_ecos,
            #         controls(Z_track)[k_mpc:k_mpc + prob_mpc.N - 2],
            #         title = "Controls between ALTRO and ECOS at iter $i")
        end

        altro_cost = sanity_check_cost(X_altro, U_altro,
                                reformat_X_track, reformat_U_track,
                                Qk, Qfk, Rk)
        ecos_cost  = sanity_check_cost(X_ecos_eval, U_ecos_eval,
                                reformat_X_track, reformat_U_track,
                                Qk, Qfk, Rk)
        ref_cost   = sanity_check_cost(reformat_X_track, reformat_U_track,
                                reformat_X_track, reformat_U_track,
                                Qk, Qfk, Rk)
        cost_diff  = altro_cost - ecos_cost

        if cost_diff > 1e-3 || ecos_not_optimal
            println("Iteration $i: ")
            println("ALTRO cost                = $altro_cost")
            println("ECOS cost                 = $ecos_cost")
            println("Cost Difference (~ zero)  = $(altro_cost - ecos_cost)")
        end

        if ref_cost != 0.0
            println("Iteration $i: ")
            println("Ref cost (should be zero) = $ref_cost")
        end

    end
    return X_traj, X_ecos, U_ecos, Dict(:time=>times, :iter=>iters,
                                :err_traj=>err_traj, :err_x0=>err_x0)
end
