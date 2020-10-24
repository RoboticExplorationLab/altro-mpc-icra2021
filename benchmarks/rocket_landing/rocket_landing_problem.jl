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

getalpha(θ, deg=true) = deg ? tan(θ) : tand(θ)

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
        θ_max = 7.0,  # deg
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
    x_min = SA[-Inf, -Inf, ground_level, -Inf, -Inf, -Inf]
    TO.add_constraint!(cons, BoundConstraint(n,m, x_min = x_min), 2:N-1)

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
    @assert perWeightMax > 1
    u_bnd = mass * abs(grav[3]) * perWeightMax
    maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
    TO.add_constraint!(cons, maxThrust, 1:N-1)

    # Max Thrust Angle
    #=
    This accounts for the gimbal limits of the thruster.

    We write this as || [u_x; u_y] || ≤ tan(θ) u_z
    To write this in the SOCP form of ||Ax|| ≤ c'x, we take the A and c below.
    =#
    α_max = getalpha(θ_max)
    ARocket = SA_F64[
        1 0 0;
        0 1 0;
        0 0 0
    ]
    cRocket = SA_F64[0, 0, α_max]
    maxAngle = NormConstraint2(n, m, ARocket, cRocket,
                                    TO.SecondOrderCone(), :control)
    TO.add_constraint!(cons, maxAngle, 1:N-1)

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

function run_Rocket_MPC(prob_mpc, opts_mpc, Z_track,
                            num_iters = length(Z_track) - prob_mpc.N)
    # First, Let's generate the ALTRO problem
    altro = ALTROSolver(prob_mpc, opts_mpc)

    # To match the altro problem, we can derive it from the altro
    ecos, X_ecos, U_ecos = gen_ECOS(prob_mpc, verbose = true)

    # Solve initial iteration
    Altro.solve!(altro)

    err_traj = zeros(num_iters,2)
    err_x0 = zeros(num_iters,2)
    iters = zeros(Int, num_iters,2)
    times = zeros(num_iters,2)

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
        x0 += (@SVector randn(n)) * norm(x0,Inf) / 100  # 1% noise
        TO.set_initial_state!(prob_mpc, x0)
        X_traj[i+1] = x0

        # Update tracking cost
        TO.update_trajectory!(prob_mpc.obj, Z_track, k_mpc)

        # Shift the initial trajectory
        RD.shift_fill!(prob_mpc.Z)

        # Shift the multipliers and penalties
        Altro.shift_fill!(TO.get_constraints(altro))

        # Again, to insure they are solving the same problem, we get regenerate
        # the ECOS problem
        ecos, X_ecos, U_ecos = gen_ECOS(prob_mpc)

        # Solve the updated problem
        Altro.solve!(altro)
        Convex.solve!(ecos, ECOS.Optimizer(verbose = 0, feastol=1e-4,
                                            abstol=1e-4, reltol=1e-4))

        iters[i,1] = iterations(altro)
        # iters[i,2] = res.info.iter
        # times[i,1] = median(b).time * 1e-6
        times[i,1] = altro.stats.tsolve
        # times[i,2] = res.info.solve_time * 1000

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

        # Get the Euclidean norm beteen the initial states.
        err_x0[i,1] = norm(X_altro[:,1] - x0)
        err_x0[i,2] = norm(X_ecos_eval[:,1] - x0)

    end
    return X_traj, Dict(:time=>times, :iter=>iters,
                                :err_traj=>err_traj, :err_x0=>err_x0)
end
