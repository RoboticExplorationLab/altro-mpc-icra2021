mutable struct ControllerParams{T,S,O}
    # initialize everything once
    mpc_torques
    swing_torques

    prev_phase::S
    cur_phase::S
    cur_phase_time::T # current time relative to beginning of phase
    last_replan_t::T # last replan time
    replan_update::T # time between replanning foot placements

    cur_foot_loc::FootstepLocation{T} # current foot location calculated by FK
    active_feet::SVector{4, S} # active feet on ground
    active_feet_12 # expanded active feet on ground

    planner_foot_loc::FootstepLocation{T} # footstep location for FootstepPlanner
    next_foot_loc::FootstepLocation{T} # actual planned next foot step for SwingLegController

    nom_foot_loc::FootstepLocation{T} # foot location with all joint angles = 0

    N::S # mpc number of time steps
    use_lqr::Bool # use LQR terminal cost to go in optimization
    vel_ctrl::Bool # velocity based control (reference integrates position)

    mpc_update::T # rate at which mpc_forces are updated
    last_t::T # last time that foot forces were calculated


    contacts::Vector{SVector{4,T}} # contact modes over optimization horizon
    foot_locs::Vector{FootstepLocation{T}} # body relative foot locations over optimization horizon
    x_ref::Vector{SVector{12,T}} # state reference trajectory over optimization horizon
    forces::SVector{12,T} # first step of mpc forces
    x_des::SVector{12,T} # desired state for mpc

    optimizer::O
    gait::GaitParams{T, S}
    swing::SwingLegParams{T}
end

function ControllerParams(T, S)
    # TODO: make sure zeros outputs type T
    data = YAML.load(open(joinpath(@__DIR__, "../MPC.yaml")))
    N = data["N"]

    α_0 = @SVector zeros(12)

    mpc_torques = zeros(12)
    swing_torques = zeros(12)

    prev_phase = 1
    cur_phase = 1
    cur_phase_time = 0.0

    last_replan_t = 0.0
    replan_update = data["footstep_replan"]

    cur_foot_loc = footstep_location_from_angles(α_0)
    active_feet = @SVector zeros(S, 4)
    active_feet_12 = zeros(S, 12)

    planner_foot_loc = footstep_location_from_angles(α_0)
    next_foot_loc = footstep_location_from_angles(α_0)

    # ensures that foot forces are calculated at start
    last_t = -1

    contacts = [@SVector zeros(S, 4) for _ = 1:(N+1)]
    foot_locs = [footstep_location_from_angles(α_0) for _ = 1:(N+1)]

    x_ref = [@SVector zeros(12) for _ = 1:(N+1)]

    forces = SVector{12}(
        [0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0, 0, 0, 1.0] *
        woofer.inertial.sprung_mass *
        9.81 / 4,
    )

    ψ = data["yaw_angle"]

    x_des = SVector{12}([
        0.0;
        0.0;
        data["stance_height"];
        zeros(2);
        ψ;
        data["xy_vel"];
        zeros(3);
        data["omega_z"]
    ])

    # TODO: use IK to make sure nominal foot location respects stance height
    nom_foot_loc_arr = ForwardKinematicsAll(zeros(12))
    offset = [1 -1 1 -1]
    Δx = data["foot_dx"]
    Δy = data["foot_dy"]
    for i = 1:4
        nom_foot_loc_arr[LegIndexToRange(i)] += [Δx, Δy * offset[i], 0]
    end

    nom_foot_loc = FootstepLocation(nom_foot_loc_arr)

    μ = data["mu"]
    min_vert_force = data["min_vert_force"]
    max_vert_force = data["max_vert_force"]

    if using_altro
        optimizer = OptimizerParams(
            data["dynamics_discretization"],
            N,
            data["q"],
            data["r"],
            x_des,
            μ,
            min_vert_force,
            max_vert_force,
        )
    else
        optimizer = OptimizerParams(
            data["dynamics_discretization"],
            N,
            data["q"],
            data["r"],
            μ,
            min_vert_force,
            max_vert_force,
        )
    end

    gait_type = data["gait"]["type"]

    if gait_type == "trot"
        gait = trot(
            stance_time = data["gait"]["stance_time"],
            swing_time = data["gait"]["swing_time"],
        )
    elseif gait_type == "pronk"
        gait = pronk(
            stance_time = data["gait"]["stance_time"],
            flight_time = data["gait"]["swing_time"],
        )
    elseif gait_type == "pace"
        gait = pace(
            stance_time = data["gait"]["stance_time"],
            swing_time = data["gait"]["swing_time"],
        )
    elseif gait_type == "flying_trot"
        gait = flying_trot(
            stance_time = data["gait"]["stance_time"],
            flight_time = data["gait"]["swing_time"],
        )
    else
        gait = stand()
    end

    swing = SwingLegParams(
        data["swing"]["step_height"],
        data["swing"]["omega"],
        data["swing"]["zeta"],
    )

    O = typeof(optimizer)

    ControllerParams{T,S,O}(
        mpc_torques,
        swing_torques,
        prev_phase,
        cur_phase,
        cur_phase_time,
        last_replan_t,
        replan_update,
        cur_foot_loc,
        active_feet,
        active_feet_12,
        planner_foot_loc,
        next_foot_loc,
        nom_foot_loc,
        N,
        data["use_lqr"],
        data["velocity_control"],
        data["update_dt"],
        last_t,
        contacts,
        foot_locs,
        x_ref,
        forces,
        x_des,
        optimizer,
        gait,
        swing,
    )
end
