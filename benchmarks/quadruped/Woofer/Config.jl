import YAML
using StaticArrays
using LinearAlgebra

struct InertialConfig
    frame_mass::Float64
    module_mass::Float64
    upper_link_mass::Float64
    lower_link_mass::Float64
    leg_mass::Float64
    robot_mass::Float64
    sprung_mass::Float64
    body_inertia::SMatrix{3,3,Float64,9}
end

struct ActuatorConfig
    max_joint_torque::Float64
    max_joint_range::Float64
end

struct GeometryConfig
    hip_center_y::Float64 # front-back distance from center line to leg axis
    hip_center_x::Float64  # left-right distance from center line to leg plane
    abduction_offset::Float64 # distance from abduction axis to leg
    foot_radius::Float64

    module_length::Float64
    module_width::Float64
    module_height::Float64

    hip_layout::SMatrix{4,3,Float64,12}
    feet_layout::SMatrix{4,3,Float64,12}
    abduction_layout::SVector{4,Float64}
    body_length::Float64
    body_width::Float64
    body_height::Float64
    upper_link_length::Float64
    lower_link_length::Float64
end


struct WooferConfig
    inertial::InertialConfig
    actuator::ActuatorConfig
    geometry::GeometryConfig
end

function WooferConfig(filename::String)
    data = YAML.load(open(filename))

    # Inertial
    frame_mass = data["inertial"]["frame_mass"]
    module_mass = data["inertial"]["module_mass"]
    upper_link_mass = data["inertial"]["upper_link_mass"]
    lower_link_mass = data["inertial"]["lower_link_mass"]
    leg_mass = (upper_link_mass + lower_link_mass) * 2
    robot_mass = frame_mass + 4 * module_mass + 4 * leg_mass
    sprung_mass = frame_mass + 4 * module_mass + 8 * upper_link_mass
    body_ix = data["inertial"]["body_ix"]
    body_iy = data["inertial"]["body_iy"]
    body_iz = data["inertial"]["body_iz"]
    body_inertia = Diagonal(@SVector [body_ix, body_iy, body_iz])
    inertial = InertialConfig(
        frame_mass,
        module_mass,
        upper_link_mass,
        lower_link_mass,
        leg_mass,
        robot_mass,
        sprung_mass,
        body_inertia,
    )

    # Actuator
    max_joint_torque = data["actuator"]["max_joint_torque"]
    joint_range = data["actuator"]["revolute_range"]
    actuator = ActuatorConfig(max_joint_torque, joint_range)

    # Geometry
    hip_center_y = data["geometry"]["hip_center_y"]
    hip_center_x = data["geometry"]["hip_center_x"]
    abduction_offset = data["geometry"]["abduction_offset"]
    foot_radius = data["geometry"]["foot_radius"]
    module_length = data["geometry"]["module_length"]
    module_width = data["geometry"]["module_width"]
    module_height = data["geometry"]["module_height"]
    hip_layout = @SMatrix ([
        hip_center_x -hip_center_y 0
        hip_center_x hip_center_y 0
        -hip_center_x -hip_center_y 0
        -hip_center_x hip_center_y 0
    ])
    abduction_layout = @SVector([
        -abduction_offset,
        abduction_offset,
        -abduction_offset,
        abduction_offset,
    ])
    feet_layout = hip_layout + ([@SVector(zeros(4)) abduction_layout @SVector(zeros(4))])
    body_length = data["geometry"]["body_length"]
    body_width = data["geometry"]["body_width"]
    body_height = data["geometry"]["body_height"]
    upper_link_length = data["geometry"]["upper_link_length"]
    lower_link_length = data["geometry"]["lower_link_length"]
    geometry = GeometryConfig(
        hip_center_y,
        hip_center_x,
        abduction_offset,
        foot_radius,
        module_length,
        module_width,
        module_height,
        hip_layout,
        feet_layout,
        abduction_layout,
        body_length,
        body_width,
        body_height,
        upper_link_length,
        lower_link_length,
    )
    return WooferConfig(inertial, actuator, geometry)
end

function WooferConfig()
    return WooferConfig(joinpath(@__DIR__,"Woofer.yaml"))
end

const woofer = WooferConfig()
