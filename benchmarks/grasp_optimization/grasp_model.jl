import RobotDynamics: dynamics, AbstractModel

# set up model
struct SquareObject{T} <: AbstractModel
    mu::T  # coefficient of friction
    m::T  # mass
    j::T   # inertia
    g::Array{T, 1}   # gravity
    p::Array
    v::Array
    B::Array
end

function dynamics(model::SquareObject, x, u)
    g = model.g
    m = model.m

    q = x[ @SVector [1,2] ]
    qd = x[ @SVector [3,4] ]

    B = @SMatrix [1/m 0 1/m 0;
                  0 1/m 0 1/m]
    qdd = B*u + g
    return [qd; qdd]
end

RobotDynamics.state_dim(::SquareObject) = 4
RobotDynamics.control_dim(::SquareObject) = 4
