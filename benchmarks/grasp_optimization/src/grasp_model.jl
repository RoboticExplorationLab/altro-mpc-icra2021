include("utils.jl")
import RobotDynamics: dynamics, AbstractModel

# set up model
struct SquareObject{T} <: AbstractModel
    n::Int      # state size
    m::Int      # control size
    mu::T       # coefficient of friction
    mass::T     # mass
    j::T        # inertia
    f::T        # max grasp force
    g::SVector  # gravity
    p::Array
    v::Array
    B::Array
end

function dynamics(model::SquareObject, x, u)
    g = model.g
    m = model.mass

    nqd = Int(model.n/2)
    qd = x[nqd+1:end]

    F1 = u[1:nqd]
    F2 = u[nqd+1:end]
    qdd = F1/m + F2/m + g

    return SVector{model.n}([qd; qdd])
end

RobotDynamics.state_dim(o::SquareObject) = o.n
RobotDynamics.control_dim(o::SquareObject) = o.m
