include("utils.jl")
import RobotDynamics: dynamics, AbstractModel

mutable struct SquareObject{T} <: AbstractModel
    n::Int      # state size
    m::Int      # control size

    mu::T       # coefficient of friction
    mass::T     # mass
    j::T        # inertia
    f::T        # max grasp force
    g::SVector  # gravity

    p::Array    # position vectors of contacts wrt object
    v::Array    # inward pointing surface normals at contacts
    B::Array    # skew matrices for computing torques
    θ::Array  # orientation trajectory
    θdd::Array  # angular acceleration trajectory
end

function SquareObject(mu = 0.5,     # coefficient of friction
                    mass = 0.2,     # kg
                    inertia = 1.,   # kg/m^2
                    f_max = 3.0     # max grasp force
                    )
    n = 6           # state size
    m = 6           # control size
    g = @SVector [0, 0, -9.81]  # gravity

    return SquareObject(n, m, mu, mass, inertia, f_max, g, [], [], [], [], [])
end

function set_orientation_traj!(o::SquareObject, dt, tf;
                             t0 = 0, θ0 = 0, θf= pi/4, θd0 = 0, θdf = .15)
    n = o.n
    m = o.m

    # rotational trajectory
    c = compute_rot_traj_coeffs(t0, tf, [θ0; θf; θd0; θdf])
    o.θ = [dot(c, [t^3,t^2,t,1]) for t = 0:dt:tf]
    o.θdd = [dot(c, [6t,2,0,0]) for t = 0:dt:tf]

    # generate p v B matrices
    p1_0 = [.0,1, 0]; v1_0 = [.0,-1, 0]
    p2_0 = [.0,-1, 0]; v2_0 = [.0,1, 0]
    p2_0 = [.0,1, 0]; v2_0 = [.0,-1, 0]
    p1_0 = [.0,-1, 0]; v1_0 = [.0,1, 0]
    p1, v1, B1 = generate_pvB_3D(p1_0, v1_0, o.θ)
    p2, v2, B2 = generate_pvB_3D(p2_0, v2_0, o.θ)

    # set pvB arrays
    o.p = [p1, p2]
    o.v = [v1, v2]
    o.B = [B1, B2]
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

function RobotDynamics.discrete_dynamics(::Type{PassThrough}, model::SquareObject,
        x::StaticVector, u::StaticVector, t, dt)

    g = model.g
    m = model.mass

    nqd = Int(model.n/2)
    q  = x[SA[1,2,3]]
    qd = x[SA[4,5,6]]

    F1 = u[SA[1,2,3]]
    F2 = u[SA[1,2,3]]
    u = 1/m * (F1 + F2) + g

    q⁺ = q + qd*dt + .5*u*dt^2
    qd⁺ = qd + u*dt

    return [q⁺; qd⁺] 
end
