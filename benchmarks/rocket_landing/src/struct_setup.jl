#=
This file defines the key structs to ensure that ALTRO.jl and Convex.jl are
indeed solving the same problem
=#

using LinearAlgebra, SparseArrays
include("utils.jl")

"""
    TrajectoryOptions

Convenient struct that holds the necessary constants pertinent to the
trajectory and its discretization.
"""
struct TrajectoryOptions
    N::Int64
    dt::Float64
    x0
    xf
end


"""
    Rocket

Convenient stuct that holds the necessary constants pertinent to the rocket
itself. (Many fields are constructed from a few). We assume a linear dynamics
model.


Acont           Continuous A Matrix (constructed)
Bcont           Continuous B Matrix (constructed)
Adis            Discrete A Matrix (constructed)
Bdis            Discrete B Matrix (constructed)
mass            Mass
grav            Gravity
isp             Specific Impulse
num_controls    Number of Controls (assumed 3)
num_states      Number of States (assumed 6)
theta           Max Thruster Gimbal Angle
u_max           Max Thrust (control)
"""
struct Rocket
    Acont           # Continuous A Matrix (constructed)
    Bcont           # Continuous B Matrix (constructed)
    Gcont           # Continuous G Matrix (constructed)
    Adis            # Discrete A Matrix (constructed)
    Bdis            # Discrete B Matrix (constructed)
    Gdis            # Discrete G Matrix (constructed)
    mass            # Mass
    grav            # Gravity
    isp             # Specific Impulse
    num_controls    # Number of Controls (assumed 3)
    num_states      # Number of States (assumed 6)
    theta           # Max Thruster Gimbal Angle
    u_max           # Max Thrust (control)
end

"""
    Rocket(mass::Float64, grav, theta, u_max, omega = [0.0; 0.0; 0.0];
                        isp = 0.0, num_controls = 3, num_states = 6)

Simplified Constructor for Rocket Struct.

Omega is the rotational velocity of the planet / moon.
"""
function Rocket(mass::Float64, grav, theta, u_max,
                        dt = 1, omega = [0.0; 0.0; 0.0]; isp = 0.0)

    # Write the Dynamics in Continuous Form
    A = [zeros(3,3)     I;
         -S(omega)^2    -2*S(omega)]
    B = [zeros(3,3);    -(1/mass) * I]
    G = [zeros(3,1);    grav]

    rows = size(A, 1)
    cols = size(A, 2) + size(B, 2) + size(G, 2)
    ABG_tot = [A B G; zeros((cols - rows), cols)]

    discrete_mat = exp(ABG_tot * dt)
    rowsA, colsA, colsB = size(A, 1), size(A, 2), size(B, 2)
    # colsA = size(aMat, 2)
    # colsB = size(bMat, 2)

    Ad = sparse(discrete_mat[1:rowsA, 1:colsA])
    Bd = sparse(discrete_mat[1:rowsA, colsA + 1:colsA + colsB])
    Gd = sparse(discrete_mat[1:rowsA, colsA + colsB + 1:end])

    return Rocket(A, B, G, Ad, Bd, Gd, mass, grav, isp, 3, 6, theta, u_max)
end

"""
    ObjectiveOptions

Convenient stuct that holds the necessary constants pertinent to the objective
itself.
"""
struct ObjectiveOptions
    Qk
    Qfk
    Rk
end
