#=
This file defines the key structs to ensure that ALTRO.jl and Convex.jl are
indeed solving the same problem
=#

using LinearAlgebra, SparseArrays, StaticArrays
include("utils.jl")

# Select between solvers
@enum SOLVER_TYPE USE_ALTRO=1 USE_CONVEX=2
@enum WARMCOLD WARM=1 COLD=2


"""
    selection

Convenient struct that holds the user selection in solver.
"""
struct selection
    st::SOLVER_TYPE
    wc::WARMCOLD
end

"""
    output_options

Convenient struct that holds the printing and plotting options.
"""
struct OutputOptions
    save_data::Bool
    plot_during::Bool
    verbose::Bool
    benchmark_per_solve::Bool
    benchmark_per_outer_loop::Bool
end

"""
    TrajectoryOptions

Convenient struct that holds the necessary constants pertinent to the
trajectory and its discretization.
"""
abstract type TrajectoryOptions end

"""
    TrajectoryOptionsCOLD

TrajectoryOptions for Cold Start Cases.
"""
struct TrajectoryOptionsCOLD <: TrajectoryOptions
    N::Int64
    dt::Float64
    x0
    xf
end

"""
    TrajectoryOptionsWARM

TrajectoryOptions for Warm Start Cases.
"""
struct TrajectoryOptionsWARM <: TrajectoryOptions
    Horizon::Int64
    N::Int64
    dt::Float64
    x0
    xf
    ref_traj_x
    ref_traj_u
end

get_tf(t_opt::TrajectoryOptions) = t_opt.N * t_opt.dt
get_tf_horizon(t_opt::TrajectoryOptionsWARM) = t_opt.Horizon * t_opt.dt

convertToWARM(t_COLD::TrajectoryOptions, Horizon, ref_traj_x, ref_traj_u) =
                            TrajectoryOptionsWARM(Horizon, t_COLD.N, t_COLD.dt,
                                t_COLD.x0, t_COLD.xf, ref_traj_x, ref_traj_u)


"""
    ObjectiveOptions

Convenient stuct that holds the necessary constants pertinent to the objective
itself (specifically for LQR-like objectives).
"""
struct ObjectiveOptions
    Qk::Float64
    Qfk::Float64
    Rk::Float64
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

Theta is in degrees. Omega is the rotational velocity of the planet / moon.
"""
function Rocket(mass::Float64, grav, theta, dt = 1, omega = [0.0; 0.0; 0.0];
                                isp = 0.0, perWeightMax = 2)

    # Write the Dynamics in Continuous Form
    # Alternatively, we could include gravity with u (and hence B)
    A = [zeros(3,3)     I;
         -S(omega)^2    -2*S(omega)]
    B = [zeros(3,3);    -(1/mass) * I]
    G = [zeros(3,1);    grav]

    # Form as a single large matrix
    rows = size(A, 1)
    cols = size(A, 2) + size(B, 2) + size(G, 2)
    ABG_tot = [A B G; zeros((cols - rows), cols)]

    # Exponentiate the Matrix and convert to discrete
    discrete_mat = exp(ABG_tot * dt)

    # Extract Discrete Matrices
    rowsA, colsA, colsB = size(A, 1), size(A, 2), size(B, 2)
    Ad = sparse(discrete_mat[1:rowsA, 1:colsA])
    Bd = sparse(discrete_mat[1:rowsA, colsA + 1:colsA + colsB])
    Gd = sparse(discrete_mat[1:rowsA, colsA + colsB + 1:end])

    u_max = mass * grav[3] * perWeightMax

    # We assume 3 controls and 6 states, otherwise the above does not hold.
    return Rocket(A, B, vec(G), Ad, Bd, Gd, mass, grav, isp, 3, 6, theta, u_max)
end

"""
    propogate_dynamics(r::Rocket, state, controls, disturbance = zeros(6, 1))

Propogate the system forward by one time step.
"""
function propogate_dynamics(r::Rocket, state, controls,
                                disturbance = zeros(6, 1))
    return r.Adis * state + r.Bdis * controls + r.Gdis + disturbance
end

function get_traj_size(r::Rocket, t_opt::TrajectoryOptions)
    # Last step does not have a control
    return t_opt.N * (r.num_controls + r.num_states) + r.num_states
end
