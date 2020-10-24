# Convert from ALTRO to ECOS
#=
Since ALTRO has more structure, converting from ALTRO to ECOS is much more
straight forward.
=#

using LinearAlgebra, SparseArrays, StaticArrays
using RobotDynamics, TrajectoryOptimization, Altro
using Convex, ECOS


function get_constraint_from_type(cList, this_type)
    return findall([typeof(c) == this_type for c in cList])
end

function getX_toECOS(Z::Traj)
    return hcat(Vector.(states(Z))...)
end

function getX_toECOS(prob::TO.Problem)
    return getX_toECOS(prob.Z)
end

function getU_toECOS(Z::Traj)
    return hcat([vcat(u...) for u in controls(Z)]...)
end

function getU_toECOS(prob::TO.Problem)
    return getU_toECOS(prob.Z)
end



"""
    gen_ECOS(prob_altro::Problem, opts::SolverOptions{Float64})

Convert ALTRO problem to Convex (ECOS). Note that ALTRO -> ECOS is fairly
straight forward, but ECOS -> ALTRO is non trivial (depending on how the
ECOS constraints are written).
"""
function gen_ECOS(prob_altro::TrajectoryOptimization.Problem,
                        opts::SolverOptions{Float64}; verbose::Bool = false,
                        setStates::Bool = true, setControls::Bool = true,
                        track::Bool = true)

    # Copy the problem to prevent overwritting it
    prob_copy = copy(prob_altro)
    # Move the dynamics into the constraints section
    # TrajectoryOptimization.add_dynamics_constraints!(prob_copy)
    # Convert to NLP to improve a
    # nlp = TrajOptNLP(prob_c, remove_bounds = true);

    # First we recreate the cost function

    n, m, N = size(prob_altro)
    X = Variable(n, N) # State Trajectory
    U = Variable(m, N - 1) # State Trajectory
    dt = prob_copy.Z[1].dt

    # First, we build the cost function.
    # Qk, Rk, and Qfk are the Diagonal coefficients of the objective.
    Qk = prob_copy.obj[1].Q[1]
    Rk = prob_copy.obj[1].R[1]
    Qfk = prob_copy.obj[end].Q[1]

    if track
        # We make a cost the penalizes deviations from a reference trajectory
        XTrack = getX_toECOS(prob_altro)
        UTrack = getU_toECOS(prob_altro)

        objective = Qk * sumsquares(X[:,1:N-1] - XTrack[:,1:N-1]) +
                        Qfk * sumsquares(X[:,N] - XTrack[:,N]) +
                        Rk * sumsquares(U - UTrack)
    else
        objective = Qk * sumsquares(X[:,1:N-1]) +
                        Qfk * sumsquares(X[:,N]) +
                        Rk * sumsquares(U)
    end

    prob_ecos = minimize(objective)

    if verbose
        println("Objective Set")
    end

    # Next, we put the constraints

    # First up are the dynamics constraints
    if RobotDynamics.is_discrete(prob_copy.model)
        if isa(prob_copy.model, LinearizedModel)
            mod = prob_copy.model.linmodel
        else
            mod = prob_copy.model
        end
    else
        lm = LinearizedModel(prob_copy.model, dt = dt, is_affine = true)
        mod = lm.linmodel
    end

    Ad, Bd, Gd = hcat(mod.A...), hcat(mod.B...), [mod.d[1][k] for k in 1:n]

    constraints = Constraint[ X[:,k+1] == Ad * X[:,k] + Bd * U[:,k] + Gd
                                                            for k in 1:N-1 ]

    # Including the initial conditions
    push!(constraints, X[:,1] == prob_copy.x0)

    if verbose
        println("Dynamics Constraint Set")
    end

    # Second up is the goal constraint
    push!(constraints, X[:, end] == prob_copy.xf)
    if verbose
        println("Goal Constraint Set")
    end

    # Third up is the ground constraint
    inds = get_constraint_from_type(prob_copy.constraints,
                                        BoundConstraint{1,n + m,Float64})
    ground_level = prob_copy.constraints[inds[1]].z_min[3]
    push!(constraints, X[3, :] >= ground_level)

    if verbose
        println("Ground Constraint Set")
    end

    # Fourth up is the max thrust constraint
    inds = get_constraint_from_type(prob_copy.constraints,
            NormConstraint{TrajectoryOptimization.SecondOrderCone,m,Float64})
    u_max = prob_copy.constraints[inds[1]].val
    [push!(constraints, norm(U[:,i]) <= u_max) for i in N - 1]

    if verbose
        println("Max Thrust Constraint Set")
    end

    # Fifth up is the max thrust angle constraint
    inds = get_constraint_from_type(prob_copy.constraints,
            NormConstraint2{TrajectoryOptimization.SecondOrderCone,m,m,m})
    maxTAalpha = prob_copy.constraints[inds[1]].c[3]
    [push!(constraints, norm(U[1:2, i]) <= maxTAalpha * U[3, i]) for i in N - 1]

    if verbose
        println("Max Thrust Angle Constraint Set")
    end

    # Now we are done with the constraints!


    # Lastly, we set the initial states and controls
    if setStates
        set_value!(X, getX_toECOS(prob_copy))
        if verbose
            println("Reference State Trajectory Set")
        end
    end

    if setControls
        set_value!(U, getU_toECOS(prob_copy))
        if verbose
            println("Reference Control Trajectory Set")
        end
    end


    return prob_ecos


end
