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

returns (1) the ecos problem, (2) the state trajectory variable, and
(3) the control trajectory variable
"""
function gen_ECOS_Rocket(prob_altro::TrajectoryOptimization.Problem,
                        Z_track=prob_altro.Z,
                        k_track=1;
                        verbose::Bool = false,
                        setStates::Bool = true,
                        setControls::Bool = true,
                        track::Bool = true,
                        includeGoal::Bool = false)

    # Copy the problem to prevent overwritting it
    prob_copy = copy(prob_altro)
    # Move the dynamics into the constraints section
    # TrajectoryOptimization.add_dynamics_constraints!(prob_copy)
    # Convert to NLP to improve a
    # nlp = TrajOptNLP(prob_c, remove_bounds = true);

    # First we recreate the cost function

    # Initialize decision variables
    n, m, N = size(prob_altro)
    X = Variable(n, N)       # State Trajectory
    U = Variable(m, N - 1)   # Control Trajectory
    dt = prob_copy.Z[1].dt

    # First, we build the cost function.
    # Qk, Rk, and Qfk are the Diagonal coefficients of the objective.
    Qk = prob_copy.obj[1].Q[1] * dt
    Rk = prob_copy.obj[1].R[1] * dt
    Qfk = prob_copy.obj[end].Q[1]

    # TO.add_dynamics_constraints!(prob_copy)
    # n,m,N = size(prob_copy)
    # nlp = TrajOptNLP(prob_copy, remove_bounds=true)
    # NN = N*n + (N-1)*m
    #
    # # Cost function
    # TO.hess_f!(nlp)
    # P = nlp.data.G
    # q = zeros(n+m, N)
    # for k = 1:N
    #     q[1:n,k] .= prob.obj[k].q
    #     q[n+1:end,k] .= prob.obj[k].r
    # end
    # dt = prob.Z[1].dt
    # q[:,1:N-1] .*= dt
    # q = q[1:end-m]

    if track
        # We make a cost the penalizes deviations from a reference trajectory
        XTrack = getX_toECOS(Z_track)[:,(k_track-1) .+ (1:N)]
        UTrack = getU_toECOS(Z_track)[:,(k_track-1) .+ (1:N-1)]

        objective = Qk * sumsquares(X[:,1:N-1] - XTrack[:,1:N-1]) +
                        Qfk * sumsquares(X[:,N] - XTrack[:,N]) +
                        Rk * sumsquares(U - UTrack)
    else
        objective = Qk * sumsquares(X[:,1:N-1]) +
                        Qfk * sumsquares(X[:,N]) +
                        Rk * sumsquares(U)
    end

    Z = [vec([X[:,1:end-1]; U]); X[:,N]]
    # objective = quadform(Z, P) + dot(Z,q)
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
    verbose && println("Dynamics Constraint Set")

    # Second up is the goal constraint
    if includeGoal && any(x->x isa GoalConstraint, prob_copy.constraints)
        push!(constraints, X[:, end] == prob_copy.xf)
        verbose && println("Goal Constraint Set")
    end

    # Third up is the ground constraint
    # inds = get_constraint_from_type(prob_copy.constraints,
    #                                     BoundConstraint{1,n + m,Float64})
    # ground_level = prob_copy.constraints[inds[1]].z_min[3]
    # [push!(constraints, X[3, i] >= ground_level) for i in 1:N - 1]
    # verbose && println("Ground Constraint Set at inds: $inds")

    # Fourth up is the max thrust constraint
    inds = get_constraint_from_type(prob_copy.constraints,
            NormConstraint{TrajectoryOptimization.SecondOrderCone,m,Float64})
    if !isempty(inds)
        u_max = prob_copy.constraints[inds[1]].val
        [push!(constraints, norm(U[:,i]) <= u_max) for i in 1:N - 1]
        verbose && println("Max Thrust Constraint Set at inds: $inds")
    else
        verbose && println("Missing Max Thrust Constraint")
    end

    # Fifth up is the max thrust angle constraint
    inds = get_constraint_from_type(prob_copy.constraints,
            NormConstraint2{TrajectoryOptimization.SecondOrderCone,m,m,m})
    if !isempty(inds)
        maxTAalpha = prob_copy.constraints[inds[1]].c[3]
        [push!(constraints, norm(U[1:2, i]) <= maxTAalpha * U[3, i])
                                                            for i in 1:N - 1]
        verbose && println("Max Thrust Angle Constraint Set at inds: $inds")
    else
        verbose && println("Missing Max Thrust Angle Constraint")
    end

    # Sixth up is the glideslope constraint
    inds = get_constraint_from_type(prob_copy.constraints,
            NormConstraint2{TrajectoryOptimization.SecondOrderCone,n,n,n})
    if !isempty(inds)
        glideslope = prob_copy.constraints[inds[1]].c[3]
        [push!(constraints, norm(X[1:2, i]) <= glideslope * X[3, i])
                                                            for i in 1:N - 1]
        # [push!(constraints, cosd(70) * norm(U[1:3, i]) <= U[3, i])
        #                                                     for i in 1:N - 1]
        verbose && println("Glideslope Constraint Set at inds: $inds")
    else
        verbose && println("Missing Glideslope Constraint")
    end
    # glideslope = prob_copy.constraints[inds[2]].c[3]
    # [push!(constraints, cosd(70) * norm(U[1:3, i]) <= U[3, i])
    #                                                         for i in 1:N - 1]

    # Now we are done with the constraints!
    prob_ecos.constraints += constraints


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


    return prob_ecos, X, U


end
