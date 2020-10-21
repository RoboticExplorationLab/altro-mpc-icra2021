# Helpful helper functions

using LinearAlgebra

# Select between solvers
@enum SOLVER_TYPE USE_ALTRO=1 USE_CONVEX=2
@enum WARMCOLD WARM=1 COLD=2

"""
    get_umax(grav, num_gs)

Get the max thrust from the amount of gs
"""
function get_umax(grav, num_gs)
    return grav[3] * num_gs
end

"""
    getAngle3D(u)

Helper function for max thrust angle constraint.
"""
function getAngle3D(u)
    if norm(u[3]) == 0.0
        # hit the tan singularity
        return 0.0
    else
        return atand(norm([u[1]; u[2]]), u[3])
    end
end

"""
    getAngleXY(u)

Helper function for thrust polar plots
"""
function getAngleXY(u)
    return atand(u[2], u[1])
end

"""
    getArrAtInd(arr, ind)

Disaggregate an array
"""
function getArrAtInd(arr, ind)
    return [x[ind] for x in arr]
end


"""
    S(vec)

Skew Symmetric Matrix
"""
function S(vec)
    return [0       -vec[3] vec[2];
            vec[3]  0       -vec[1];
            -vec[2] vec[1]  0     ]
end
