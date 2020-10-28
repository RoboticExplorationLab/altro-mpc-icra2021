using Altro
using TrajectoryOptimization

const TO = TrajectoryOptimization

"""
Create a Trajectory Optimization problem that tracks the trajectory in `prob`,
using the same constraints, minus the goal constraint. Tracks the first `N`
time steps.
"""
function gen_tracking_problem(prob::TO.Problem, N;
        Qk = 10.0,
        Rk = 0.1,
        Qfk = Qk,
    )
    n,m = size(prob)
    dt = prob.Z[1].dt
    tf = (N-1)*dt

    # Get sub-trajectory
    Z = Traj(prob.Z[1:N])
    x0 = state(Z[1])
    xf = state(Z[N])  # this actually doesn't effect anything

    # Generate a cost that tracks the trajectory
    Q = Diagonal(@SVector fill(Qk, n))
    R = Diagonal(@SVector fill(Rk, m))
    Qf = Diagonal(@SVector fill(Qfk, n))
    obj = TO.TrackingObjective(Q, R, Z, Qf=Qf)

    # Use the same constraints, except the Goal constraint
    cons = ConstraintList(n,m,N)
    for (inds, con) in zip(prob.constraints)
        if !(con isa GoalConstraint)
            if inds.stop > N
                inds = inds.start:N-(prob.N - inds.stop)
            end
            length(inds) > 0 && TO.add_constraint!(cons, con, inds)
        end
    end

    prob = TO.Problem(prob.model, obj, xf, tf, x0=x0, constraints=cons,
        integration=TO.integration(prob)
    )
    initial_trajectory!(prob, Z)
    return prob
end
