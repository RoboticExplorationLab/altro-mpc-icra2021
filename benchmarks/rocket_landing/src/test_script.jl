# Test Script

using Pkg
Pkg.activate(".\\.")
Pkg.instantiate()

using LinearAlgebra, SparseArrays, StaticArrays
using RobotDynamics, TrajectoryOptimization, Altro

include("make_problem.jl")
include("..\\plotting\\plot_ALTRO.jl")
include("mpc_altro.jl")


x0_new = @SVector [4.0, 2.0, 20.0, -3.0, 2.0, -5.0]
xf_new = @SVector zeros(6)
t1 = TrajectoryOptionsCOLD(301, 0.05, x0_new, xf_new)

r1 = Rocket(10.0, [0; 0; 9.81], 20, t1.dt)

s1 = selection(USE_ALTRO, COLD)
obj1 = ObjectiveOptions(1e-2, 100.0, 1e-1)

out1 = OutputOptions(true, true, true, true, true)
prob = make_problem_ALTRO_COLD(r1, obj1, t1, out1)

println("Problem Packaged")

opts = SolverOptions(
        cost_tolerance_intermediate=1e-2,
        penalty_scaling=10.,
        penalty_initial=1.0,
        verbose = 0,
        projected_newton = false,
        constraint_tolerance = 1.0e-8
    )

altro = ALTROSolver(prob, opts)
set_options!(altro, show_summary=true)
solve!(altro)

plotALTRO_traj(altro, t1)
plotALTRO_controls(altro, r1)


s2 = selection(USE_ALTRO, WARM)
t2 = convertToWARM(t1, 40, states(altro), controls(altro))

prob = make_problem_ALTRO_WARM(r1, obj1, t2, out1)

runLoop = true

if !runLoop
        altro = ALTROSolver(prob, opts)
        set_options!(altro, show_summary=true)
        solve!(altro)

        plotALTRO_trajMPC(altro, t2)

        # Add some large disturbance and show that it recovers
        x_new_test = states(altro)[1] + [1; 0.2; 0.5; 0.01; 0.02; 0.03]

        println(x_new_test)
        display(states(prob)[1:5])
        # display(controls(prob)[1:5])
        update_MPC_controller!(prob, 2, x_new_test, obj1, t2, out1)
        display(states(prob)[1:5])
        # display(controls(prob)[1:5])

        altro = ALTROSolver(prob, opts)
        set_options!(altro, show_summary=true)
        solve!(altro)

        plotALTRO_trajMPC(altro, t2, k_start = 2)
else
        disturbance = [0.0; 0.0; 0.0; 0.01; 0.02; 0.003]

        out2 = OutputOptions(true, true, false, true, true)
        X_MPC, U_MPC = loop_MPC_controller!(r1, prob, opts, obj1, t2, out2,
                                                disturbance = disturbance)
        plotALTRO_trajMPCLOOP(X_MPC, t2)
end
