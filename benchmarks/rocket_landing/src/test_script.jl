# Test Script

using Pkg
Pkg.activate(".\\.")
Pkg.instantiate()

using LinearAlgebra, SparseArrays, StaticArrays
using RobotDynamics, TrajectoryOptimization, Altro

include("make_problem.jl")

r1 = Rocket(10.0, [0; 0; 9.81], 20, 20)
s1 = selection(USE_ALTRO, WARM)
obj1 = ObjectiveOptions(1e-2, 1.0, 1e-1)

x0_new = @SVector [1.0, 1.0, 20.0, 0.0, 0.0, -5.0]
xf_new = @SVector zeros(6)
t1 = TrajectoryOptions(301, 0.5, x0_new, xf_new)

out1 = OutputOptions(true, true, true, true, true)
prob = make_problem_ALTRO_COLD(r1, obj1, t1, out1)
