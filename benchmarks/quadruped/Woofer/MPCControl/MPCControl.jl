module MPCControl

export
	control!

using LinearAlgebra
using StaticArrays
import StaticArrays: SUnitRange
using Rotations
using ..QuadrupedDynamics
using ForwardDiff
using BenchmarkTools
import YAML

include("../Config.jl")
include("../Utilities.jl")

include("find_solver.jl") # uses MPC.yaml to determine which solver is being used and import correct structs and packages
include("control.jl")
include("swing_leg.jl")
include("gait.jl")
include("footsteps.jl")

end
