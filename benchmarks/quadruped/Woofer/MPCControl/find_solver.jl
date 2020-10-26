import YAML

const data = YAML.load(open(joinpath(@__DIR__, "MPC.yaml")))

if data["solver"] == "ALTRO"
	const using_altro = true

else
	const using_altro = false
end

include("Structs/FootstepLocation.jl")
include("Structs/ALTROParams.jl")
include("Structs/LinearizedFrictionConstraint.jl")
include("Structs/FrictionConstraint.jl")
include("Structs/OSQPParams.jl")
include("Structs/ECOSParams.jl")
include("Structs/SwingLegParams.jl")
include("Structs/GaitParams.jl")
include("Structs/ControllerParams.jl")
include("linearized_dynamics.jl")
include("osqp_solver.jl")
include("altro_solver.jl")
include("ecos_solver.jl")