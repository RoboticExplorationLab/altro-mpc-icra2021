import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
import YAML

function change_solver_options(solver="ALTRO", linearized_friction_constraint=true)
    yaml_path = "Woofer/MPCControl/MPC.yaml"
    data = YAML.load(open(yaml_path))
    data["solver"] = solver
    data["linearized_friction_constraint"] = linearized_friction_constraint
    YAML.write_file(yaml_path, data)

    nothing
end

change_solver_options("ALTRO", false)
include("test_solver.jl")
using .SolverTest

tf = 2.0

(altro_times_socp, x_altro_socp, t_altro_socp) = SolverTest.test_solver(tf)
altro_times_socp = altro_times_socp[altro_times_socp .!= 0.0] 

change_solver_options("ALTRO")
include("test_solver.jl")
using .SolverTest

(altro_times, x_altro, t_altro) = SolverTest.test_solver(tf)
altro_times = altro_times[altro_times .!= 0.0] 

change_solver_options("OSQP")
include("test_solver.jl")
using .SolverTest

(osqp_times, x_osqp, t_osqp) = SolverTest.test_solver(tf)
osqp_times = osqp_times[osqp_times .!= 0.0] 

using JLD2

@save "plots/timing_data.jld2" altro_times osqp_times altro_times_socp


using Plots

bins = collect(0:100:2000)

histogram(altro_times_socp ./ 1000, bins=bins, legend=false, xlabel="Solve Times (microsecond)", title="Altro SOCP Solve Times")
png("plots/altro_socp_hist")

histogram(altro_times ./ 1000, bins=bins, legend=false, xlabel="Solve Times (microsecond)", title="Altro Solve Times")
png("plots/altro_hist")

histogram(osqp_times ./ 1000, bins=bins, legend=false, xlabel="Solve Times (microsecond)", title="OSQP Solve Times")
png("plots/osqp_hist")

histogram(altro_times ./ 1000, bins=bins, fillalpha=0.33, xlabel="Solve Times (microsecond)", title="Altro vs OSQP Solve Times", label="Altro")
histogram!(osqp_times ./ 1000, bins=bins, fillalpha=0.33, label="OSQP")
histogram!(altro_times_socp ./ 1000, bins=bins, fillalpha=0.33, label="Altro w/ SOCP")
png("plots/combined_hist")