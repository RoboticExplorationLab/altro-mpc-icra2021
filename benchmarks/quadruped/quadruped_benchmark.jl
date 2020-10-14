import YAML

function change_solver(solver="ALTRO")
    yaml_path = "Woofer/MPCControl/MPC.yaml"
    data = YAML.load(open(yaml_path))
    data["solver"] = solver
    YAML.write_file(yaml_path, data)

    nothing
end

change_solver("ALTRO")
include("test_solver.jl")
using .SolverTest

tf = 2.0

(altro_times, x_altro, t_altro) = SolverTest.test_solver(tf)

change_solver("OSQP")
include("test_solver.jl")
using .SolverTest

(osqp_times, x_osqp, t_osqp) = SolverTest.test_solver(tf)

using JLD2

@save "timing_data.jld2" altro_times osqp_times


using Plots

bins = collect(0:100:2000)

histogram(altro_times ./ 1000, bins=bins, legend=false, xlabel="Solve Times (microsecond)", title="Altro Solve Times")
png("altro_hist")

histogram(osqp_times ./ 1000, bins=bins, legend=false, xlabel="Solve Times (microsecond)", title="OSQP Solve Times")
png("osqp_hist")
