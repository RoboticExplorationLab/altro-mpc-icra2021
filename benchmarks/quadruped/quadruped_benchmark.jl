import Pkg; Pkg.activate(joinpath(@__DIR__,".."));
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/Lyceum/LyceumRegistry.git"))
Pkg.instantiate()
cd(@__DIR__)
import YAML

tf = 2.0

include("mujoco_test.jl")

# ALTRO w/ SOCP Benchmark:
altro_socp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = false)
mpc_dt = altro_socp_controller.mpc_update
(altro_times_socp, x_altro_socp, t_altro_socp) = mujoco_controller_test(
    altro_socp_controller, tf, mpc_dt
)
altro_times_socp = altro_times_socp[altro_times_socp .!= 0.0] 
altro_times_socp .*= 1e-6 # put time into ms
println("Mean ALTRO SOCP solve time: ", mean(altro_times_socp), " ms.")


# ECOS w/ SOCP Benchmark:
ecos_controller = MPCControl.ControllerParams(solver = "ECOS", linearized_friction = false)
(ecos_times, x_ecos, t_ecos) = mujoco_controller_test(ecos_controller, tf, mpc_dt, true)
ecos_times = ecos_times[ecos_times .!= 0.0] 
ecos_times .*= 1e3 # put time into ms
println("Mean ECOS solve time: ", mean(ecos_times), " ms.")


# ALTRO w/ Linearized Friction Cone Benchmark:
altro_qp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = true)
(altro_times, x_altro, t_altro) = mujoco_controller_test(altro_qp_controller, tf, altro_qp_controller.mpc_update)
altro_times = altro_times[altro_times .!= 0.0] 
altro_times .*= 1e-6 # put time into ms
println("Mean ALTRO solve time: ", mean(altro_times), " ms.")


# OSQP w/ Linearized Friction Cone Benchmark:
osqp_controller = MPCControl.ControllerParams(solver = "OSQP", linearized_friction = false)
(osqp_times, x_osqp, t_osqp) = mujoco_controller_test(osqp_controller, tf, mpc_dt)
osqp_times = osqp_times[osqp_times .!= 0.0] 
osqp_times .*= 1e-6 # put time into ms
println("Mean OSQP solve time: ", mean(osqp_times), " ms.")

# using JLD2

# @save "plots/timing_data.jld2" altro_times osqp_times altro_times_socp ecos_times


# using Plots

# bins = collect(0:.10:5)

# histogram(altro_times_socp, bins=bins, legend=false, xlabel="Solve Times (ms)", title="Altro SOCP Solve Times")
# png("plots/altro_socp_hist")

# histogram(altro_times, bins=bins, legend=false, xlabel="Solve Times (ms)", title="Altro Solve Times")
# png("plots/altro_hist")

# histogram(osqp_times, bins=bins, legend=false, xlabel="Solve Times (ms)", title="OSQP Solve Times")
# png("plots/osqp_hist")

# histogram(ecos_times, bins=bins, legend=false, xlabel="Solve Times (ms)", title="ECOS Solve Times")
# png("plots/ecos_hist")

# histogram(altro_times, bins=bins, fillalpha=0.25, xlabel="Solve Times (microsecond)", title="Altro vs OSQP Solve Times", label="Altro")
# histogram!(osqp_times, bins=bins, fillalpha=0.25, label="OSQP")
# histogram!(altro_times_socp, bins=bins, fillalpha=0.25, label="Altro w/ SOCP")
# histogram!(altro_times_socp, bins=bins, fillalpha=0.25, label="ECOS")
# png("plots/combined_hist")