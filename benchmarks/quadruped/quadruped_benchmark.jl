import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
import YAML

tf = 2.0

include("mujoco_test.jl")

# altro_socp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = false, tol=1e-8)
# ecos_controller = MPCControl.ControllerParams(solver = "ECOS", linearized_friction = false, tol=1e-8)
# test_same_solution(altro_socp_controller, ecos_controller, altro_socp_controller.mpc_update, linearized_friction=false, direct_time = false, test_tol=1e-4)
# altro_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = true, tol=1e-6)
# osqp_controller = MPCControl.ControllerParams(solver = "OSQP", linearized_friction = true, tol=1e-10)
# test_same_solution(altro_controller, osqp_controller, osqp_controller.mpc_update, linearized_friction=true, direct_time = false)

# ALTRO w/ SOCP Benchmark:
altro_socp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = false)
(altro_times_socp, x_altro_socp, t_altro_socp) = mujoco_simulate(altro_socp_controller, tf, altro_socp_controller.mpc_update)
altro_times_socp = altro_times_socp[altro_times_socp .!= 0.0] 
altro_times_socp .*= 1e-6 # put time into ms
println("Mean ALTRO SOCP solve time: ", mean(altro_times_socp), " ms.")


# ECOS w/ SOCP Benchmark:
ecos_controller = MPCControl.ControllerParams(solver = "ECOS", linearized_friction = false)
(ecos_times, x_ecos, u_ecos, t_ecos) = mujoco_simulate(ecos_controller, tf, ecos_controller.mpc_update)
ecos_times = ecos_times[ecos_times .!= 0.0] 
ecos_times .*= 1e3 # put time into ms
println("Mean ECOS solve time: ", mean(ecos_times), " ms.")


# ALTRO w/ Linearized Friction Cone Benchmark:
altro_qp_controller = MPCControl.ControllerParams(solver = "Altro", linearized_friction = true)
(altro_times, x_altro, t_altro) = mujoco_simulate(altro_qp_controller, tf, altro_qp_controller.mpc_update)
altro_times = altro_times[altro_times .!= 0.0] 
altro_times .*= 1e-6 # put time into ms
println("Mean ALTRO solve time: ", mean(altro_times), " ms.")


# OSQP w/ Linearized Friction Cone Benchmark:
osqp_controller = MPCControl.ControllerParams(solver = "OSQP", linearized_friction = false)
(osqp_times, x_osqp, t_osqp) = mujoco_simulate(osqp_controller, tf, osqp_controller.mpc_update)
osqp_times = osqp_times[osqp_times .!= 0.0] 
osqp_times .*= 1e-6 # put time into ms
println("Mean OSQP solve time: ", mean(osqp_times), " ms.")

using JLD2

@save "plots/timing_data.jld2" altro_times osqp_times altro_times_socp ecos_times

# using Plots

# get_states(x,j) = [x[i][j] for i=1:length(x)]

# plot(get_states(x_ecos, 3))
# plot(get_states(x_altro_socp, 3))