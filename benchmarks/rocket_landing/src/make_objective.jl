# Create objective for ALTRO or Convex

include("utils.jl")
include("struct_setup.jl")


function make_objective(st::SOLVER_TYPE, wc::WARMCOLD, r::Rocket,
                        o_opts::ObjectiveOptions, t_opts::TrajectoryOptions)
    if st == USE_ALTRO
        return make_objective_ALTRO(wc, r, o_opts, t_opts)
    elseif st == USE_CONVEX
        return make_objective_CONVEX(wc, r, o_opts, t_opts)
    end

    error("Solver Not Implemented")
end

function make_objective_ALTRO(wc, r, o_opts, t_opts)
    println("Chose ALTRO")
end

function make_objective_CONVEX(wc, r, o_opts, t_opts)
    println("Chose CONVEX")
end
