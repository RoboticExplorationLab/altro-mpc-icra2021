# Create objective for ALTRO or Convex

include("utils.jl")
include("struct_setup.jl")
include("make_problem_ALTRO.jl")



function make_problem(s::selection, r::Rocket, obj_opts::ObjectiveOptions,
                            t_opts::TrajectoryOptions, out_opts::OutputOptions)
    if s.st == USE_ALTRO
        if s.wc == COLD
            return make_problem_ALTRO_COLD(r, o_opts, t_opts, out_opts)
        else
            return make_problem_ALTRO_WARM(r, o_opts, t_opts, out_opts)
        end
    elseif s.st == USE_CONVEX
        if s.wc == COLD
            return make_problem_CONVEX_COLD(r, o_opts, t_opts, out_opts)
        else
            return make_problem_CONVEX_WARM(r, o_opts, t_opts, out_opts)
        end
    end

    error("Solver Not Implemented")
end
