#=
This file defines the key structs to ensure that ALTRO.jl and Convex.jl are
indeed solving the same problem
=#

struct Rocket
    Acont
    Bcont
    Adis
    Bdis
    m
    grav
    isp
end
