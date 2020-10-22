using JLD2, Plots
@load string(@__DIR__,"/grasp_ref_traj.jld2") y z F1 F2 θ o_p o_mass o_g

# static version of gif
plot([])
for t = 1:2:N-1
    global y, z, F1, F2, o_p, θ
    local p, F
    p = [o_p[1][t][2:3], o_p[2][t][2:3]]
    F = [F1[t], F2[t]]
    visualize_square([y[t],z[t]], θ[t], p, F, o_mass*o_g[2:3], fa=t/(N-1))
end
plot!([])
