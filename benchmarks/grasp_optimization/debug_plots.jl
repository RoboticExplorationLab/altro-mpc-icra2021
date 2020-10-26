using Plots
altro_traj = results[1][2]
ecos_controls = results[1][3]

X = altro_traj[:states]
U = altro_traj[:controls]
Ue = ecos_controls
Uc = controls(prob_cold)[1:18]
N = 15

## COMPARE F1
u1 = [U[t][1] for t = 1:N-1]
u2 = [U[t][2] for t = 1:N-1]
u3 = [U[t][3] for t = 1:N-1]
u1e = [Ue[t][1] for t = 1:N-1]
u2e = [Ue[t][2] for t = 1:N-1]
u3e = [Ue[t][3] for t = 1:N-1]
u1c = [Uc[t][1] for t = 1:N-1]
u2c = [Uc[t][2] for t = 1:N-1]
u3c = [Uc[t][3] for t = 1:N-1]

plot([u1 u2 u3 u1e u2e u3e u1c u2c u3c],
    xlabel="Time Step",
    ylabel="Coordinate",
    label = ["u1" "u2" "u3" "u1e" "u2e" "u3e" "u1c" "u2c" "u3c"])

## COMPARE F2
# u1 = [U[t][4] for t = 1:N-1]
# u2 = [U[t][5] for t = 1:N-1]
# u3 = [U[t][6] for t = 1:N-1]
# u1e = [Ue[t][4] for t = 1:N-1]
# u2e = [Ue[t][5] for t = 1:N-1]
# u3e = [Ue[t][6] for t = 1:N-1]
# u1c = [Uc[t][4] for t = 1:N-1]
# u2c = [Uc[t][5] for t = 1:N-1]
# u3c = [Uc[t][6] for t = 1:N-1]
#
# plot([u1 u2 u3 u1e u2e u3e u1c u2c u3c],
#     xlabel="Time Step",
#     ylabel="Coordinate",
#     label = ["u1" "u2" "u3" "u1e" "u2e" "u3e" "u1c" "u2c" "u3c"])

## PLOT STATE TRAJ
# x = [X[t][1] for t = 1:N]
# y = [X[t][2] for t = 1:N]
# z = [X[t][3] for t = 1:N]
# xd = [X[t][4] for t = 1:N]
# yd = [X[t][5] for t = 1:N]
# zd = [X[t][6] for t = 1:N]

# plot([x y z xd yd zd],
#     xlabel="Time Step",
#     ylabel="Coordinate",
#     label = ["x" "y" "z" "xd" "yd" "zd"])

## Plot Timing Results
# altro_times = res[:time][:,1]
# ecos_times = res[:time][:,2]
# bounds = extrema([altro_times; ecos_times])
# bin_min = floor(Int, bounds[1]) - 1
# bin_max = ceil(Int, bounds[2]) + 1
# bins = collect(bin_min:bin_max)
# histogram(altro_times, bins=bins, fillalpha=.5, label="ALTRO")
# histogram!(ecos_times, bins=bins, fillalpha=.5, label="ECOS")
# xlabel!("Solve Time (ms)")
# ylabel!("Counts")
