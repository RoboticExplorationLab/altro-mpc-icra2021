X = altro_states
N = 15

x = [X[t][1] for t = 1:N]
y = [X[t][2] for t = 1:N]
z = [X[t][3] for t = 1:N]
xd = [X[t][4] for t = 1:N]
yd = [X[t][5] for t = 1:N]
zd = [X[t][6] for t = 1:N]

using Plots
plot([x y z xd yd zd],
    xlabel="Time Step",
    ylabel="Coordinate",
    label = ["x" "y" "z" "xd" "yd" "zd"])
