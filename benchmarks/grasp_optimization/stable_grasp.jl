using LinearAlgebra
using Convex
using ECOS
using Plots

# parameters
p = [[-1, 0], [1, 0]]   # locations of contact
v = [[1, 0], [-1, 0]]   # normal direction into object
Fg = [0, -1]            # gravity
mu = .1                 # friction constant
m = length(p)           # number of contact points

# variables
t = Variable()
F = Variable(2, m)

# objective
prob = minimize(t)

# friction constraint
for i = 1:m
    global p, F, mu
    prob.constraints += norm((I - v[i]*v[i]')*F[:, i]) <= mu*v[i]'*F[:, i] # friction cone
    prob.constraints += v[i]'*F[:, i] <= t
end

# force balance
prob.constraints += 0 == F[:, 1] + F[:, 2] + Fg

# torque balance
prob.constraints += F[2, 1] == F[2, 2]

# solve
solve!(prob, ECOS.Optimizer)

println(F.value)
println(t.value)
