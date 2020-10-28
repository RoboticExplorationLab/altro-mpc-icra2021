using Convex, LinearAlgebra, ECOS #, Mosek, MosekTools


function S(vec)
    """Skew from vec"""
    return [0      -vec[3] vec[2];
            vec[3] 0      -vec[1];
            -vec[2] vec[1] 0     ]
end

println("Starting Soft-Landing with ECOS")

# angular velocity of mars
w_mars = [  0 ;  6.62e-5   ;  2.53e-5  ]

# A dynamics matrix for state x = [position;velocity]
A = [zeros(3,3) I;
    -S(w_mars)^2  -2*S(w_mars)]

# B control input matrix
B = [zeros(3,3);I]

# spacecraft mas
sc_mass = 1800.0

# maximym thrust 24,000 N
T_max = 24000

# gravity in state space on Mars
g = 3.71
G = [-g;0;0]

# time step
dt = 1.0

# get discrete dynamics model (this is like c2d in matlab)
exp_discrete = exp([A B; zeros(3,9)]*dt)
Ad = exp_discrete[1:6,1:6]
Bd = exp_discrete[1:6,7:9]

# initial position
initial_position = [2400;450;-330.0]

# initial velocity
initial_velocity = [-10;-40;10]

# goals
final_position = zeros(3)
final_velocity = zeros(3)

T = 60# The number of timesteps
h = copy(dt) # The time between time intervals

# Declare the variables we need
position = Variable(3, T)
velocity = Variable(3, T)
force = Variable(3, T - 1)

# Add dynamics constraints on our variables
constraints = Constraint[ [position[:,i+1];velocity[:, i + 1]] ==  Ad*[position[:, i];velocity[:, i]] + Bd*(force[:,i]/sc_mass + G) for i in 1 : T - 1]

# this is the engine throttle variable
Γ = Variable(T-1)
for i in 1:T-1
    # here we say we want thrust to be in between 20% and 80% power
    push!(constraints, norm(force[:,i]) <= Γ[i])
    push!(constraints,Γ[i] <= .8*T_max)
    push!(constraints,Γ[i] >= .2*T_max)
end

# keep the rocket above the ground
for i in 1:T
    push!(constraints, position[1,i] >= 0.0 )
end

# glide slope angle
glide_slope_angle = deg2rad(45)
tanglide = (tan(glide_slope_angle))

# ensure rocket is within glideslope
E = [0 1 0 ; 0 0 1]
c_cone = [1;0;0]/tanglide
for i in 1:T-1
    push!(constraints, norm(E*(position[:,i] - position[:,T]))-c_cone'*(position[:,i] - position[:,T])<= 0)
end

# pointing constraint, we want the thrust vector to always be within 35
# degrees of vertical (this is so the rocket doesn't flip or something)
n_vec = [1;0;0]
for i in 1:T-1
    push!(constraints, n_vec'*force[:,i] >= cosd(35)*Γ[i])
end

# Add position and velocity constraints
push!(constraints, position[:, 1] == initial_position)
push!(constraints, velocity[:, 1] == initial_velocity)
push!(constraints, velocity[:, T] == final_velocity)
push!(constraints, position[1,T] == 0.0)

# minimze distance from landing spot
# solver = COSMO.Optimizer(max_iter = 9000)
ecos_optimizer = ECOS.Optimizer(
    verbose=0,
    feastol=1e-6,
    abstol=1e-6,
    reltol=1e-6
)

problem = minimize(sumsquares(position[2:3,T]), constraints)
# solve!(problem, () -> Mosek.Optimizer())
# solve!(problem, () -> SDPT3.Optimizer())
# solve!(problem,  () -> COSMO.Optimizer(max_iter = 9000))
Convex.solve!(problem, ecos_optimizer)
println("Solve Time = $(ecos_optimizer.sol.solve_time)")

pos = evaluate(position)

# mat"
# figure
# title('3D position')
# xlabel('X')
# ylabel('Y')
# zlabel('Z')
# hold on
# plot3($pos(1,:),$pos(2,:),$pos(3,:))
# hold off
# "



vel = evaluate(velocity)
u = evaluate(force)

rnorm = zeros(size(pos,2))
vnorm = zeros(size(pos,2))
unorm = zeros(size(force,2))
for i = 1:size(pos,2)
    rnorm[i] = norm(pos[:,i])
    vnorm[i] = norm(vel[:,i])
end
for i = 1:size(force,2)
    unorm[i] = norm(u[:,i])
end


# mat"
# figure
# hold on
# title('Thrust')
# plot($u')
# legend('T_x,T_y,T_z')
# hold off
# "

t_lower = ones(length(unorm))*.2*T_max
t_upper = ones(length(unorm))*.8*T_max

# mat"
# figure
# hold on
# title('Thrust Norm')
# plot($unorm)
# plot(1:length($unorm), $t_lower,'r')
# plot(1:length($unorm), $t_upper,'r')
# legend('Engine Thrust','Minimum Throttle','Maximum Throttle')
# hold off
# "
