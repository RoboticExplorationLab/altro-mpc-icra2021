N = 31          # horizon
tf = 3.        # final time
dt = tf/(N-1)   # time step
n = 6           # state size
m = 6           # control size

g = @SVector [0, 0, -9.81]  # gravity
mu = .5                     # friction constant
mass = .2                   # mass
j = 1.                      # inertia
f = 3.                      # max grasp force

# rotational trajectory
θ0 = 0; θf= pi/4; θd0 = 0; θdf = .15; t0 = 0;
c = compute_rot_traj_coeffs(t0, tf, [θ0; θf; θd0; θdf])
θ = [dot(c, [t^3,t^2,t,1]) for t = 0:dt:tf]
θdd = [dot(c, [6t,2,0,0]) for t = 0:dt:tf]

# generate p v B matrices
p1_0 = [.0,1, 0]; v1_0 = [.0,-1, 0]
p2_0 = [.0,-1, 0]; v2_0 = [.0,1, 0]
p1, v1, B1 = generate_pvB_3D(p1_0, v1_0, θ)
p2, v2, B2 = generate_pvB_3D(p2_0, v2_0, θ)

# model
o = SquareObject(n, m, mu, mass, j, f, g, [p1, p2], [v1, v2], [B1, B2])

# initial and final positions
x0 = [0.,3.,3.,0.,0.,0.]
xf = zeros(n)
