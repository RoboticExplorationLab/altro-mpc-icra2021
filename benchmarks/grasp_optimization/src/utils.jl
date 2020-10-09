rot(θ) = [cos(θ) -sin(θ);
        sin(θ) cos(θ)]

rot3(θ) = [1    0       0;
           0 cos(θ) -sin(θ);
           0 sin(θ) cos(θ)]

skew(a) = [  0  -a[3] a[2];
           a[3]    0 -a[1];
          -a[2]  a[1]   0  ]

function generate_pvB_3D(p0, v0, θ)
    T = length(θ)
    p = [rot3(θ[t])*p0 for t = 1:T]
    v = [rot3(θ[t])*v0 for t = 1:T]
    B = [skew(p[t]) for t = 1:T]
    return p, v, B
end

function compute_rot_traj_coeffs(t0,tf,c)
    A = [t0^3 t0^2 t0 1;
         tf^3 tf^2 tf 1;
         3t0^2 2t0 1 0;
         3tf^2 2tf 1 0]
    A\c
end

# rotational trajectory
# θ0 = 0; θf= pi/4; θd0 = 0; θdf = .1
# t0 = 0; tf = (T-1)*dt
# c = compute_rot_traj_coeffs(t0, tf, [θ0; θf; θd0; θdf])
# θ = [dot(c, [t^3,t^2,t,1]) for t = 0:dt:tf]
# θdd = [dot(c, [6t,2,0,0]) for t = 0:dt:tf]
