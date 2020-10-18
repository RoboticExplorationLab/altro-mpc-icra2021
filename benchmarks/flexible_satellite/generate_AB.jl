using LinearAlgebra

function c2d(A,B,dt)
    n = size(A,1)
    p = size(B,2)

    expAB = exp([A*dt B*dt; zeros(p,n+p)])

    A_d = expAB[1:n,1:n]
    B_d = expAB[1:n, (n+1):end ]

    return A_d, B_d
end

function generate_AB()
# inertia matrix
J = diagm([1;2;3])

# reaction wheel jacobian
B_sc = diagm(ones(3))


# linear momentum coupling matrix
phi = [0 1 0;
       1 0 0;
       0 .2 -.8];

# angular momentum coupling matrix
delta = [0 0 1;
         0 1 0;
        -.7 .1 .1]

# store this matrix for faster computations
T = inv(J-delta'*delta)

j = 3; # 3 modes

# damping and stiffness
zeta = [.001;.001;.001]
Delta = [.05; .2; .125] * (2*pi)

# damping and stiffness matrices
C = zeros(j,j)
K = zeros(j,j)
for i =1:j
    C[i,i] = 2*zeta[i]*Delta[i];
    K[i,i] = Delta[i]^2;
end


           #   mrp        w                  n                       ndot
pdot_row = [zeros(3,3) .25*eye(3)       zeros(3,j)                 zeros(3,j)];
wdot_row = [zeros(3,3) zeros(3,3)     T*delta'*K                  T*delta'*C];
ndot_row = [zeros(j,3) zeros(j,3)     zeros(j,j)                  eye(j)];
nddot_row = [zeros(j,3) zeros(j,3) (-K - delta*T*delta'*K)    (-C - delta*T*delta'*C)];

# analytical A
A_analytical = [pdot_row;wdot_row;ndot_row;nddot_row];

# analytical B
B_analytical = [zeros(3,3);
          -T*B_sc;
          zeros(j,3);
          delta*T*B_sc];

# sample time
dt = .5;
Ad, Bd = c2d(A_analytical,B_analytical,dt)

return Ad, Bd
end

# @save "A_B_flexsat" Ad Bd
