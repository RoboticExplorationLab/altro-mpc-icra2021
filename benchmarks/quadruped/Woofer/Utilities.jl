function LegIndexToRange(i::Int)
    return 3*(i-1)+1:3*(i-1)+3
end

function SLegIndexToRange(i::Int)
    return SVector{3}(3*(i-1)+1:3*(i-1)+3)
end

function Vec12ToSVector(q::Vector)
    @assert length(q) == 12
    return SVector{12}(
        q[1],
        q[2],
        q[3],
        q[4],
        q[5],
        q[6],
        q[7],
        q[8],
        q[9],
        q[10],
        q[11],
        q[12],
    )
end

function Vec3ToSVector(q::Vector)
    @assert length(q) == 3
    return SVector{3}(q[1], q[2], q[3])
end

function Vec4ToSVector(q::Vector)
    @assert length(q) == 4
    return SVector{4}(q[1], q[2], q[3], q[4])
end

"""
Discrete Algebraic Ricatti Equation
- taken from ControlSystems.jl
"""
function dare(A, B, Q, R)
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Ait = try
        inv(A)'
    catch
        error("A must be non-singular.")
    end

    Z = [A + G*Ait*Q   -G*Ait;
         -Ait*Q        Ait]

    S = schur(Z)
    S = ordschur(S, abs.(S.values).<=1)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end

"""
Returns the Discrete LQR gain matrix
- taken from ControlSystems.jl
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end
