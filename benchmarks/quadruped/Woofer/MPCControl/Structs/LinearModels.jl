abstract type LinearModel <: TO.RobotDynamics.AbstractModel end

abstract type AffineModel <: TO.RobotDynamics.AbstractModel end

abstract type AffineTimeInvariantModel <: TO.RobotDynamics.AffineModel end

function discrete_dynamics(::Type{Euler}, model::TO.RobotDynamics.AbstractModel, x::StaticVector, u::StaticVector, t, dt)
	T = eltype(x)
	A_d = oneunit(SizedMatrix{model.n, model.n, T}) + model.A*dt
	B_d = model.B*dt
	d_d = model.d*dt
	A_d * x + B_d * u + d_d
end

function discrete_jacobian!()

abstract type AffineTimeVaryingModel <: TO.RobotDynamics.AffineModel end

function discrete_dynamics(model::AffineTimeVaryingModel, x, u, t)
	k = continuous_to_discrete(t)
	ẋ = model.A[k] * x + model.B[k] * u + model.d[k]
	ẋ
end

continuous_to_discrete(t) = throw(ErrorException("continuous_to_discrete not implemented"))
