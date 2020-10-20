struct OptimizerParams{T, S, L, P, A}
	# discretization length
	dt::T

	# state and control size
	n::S
	m::S

	# ALTRO Variables:
	model::L
	objective::Objective
	constraints::ConstraintList
	problem::P
	solver::A

	X0::Vector{SVector{12, T}}
	U0::Vector{SVector{12, T}}

	u_ref::Vector{SVector{12,T}}
	J::SMatrix{3,3,T,9}
	sprung_mass::T

	function OptimizerParams(
							    dt::T,
							    N::S,
							    q::AbstractVector{T},
							    r::AbstractVector{T},
							    x_des::AbstractVector{T},
							    μ::T,
							    min_vert_force::T,
							    max_vert_force::T;
							    n::Integer = 12,
							    m::Integer = 12,
							) where {T<:Number, S<:Integer}
		Q = Diagonal(SVector{n}(q))
		R = Diagonal(SVector{m}(r))


		# TODO: better initialization here
		u_guess = 9.81*woofer.inertial.sprung_mass/4*(@SVector [0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1.])

		X0 = [x_des for i=1:N]
		U0 = [u_guess for i=1:(N-1)]

		u_ref = [@SVector zeros(T, 12) for i = 1:(N)]

		tf = (N - 1) * dt

		times = collect(0.0:dt:(tf))
		model = RD.LinearModel(n, m, is_affine=true, times=times, dt=times[2])

		constraints = ConstraintList(n,m,N)

		for i=1:4
			add_constraint!(constraints, FrictionConstraint(m, μ, i) , 1:N-1)
		end

		u_min = @SVector [-Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force, -Inf, -Inf, min_vert_force]
		u_max = @SVector [Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force, Inf, Inf, max_vert_force]
		bound = BoundConstraint(n,m, u_min=u_min, u_max=u_max)
		add_constraint!(constraints, bound, 1:N)

		# objective
		objective = LQRObjective(Q, R, Q, x_des, N)

		tf = dt*(N-1)
		problem = Problem(model, objective, x_des, tf, x0=zeros(n), constraints=constraints, integration=RD.PassThrough)
		solver = ALTROSolver(problem)
		set_options!(solver, projected_newton=false, dJ_counter_limit=20)
		set_options!(solver, reset_duals=false, penalty_scaling=10., penalty_initial=10.0)
		set_options!(solver, constraint_tolerance=1e-4, cost_tolerance=1e-4)

		solve!(solver)

		P = typeof(problem)
		A = typeof(solver)
		L = typeof(model)

		new{T,S,L,P,A}(dt, n, m, model, objective, constraints, problem, solver, X0, U0, u_ref, woofer.inertial.body_inertia, woofer.inertial.sprung_mass)
	end
end
