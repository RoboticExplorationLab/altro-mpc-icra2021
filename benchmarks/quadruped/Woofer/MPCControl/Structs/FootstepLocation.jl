mutable struct FootstepLocation{T} <: AbstractArray{T, 1}
	fr::SVector{3, T}
	fl::SVector{3, T}
	br::SVector{3, T}
	bl::SVector{3, T}
end

function FootstepLocation(r::AbstractVector{T}) where T
	@assert length(r) == 12

	fr = SVector{3}(r[1], r[2], r[3])
	fl = SVector{3}(r[4], r[5], r[6])
	br = SVector{3}(r[7], r[8], r[9])
	bl = SVector{3}(r[10], r[11], r[12])

	FootstepLocation{T}(fr, fl, br, bl)
end

function footstep_location_from_angles(α::AbstractVector{T}) where T
	r = ForwardKinematicsAll(α)

	return FootstepLocation(r)
end

function Base.getindex(fsl::FootstepLocation, i::Integer)
	if i==1
		return fsl.fr
	elseif i==2
		return fsl.fl
	elseif i==3
		return fsl.br
	elseif i==4
		return fsl.bl
	else
		throw(BoundsError(fsl, i))
	end
end

function Base.setindex!(fsl::FootstepLocation, loc::SVector{3}, i::Integer)
	if i==1
		fsl.fr = loc
	elseif i==2
		fsl.fl = loc
	elseif i==3
		fsl.br = loc
	elseif i==4
		fsl.bl = loc
	else
		throw(BoundsError(fsl, i))
	end
end

Base.convert(::Type{FootstepLocation{T}}, r::AbstractVector) where T = FootstepLocation(r)
Base.convert(::Type{FootstepLocation{T}}, r::FootstepLocation) where T = r
Base.eltype(::Type{FootstepLocation{T}}) where T = SVector{3,T}
Base.length(::FootstepLocation) = 4
Base.size(::FootstepLocation) = (4,)
Base.zero(::Type{FootstepLocation}) = FootstepLocation(@SVector zeros(12))

Base.:*(rot::Rotation, r::FootstepLocation{T}) where T = FootstepLocation{T}(rot*r[1], rot*r[2], rot*r[3], rot*r[4])

Base.:+(v::SVector{3}, r::FootstepLocation{T}) where T = FootstepLocation{T}(v+r[1], v+r[2], v+r[3], v+r[4])
Base.:+(r::FootstepLocation{T}, v::SVector{3}) where T = FootstepLocation{T}(v+r[1], v+r[2], v+r[3], v+r[4])