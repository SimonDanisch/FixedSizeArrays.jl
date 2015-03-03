include("core.jl")
immutable Matrix{T <: FixedArray} <: WrappedFixedArray{T}
	data::T
end

immutable F3{T} <: FixedArray{T, 1, (3,)}
	x::T
	y::T
	z::T
end
@show Vec(F3(1,2,3))