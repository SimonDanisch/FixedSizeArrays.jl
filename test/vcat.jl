using FixedSizeArrays

immutable Vec3{T}
	x::T
	y::T
	z::T
end


Vec3{Float32}[1f0, 2f0, 3f0]
Vec3{Float32}[
	1f0, 2f0, 3f0, 
	1f0, 2f0, 3f0
]