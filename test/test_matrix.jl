using FixedSizeArrays

type A{T} <: FixedVector{T, 3}
	x::T
	y::T 
	z::T
end
type B{T} <: FixedVector{T, 3}
	f::T
	d::T 
	g::T
end

const a = A(1,2,3)
const b = A(1,2,3)

@time a+a
@time a+a

@time b+b
@time b+b
