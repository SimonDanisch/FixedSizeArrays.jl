tic()
using FixedSizeArrays
toc()
println("#########################################")
immutable Vec{FSV <: FixedVector} <: FixedArrayWrapper{FSV}
	vec::FSV
end
immutable Point{FSV <: FixedVector} <: FixedArrayWrapper{FSV}
	vec::FSV
end
@time Point(1,2,3)
@time Point(1,2,3)

@time Point(1f0,2f0,3f0)
@time Point(1f0,2f0,3f0)

@time Point(1,2)
@time Point(1,2)
const a = Point(1,2)
@show a+a
const b = Point(1f0, 2f0, 4f0, 1f0)
@show sin(a)