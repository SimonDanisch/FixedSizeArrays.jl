using GeometryTypes
import ImmutableArrays
using Benchmarks

macro mybench(expr)
	quote
		val = $(esc(expr))
		stats = Statistics(@benchmark for i=1:10
			val = $(esc(expr))
		end)
		println(stats.average_time, ": ", stats.interval[1], " - ", stats.interval[2])
		val
	end
end
function test()
	const a = Float64[1 1 1 1; 2 2 2 2; 3 3 3 3; 4 4 4 4]
	const b = Float64[1,2,3,4]
	const a2 = Mat{4,4,Float64}((
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4)
	))
	const b2 = Vec{4, Float64}(1,2,3,4)
	const a3 = ImmutableArrays.Matrix4x4{Float64}(
		ImmutableArrays.Vector4{Float64}(1,2,3,4),
		ImmutableArrays.Vector4{Float64}(1,2,3,4),
		ImmutableArrays.Vector4{Float64}(1,2,3,4),
		ImmutableArrays.Vector4{Float64}(1,2,3,4)
	)
	const b3 = ImmutableArrays.Vector4{Float64}(1,2,3,4)


	println("matmul 4x4 * 4: ")
	@mybench a*b
	@mybench a2*b2
	@mybench a3*b3

	println("matmul 4x4 * 4x4: ")
	@mybench a*a
	@mybench a2*a2
	@mybench a3*a3

	println("sum: ")
	@mybench sum(a)
	@mybench sum(a2)
	@mybench sum(a3)

	println("dot:")
	@mybench for i=1:100 dot(b,b) end
	@mybench for i=1:100 dot(b2,b2) end
	@mybench for i=1:100 dot(b3,b3) end

	println("column:")
	@mybench for i=1:1000 column(a2,1) end
	@mybench for i=1:1000 ImmutableArrays.column(a3,1) end

	println("row:")
	@mybench for i=1:1000 row(a2,1) end
	@mybench for i=1:1000 ImmutableArrays.row(a3,1) end
end
test()