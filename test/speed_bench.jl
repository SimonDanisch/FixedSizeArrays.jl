using GeometryTypes
import ImmutableArrays
using Benchmarks

macro mybench(expr)
	funsym = gensym()
	quote
		function $funsym()
			val = $(esc(expr))
			for i=1:10000
				val += $(esc(expr))
			end
			val
		end
		res = @benchmark $funsym()
		stats = Benchmarks.SummaryStatistics(res).elapsed_time_center
		println(stats)
	end
end
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
	@mybench dot(b,b)
	@mybench dot(b2,b2)
	@mybench dot(b3,b3)

	println("column:")
	@mybench column(a2,1)
	@mybench ImmutableArrays.column(a3,1)

	println("row:")
	@mybench row(a2,1)
	@mybench ImmutableArrays.row(a3,1)
