using GeometryTypes
import ImmutableArrays
using Base.Test


function test(n, a, b)
	r = a*b
	for i=1:n
		r = a*r
	end
	r
end

function test_row(n, a::ImmutableArrays.ImmutableArray)
	r = ImmutableArrays.row(a,1)
	for i=1:n
		for j=1:4
			r = ImmutableArrays.row(a,j)
		end
	end
	r
end
function test_row(n, a)
	r = row(a,1)
	for i=1:n
		for j=1:4
			r = row(a,j)
		end
	end
	r
end
function test_column(n, a::ImmutableArrays.ImmutableArray)
	r = ImmutableArrays.column(a,1)
	for i=1:n
		for j=1:4
			r = ImmutableArrays.column(a,j)
		end
	end
	r
end
function test_column(n, a)
	r = column(a,1)
	for i=1:n
		for j=1:4
			r = column(a,j)
		end
	end
	r
end
function test_sum(n, a)
	r = sum(a)
	for i=1:n
		for j=1:4
			r = sum(a)+j
		end
	end
	r
end

FixedSizeArrays.row(x::Array, i) = vec(x[i, 1:end])
FixedSizeArrays.column(x::Array, i) = x[1:end, i]

function test_dot(n, a)
	r = dot(row(a, 1),column(a, 4))
	for i=1:n
		for j=1:4
			r += dot(row(a, j),column(a, 5-j))
		end
	end
	r
end
function test_dot(n, a::ImmutableArrays.ImmutableArray)
	r = dot(ImmutableArrays.row(a, 1),ImmutableArrays.column(a, 4))
	for i=1:n
		for j=1:4
			r += dot(ImmutableArrays.row(a, j),ImmutableArrays.column(a, 5-j))
		end
	end
	r
end

function test()
	N = 10^6
	const a = [1 1 1 1; 2 2 2 2; 3 3 3 3; 4 4 4 4]
	const b = [1,2,3,4]
	println("matmul: ")
	mmul1 = test(N, a, b)
	@time test(N, a, b)
	test(N, a, a)
	@time test(N, a, a)
	println("sum: ")
	sum1 = test_sum(N, a)
	test_sum(N, a)
	@time test_sum(N, a)
	println("dot:")
	#test_dot(N, a)
	#@time test_dot(N, a)
	println("GeometryTypes: ")
	const a2 = Mat((
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4)
	))
	const b2 = Vec(1,2,3,4)
	for i=1:length(a)
		@test a[i] == a2[i]
	end
	for i=1:4, j=1:4
		@test a[i,j] == a2[i,j]
	end

	println("matmul")
	mmul2 = test(N, a2, b2)
	@time test(N, a2, b2)
	test(N, a2, a2)
	@time test(N, a2, a2)
	println("column:")
	test_column(N, a2)
	@time test_column(N, a2)

	println("row:")
	test_row(N, a2)
	@time test_row(N, a2)

	sum2 = test_sum(N, a2)
	println("sum:")
	test_sum(N, a2)
	@time test_sum(N, a2)

	@test sum2 == sum1

	println("dot:")
	test_dot(N, a2)
	@time test_dot(N, a2)

	println("ImmutableArrays: ")
	const a3 = ImmutableArrays.Matrix4x4(
		ImmutableArrays.Vector4(1,2,3,4),
		ImmutableArrays.Vector4(1,2,3,4),
		ImmutableArrays.Vector4(1,2,3,4),
		ImmutableArrays.Vector4(1,2,3,4)
	)
	const b3 = ImmutableArrays.Vector4(1,2,3,4)
	for i=1:length(a)
		@test a[i] == a3[i]
	end
	for i=1:4, j=1:4
		@test a[i,j] == a3[i,j]
	end

	println("matmul")
	mmul3 = test(N, a3, b3)
	@time test(N, a3, b3)
	
	test(N, a3, a3)
	@time test(N, a3, a3)

	println("column:")
	test_column(N, a3)
	@time test_column(N, a3)

	println("row:")
	test_row(N, a3)
	@time test_row(N, a3)

	sum2 = test_sum(N, a3)
	println("sum:")
	test_sum(N, a3)
	@time test_sum(N, a3)

	@test sum2 == sum1
	@test mmul1 == mmul2
	@test mmul1 == mmul3

	println("dot:")
	test_dot(N, a3)
	@time test_dot(N, a3)

end
@inbounds test()