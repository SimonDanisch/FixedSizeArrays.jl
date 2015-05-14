importall Base

@generated function getindex{N,T, Tupl <: Tuple}(a::NTuple{N,T}, i::Type{Tupl})
	range = Tupl.parameters[1] #Of form Tuple{1:3}
	:(tuple($([:(a[$k]) for k in range]...)))
end

const index = Tuple{1:2}
const a=(1,2,3,4)
@code_llvm a[Tuple{1:2}]
function test(N)
	a=(1,2,3,4)
	result = 0.0
	for i=1:N
		tic()
		rsult = a[index]
		b = toq()
		result += b
	end
	result
end
function test2(N)
	a=(1,2,3,4)
	result = 0.0
	for i=1:N
		tic()
		rsult = a[1:2]
		b = toq()
		result += b
	end
	result
end
function test3(N)
	a=(1,2,3,4)
	result = 0.0
	for i=1:N
		tic()
		rsult = (a[1], a[2])
		b = toq()
		result += b
	end
	result
end

function test4(N)
	a=(1,2,3,4)
	result = 0.0
	for i=1:N
		tic()
		rsult = a[Tuple{1:2}] # allocating Tuple{1:2} seems to make a difference
		b = toq()
		result += b
	end
	result
end

println("optimum?!")
@time test3(10^6)
@time test3(10^6)
@time test3(10^6)

println("tuple standard")
@time test2(10^6)
@time test2(10^6)
@time test2(10^6)

println("const type")
@time test(10^6)
@time test(10^6)
@time test(10^6)
@time test(10^6)
@time test(10^6)
@time test(10^6)

println("type 2 not const")
@time test4(10^6)
@time test4(10^6)
@time test4(10^6)
@time test4(10^6)
@time test4(10^6)
@time test4(10^6)

