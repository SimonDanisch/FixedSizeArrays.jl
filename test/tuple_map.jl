importall Base


+{N, T}(a::NTuple{N,T}, b::NTuple{N,T}) = map(Base.AddFun(), a, b)
p{T}(a::NTuple{5,T}, b::NTuple{5,T}) = (a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4], a[5]+b[5])

function test(N)
	a=(1,2,3,3,4)
	result = 0.0
	for i=1:N
		tic()
		a + a
		b = toq()
		result += b
	end
	result
end
function test2(N)
	a=(1,2,3,3,4)
	result = 0.0
	for i=1:N
		tic()
		p(a, a)
		b = toq()
		result += b
	end
	result
end

@time test(10^6)
@time test(10^6)
@time test(10^6)

@time test2(10^6)
@time test2(10^6)
@time test2(10^6)
