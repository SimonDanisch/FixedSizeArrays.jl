import Base: +

immutable Vec{T}
	x::T
	y::T
	z::T
end
immutable NVec{N, T}
	val::NTuple{N, T}
end
getindex(a::Vec, i) = getfield(a, i)
getindex(a::NVec, i) = a.val[i]
NVec{T}(a::T...) = NVec(a)
stagedfunction map{N, T}(fun::Base.AddFun, a::NVec{N, T}, b::NVec{N, T})
	args = [:(fun(a[$i], b[$i])) for i=1:N]
	:(NVec($(args...)))
end
+{T}(a::Vec{T}, b::Vec{T}) = Vec(a[1]+b[1], a[2]+b[2], a[3]+b[3])
+{N, T}(a::NVec{N, T}, b::NVec{N, T}) = map(Base.AddFun(), a, b)


function test(n, a,b)
	result = a
	for i=1:n
		result = a+b+result
	end
	result
end

const a = NVec(1,2,3)
@time test(1,a , a)
@time test(10^6, a, a)
@time test(10^6, a, a)
println("#############################")
@time test(1, Vec(1,2,3), Vec(2,3,4))
@time test(10^6, Vec(1,2,3), Vec(2,3,4))
@time test(10^6, Vec(1,2,3), Vec(2,3,4))