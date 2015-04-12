

import Base: +
immutable Vec{T}
	x::T
	y::T
	z::T
end
immutable Point{N,T}
	val::NTuple{N,T}
end

point{T}(a::T...) = Point(a)

getindex(a::Vec, i) = getfield(a, i)
getindex(a::Point, i) = a.val[i]

function +(a::Tuple{Int64, Int64, Int64}, b::Tuple{Int64, Int64, Int64})
    Base.llvmcall("""%3 = add [4 x i32] %1, %0
      ret [4 x i32] %3""",Tuple{Int64, Int64, Int64},
      (Tuple{Int64, Int64, Int64},Tuple{Int64, Int64, Int64}),
        a,
      b)
end
#+{N, T}(a::NTuple{N, T}, b::NTuple{N, T}) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
+{T}(a::Vec{T}, b::Vec{T}) = Vec(a[1]+b[1], a[2]+b[2], a[3]+b[3])
+{N, T}(a::Point{N, T}, b::Point{N, T}) = Point(a.val+b.val)


function test(n, a,b)
	result = a
	for i=1:n
		result = a+b+result
	end
	result
end
function test2V(n)
	result = Vec(2,3,4)
	for i=1:n
		a,b,c = rand(Int, 3)
		tmp = Vec(a,b,c)
		result = tmp + result
	end
	result
end
function test2P(n)
	result = point(2,3,4)
	for i=1:n
		a,b,c = rand(Int, 3)
		tmp = point(a,b,c)
		result = tmp + result
	end
	result
end
function test2T(n)
	result = (2,3,4)
	for i=1:n
		a,b,c = rand(Int, 3)
		tmp = (a,b,c)
		result = tmp + result
	end
	result
end
@time test(1, (1,2,3), (2,3,4))
@time test(10^6, (1,2,3), (2,3,4))
@time test(10^6, (1,2,3), (2,3,4))
println("#############################")
@time test(1, Vec(1,2,3), Vec(2,3,4))
@time test(10^6, Vec(1,2,3), Vec(2,3,4))
@time test(10^6, Vec(1,2,3), Vec(2,3,4))
println("#############################")

@time test(1, point(1,2,3), point(2,3,4))
@time test(10^6, point(1,2,3), point(2,3,4))
@time test(10^6, point(1,2,3), point(2,3,4))

println("#############################")

@time test2T(1)
@time test2T(10^6)
@time test2T(10^6)
println("#############################")

@time test2P(1)
@time test2P(10^6)
@time test2P(10^6)
println("#############################")

@time test2V(1)
@time test2V(10^6)
code_llvm(+, Tuple{NTuple{3, Int}, NTuple{3, Int}})

ss() = (1,2,3)
code_llvm(ss, Tuple{})
#elapsed time: 0.223158449 seconds (282 MB allocated, 2.78% gc time in 13 pauses with 0 full sweep)
#elapsed time: 0.223116501 seconds (282 MB allocated, 2.71% gc time in 13 pauses with 0 full sweep)


#elapsed time: 0.317912366 seconds (285 MB allocated, 2.16% gc time in 13 pauses with 0 full sweep)
#elapsed time: 0.220877232 seconds (282 MB allocated, 2.74% gc time in 13 pauses with 0 full sweep)
#elapsed time: 0.219239059 seconds (282 MB allocated, 2.76% gc time in 13 pauses with 0 full sweep)