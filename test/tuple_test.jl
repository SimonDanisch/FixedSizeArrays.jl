importall Base
abstract FixedArray{T, ND, SZ}
abstract FixedVector{T, N} 		  <: FixedArray{T, 1, Tuple{N}}
abstract FixedMatrix{T, Row, Col} <: FixedArray{T, 2, Tuple{Row, Col}}

immutable Vec{N, T} <: FixedVector{T, N}
	a::NTuple{N, T}
end
immutable Mat{Row, Col, T} <: FixedMatrix{T, Row, Col}
	a::NTuple{Row, NTuple{Col, T}}
end
const IndexTypes = Union(Range, AbstractArray)
length{T <: FixedVector}(a::T) = length(a.(1))

Vec{T}(a::T...) = Vec(a)
Mat{N, T}(a::NTuple{N, T}...) = Mat(a)

getindex{T <: FixedVector}(a::T, i::Integer) 					= a.(1)[i]
getindex{T <: FixedVector}(a::T, i::IndexTypes) 				= T.name.primary(a.(1)[i])
getindex{T <: FixedMatrix}(a::T, i::Integer, j::Integer) 		= a.(1)[i][j]
getindex{T <: FixedMatrix}(a::T, i::Integer, j::IndexTypes) 	= a.(1)[i][j]
getindex{T <: FixedMatrix}(a::T, i::IndexTypes, j::Integer) 	= ntuple(k->a.(1)[k][j], length(i))
getindex{T <: FixedMatrix}(a::T, i::IndexTypes, j::IndexTypes) 	= ntuple(k->a.(1)[k][j], length(i))

@show Vec(1,2,3)

immutable Lol3{T}
	x::T
	y::T
	z::T
end
immutable Lol2{T}
	x::T
	y::T
end
getindex(a::Lol3, i) = getfield(a,i)
getindex(a::Lol2, i) = getfield(a,i)
function test(N)
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



@time a = test(1)
@time b = test(10^6)
@time c = test(10^6)
println(a)
println(b)
println(c)

#elapsed time: 0.04126415 seconds (1 MB allocated)
#elapsed time: 1.581399292 seconds (167 MB allocated, 2.00% gc time in 7 pauses with 0 full sweep)
#elapsed time: 1.506825371 seconds (167 MB allocated, 0.31% gc time in 8 pauses with 0 full sweep)

