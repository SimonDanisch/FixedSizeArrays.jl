function show{FSA <: FixedVector}(io::IO, a::Vector{FSA})
	println(io, length(a), "-element", typeof(a), ":")
	for elem in a
		print(io, " ")
		showcompact(io, elem)
		println(io)
	end
end
immutable MaxFun <: Func{2} end
immutable MinFun <: Func{2} end
immutable ExtremaFun <: Func{2} end
call(::MaxFun, a, b) = max(a, b)
call(::MinFun, a, b) = min(a, b)
call(::ExtremaFun, reducevalue, a) = min(reducevalue[1], a), max(reducevalue[2], a)


minimum{T <: FixedArray}(a::AbstractVector{T}) = reduce(MinFun(), a)
maximum{T <: FixedArray}(a::AbstractVector{T}) = reduce(MaxFun(), a)
Base.extrema{T <: FixedArray}(a::AbstractVector{T}) = reduce(ExtremaFun(), a)
function isapprox{FSA <: FixedArray, A <: Union{Array, FixedArray}}(a::FSA, b::A)
    for i=1:length(a)
        !isapprox(a[i], b[i]) && return false
    end
    true
end

(.+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = T[elem .+ x for elem in a]
(+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = a .+ x
