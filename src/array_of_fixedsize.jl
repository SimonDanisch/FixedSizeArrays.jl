function show{FSA <: FixedVector}(io::IO, a::Vector{FSA})
	println(io, length(a), "-element", typeof(a), ":")
	for elem in a
		print(io, " ")
		showcompact(io, elem)
		println(io)
	end
end

immutable MaxFunctor <: Functor{2} end
immutable MinFunctor <: Functor{2} end
immutable ExtremaFun <: Functor{2} end

@compat (::MaxFunctor)(a, b) = max(a, b)
@compat (::MinFunctor)(a, b) = min(a, b)
@compat (::ExtremaFun)(reducevalue, a) = min(reducevalue[1], a), max(reducevalue[2], a)

minimum{T <: FixedArray}(a::Vector{T}) = reduce(MinFunctor(), a)
maximum{T <: FixedArray}(a::Vector{T}) = reduce(MaxFunctor(), a)
extrema{T <: FixedArray}(a::AbstractVector{T}) = reduce(ExtremaFun(), a)


function isapprox{FSA <: FixedArray, A <: FixedArray}(x::FSA, y::A; rtol::Real=Base.rtoldefault(eltype(x),eltype(y)), atol::Real=0, norm::Function=vecnorm)
    # Same behaviour as AbstractArray: julia/base/linalg/generic.jl
	d = norm(x - y)
    return isfinite(d) ? d <= atol + rtol*max(norm(x), norm(y)) : x == y
end
function isapprox{FSA <: FixedArray, A <: Array}(x::FSA, y::A; rtol::Real=Base.rtoldefault(eltype(x),eltype(y)), atol::Real=0, norm::Function=vecnorm)
    d = norm(Array(x) - y) # operations between Arrays and FixedArrays not defined...
    return isfinite(d) ? d <= atol + rtol*max(norm(x), norm(y)) : x == y
end

(.+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = T[elem .+ x for elem in a]
(+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = a .+ x
