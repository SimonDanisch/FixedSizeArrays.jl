function show{FSA <: FixedVector}(io::IO, a::Vector{FSA})
	print(io, "Vector with: ", length(a), "x", FSA, "[")
	for elem in a
		print(io, "[")
		print(io, join(elem, ", "))
		print(io, "]")
	end
	println(io, "]")
end
immutable MaxFun <: Func{2} end
immutable MinFun <: Func{2} end
call(::MaxFun, a, b) = max(a, b)
call(::MinFun, a, b) = min(a, b)
minimum{T <: FixedArray}(a::Vector{T}) = reduce(MinFun(), a)
maximum{T <: FixedArray}(a::Vector{T}) = reduce(MaxFun(), a)
function isapprox{FSA <: FixedArray}(a::FSA, b::Array)
    for i=1:length(a)
        !isapprox(a[i], b[i]) && return false
    end
    true
end

(.+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = T[elem .+ x for elem in a]
(+){T<:FixedArray, ND}(a::Array{T, ND}, x::T) = a .+ x