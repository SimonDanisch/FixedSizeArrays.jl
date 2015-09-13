VERSION >= v"0.4.0-dev+6521" && __precompile__(true)
module FixedSizeArrays

importall Base
import Base.Func

include("core.jl")
include("functors.jl")
include("constructors.jl")

# put them here due to #JuliaLang/julia#12814
# needs to be befor indexing and ops, but after constructors
immutable Mat{Row, Column, T} <: FixedMatrix{Row, Column, T}
    _::NTuple{Column, NTuple{Row, T}}
end
size_or{R, C, T}(::Type{Mat{R, C, T}}, OR) = (R, C)
size_or{T <: Mat}(::Type{T}, OR) = OR
# most common FSA types
immutable Vec{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
end
immutable Point{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
end

include("mapreduce.jl")
include("indexing.jl")
include("ops.jl")
include("array_of_fixedsize.jl")
include("conversion.jl")


function show{R,C,T}(io::IO, m::Mat{R,C,T})
	println(io, typeof(m), "(")
	for i=1:R
		println(io, "    ", join(row(m, i), " "))
	end
	println(io, ")")
end

export FixedArray
export FixedVector
export FixedMatrix
export MutableFixedArray
export MutableFixedVector
export MutableFixedMatrix
export Mat, Vec, Point
export @fsa

export unit
export normalize
export row
export column
export MatMulFunctor
export setindex

end
