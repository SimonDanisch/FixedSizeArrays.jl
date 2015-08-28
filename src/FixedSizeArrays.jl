__precompile__(true)
module FixedSizeArrays
using Compat

importall Base
import Base.Func

include("core.jl")
include("functors.jl")
include("constructors.jl")
include("mapreduce.jl")
include("indexing.jl")
include("ops.jl")
include("array_of_fixedsize.jl")
include("conversion.jl")

# put them here due to #JuliaLang/julia#12814
# most common FSA types
immutable Vec{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
end
immutable Point{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
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
