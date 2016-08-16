__precompile__(true)
module FixedSizeArrays

using Compat

if isdefined(Base, :normalize)
    import Base: normalize
else
    import Compat: normalize
    export normalize
end

importall Base
import Base.LinAlg.chol!

if VERSION >= v"0.5.0-dev+4677"
    import Base.LinAlg._chol!
end


# for 0.5 and 0.4 compat, use our own functor type
abstract Functor{N}

if VERSION < v"0.5.0-dev+698"
    macro pure(ex)
        esc(ex)
    end
else
    import Base: @pure
end

include("core.jl")
include("functors.jl")
include("constructors.jl")

# put them here due to #JuliaLang/julia#12814
# needs to be before indexing and ops, but after constructors
immutable Mat{Row, Column, T} <: FixedMatrix{Row, Column, T}
    values::NTuple{Column, NTuple{Row, T}}
end

# most common FSA types
immutable Vec{N, T} <: FixedVector{N, T}
    values::NTuple{N, T}
end
immutable Point{N, T} <: FixedVector{N, T}
    values::NTuple{N, T}
end

include("mapreduce.jl")
include("destructure.jl")
include("indexing.jl")
include("ops.jl")
include("expm.jl")
include("array_of_fixedsize.jl")
include("conversion.jl")


function show{R,C,T}(io::IO, m::Mat{R,C,T})
	println(io, typeof(m), "(")
	for i=1:R
		println(io, "    ", join(row(m, i), " "))
	end
	println(io, ")")
end

show(io::IO, v::FixedVector{0}) = print(io, typeof(v).name.name, "()")
function show{N,T}(io::IO, v::FixedVector{N,T})
    print(io, typeof(v).name.name, "(", v[1])
    for i = 2:N
        print(io, ",", v[i])
    end
    print(io, ")")
end

export FixedArray
export FixedVector
export FixedMatrix
export MutableFixedArray
export MutableFixedVector
export MutableFixedMatrix
export Mat, Vec, Point
export @fsa
export similar_type
export construct_similar

export unit

export row
export column
export MatMulFunctor
export setindex
export eltype_or, size_or, ndims_or
export @fslice
export destructure
export push, pop, shift, unshift, deleteat, insert




end
