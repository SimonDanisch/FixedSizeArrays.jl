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


#Abstract Types
export FixedArray
export FixedVector
export FixedMatrix
export MutableFixedArray
export MutableFixedVector
export MutableFixedMatrix

#concrete types
export Vec
export Point
export Mat


#Functions outside the AbstractArray interface
export unit
export normalize
export row
export column
export MatMulFunctor

end
