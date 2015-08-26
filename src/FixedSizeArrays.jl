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



export FixedArray
export FixedVector
export FixedMatrix
export MutableFixedArray
export MutableFixedVector
export MutableFixedMatrix
export Mat

export unit
export normalize
export row
export column
export MatMulFunctor

end
