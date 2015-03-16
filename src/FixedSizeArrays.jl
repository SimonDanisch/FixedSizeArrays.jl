module FixedSizeArrays

include("core.jl")
include("constructors.jl")
include("mapreduce.jl")
include("indexing.jl")
include("ops.jl")
include("array_of_fixedsize.jl")

export FixedArray
export FixedArrayWrapper
export FixedVector
export FixedMatrix
export nvec
export normalize
export gen_fixed_sizea_array
end