module FixedSizeArrays

import Quaternions.normalize

include("core.jl")
include("constructors.jl")
include("mapreduce.jl")
include("indexing.jl")
include("ops.jl")
include("array_of_fixedsize.jl")

export FixedArray
export FixedVector
export FixedMatrix
export MutableFixedArray
export MutableFixedVector
export MutableFixedMatrix

export nvec
export normalize
export @gen_fixed_size_vector
export gen_fixed_size_matrix
export fieldname
export row
export column
end