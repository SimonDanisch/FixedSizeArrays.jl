module FixedSizeArrays

export RGB
include("core.jl")
include("staged.jl")
include("ops.jl")
include("array_of_fixedsize.jl")

export AbstractFixedArray
export AbstractFixedVector
export AbstractFixedMatrix
export nvec


end