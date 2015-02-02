# Wrapper type, for types that have an immutable as an
abstract WrappedFixedArray{T <: AbstractFixedArray}	
getindex        (A::WrappedFixedArray,            inds::Real...) = A.(1)[inds...]
