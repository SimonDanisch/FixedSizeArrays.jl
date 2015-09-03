abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, Tuple{CARDINALITY}}
typealias MutableFixedMatrix{T, M, N}        MutableFixedArray{T, 2, Tuple{M,N}}

typealias FixedVector{CARDINALITY, T}        FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{Row, Column, T}        FixedArray{T, 2, Tuple{Row, Column}}

abstract FixedVectorNoTuple{CARDINALITY, T} <: FixedVector{CARDINALITY, T}
export FixedVectorNoTuple


# Get the abstract FixedSizeArray type, even for complex type hirarchies
function fixedsizearray_type{FSA <: FixedArray}(::Type{FSA})
    ff = FSA
    while ff.name.name != :FixedArray
        ff = super(ff)
    end
    ff
end
isfullyparametrized{T}(::Type{T}) = !any(x-> isa(x, TypeVar), T.parameters)

eltype{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})         = T
eltype{T <: FixedArray}(A::Type{T})                 = eltype(fixedsizearray_type(T))
eltype{T <: FixedArray}(A::T)                       = eltype(T)

length{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})         = prod(SZ.parameters)::Int
length{T <: FixedArray}(A::Type{T})                 = length(fixedsizearray_type(T))
length{T <: FixedArray}(A::T)                       = length(T)


endof{T,N,SZ}(A::FixedArray{T,N,SZ})                = length(A)


ndims{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})          = N
ndims{T <: FixedArray}(A::Type{T})                  = ndims(fixedsizearray_type(T))
ndims{T <: FixedArray}(A::T)                        = ndims(T)

size{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})           = (SZ.parameters...)::NTuple{N, Int}
size{T <: FixedArray}(A::Type{T})                   = size(fixedsizearray_type(T))
size{T <: FixedArray}(A::T)                         = size(T)

size{T <: FixedArray}(A::T, d::Integer)             = size(T, d)
size{T <: FixedArray}(A::Type{T}, d::Integer)       = size(T)[d]::Int

# Iterator
start(A::FixedArray)                                = 1
next(A::FixedArray, state::Integer)                 = (A[state], state+1)
done(A::FixedArray, state::Integer)                 = length(A) < state



