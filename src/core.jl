abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, Tuple{CARDINALITY}}
typealias MutableFixedMatrix{T, M, N}        MutableFixedArray{T, 2, Tuple{M,N}}

typealias FixedVector{CARDINALITY, T}        FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{Row, Column, T}        FixedArray{T, 2, Tuple{Row, Column}}

abstract FixedVectorNoTuple{CARDINALITY, T} <: FixedVector{CARDINALITY, T}
export FixedVectorNoTuple


_length{T <: Tuple}(::Type{T})						= *(T.parameters...)
_length{N, N2}(::Type{Tuple{N, N2}})				= N*N2
_length{N}(::Type{Tuple{N}})						= N

_size{T <: Tuple}(::Type{T})						= (T.parameters...)
_size{N, N2}(::Type{Tuple{N, N2}})					= (N,N2)
_size{N}(::Type{Tuple{N}})							= (N,)

eltype{T <: FixedArray}(A::Type{T})                 = eltype_or(T, Any)
eltype{T <: FixedArray,N,SZ}(A::FixedArray{T,N,SZ}) = T


length{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})         = _length(SZ)
length{T,N,SZ}(::FixedArray{T,N,SZ})                = _length(SZ)
length{T <: FixedArray}(A::Type{T})                 = length(super(T))

endof{T,N,SZ}(A::FixedArray{T,N,SZ})                = length(A)


ndims{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})          = N
ndims{T <: FixedArray}(A::Type{T})                  = ndims(super(T))
ndims{T <: FixedArray}(A::T)                        = ndims(T)


size{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})           = _size(SZ)
size{T <: FixedArray}(A::Type{T})                   = size_or(T, ())
size{T <: FixedArray}(A::T)                         = size(T)

size{T <: FixedArray}(A::Type{T}, d::Integer)       = size(T)[d]
size{T <: FixedArray}(A::T, d::Integer)             = size(T, d)

# Iterator
start(A::FixedArray)                                = 1
next(A::FixedArray, state::Integer)                 = (A[state], state+1)
done(A::FixedArray, state::Integer)                 = length(A) < state




