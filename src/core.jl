importall Base
import Base.Func

# Alot of workarounds for not having triangular dispatch
const TYPE_PARAM_POSITION = 1
const NDIM_PARAM_POSITION = 2
const SIZE_PARAM_POSITION = 3

abstract AbstractFixedArray{T, NDim, SIZE} <: DenseArray{T, NDim}

typealias AbstractFixedVector{T, CARDINALITY} AbstractFixedArray{T, 1, (CARDINALITY,)}
typealias AbstractFixedMatrix{T, M, N} 		  AbstractFixedArray{T, 2, (M, N)}

	 


eltype{T,N,SZ}(A::AbstractFixedArray{T,N,SZ}) 				= T
eltype{T <: AbstractFixedArray}(A::Type{T})                 = first(T.types) 

length{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})           	= prod(SZ)
length{T <: AbstractFixedArray}(A::Type{T})                 = prod(super(super(T)).parameters[SIZE_PARAM_POSITION])

endof{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})                = length(A)


ndims{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})            	= N
ndims{T <: AbstractFixedArray}(A::Type{T})            		= super(T).parameters[NDIM_PARAM_POSITION]

size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})             	= SZ
size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ}, d::Integer) 	= SZ[d]

size{T <: AbstractFixedArray}(A::Type{T})            		= super(T).parameters[SIZE_PARAM_POSITION]
size{T <: AbstractFixedArray}(A::Type{T}, d::Integer) 		= super(T).parameters[SIZE_PARAM_POSITION][d]

# Iterator 
start(A::AbstractFixedArray)            					= 1
next (A::AbstractFixedArray, state::Integer) 				= (A[state], state+1)
done (A::AbstractFixedArray, state::Integer) 				= length(A) < state


getindex{T,N,SZ}(A::AbstractFixedArray{T, N, SZ}, inds::Real...) = A.(sub2ind(SZ, inds...))


immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = getfield(f.args1, i) 
getindex(A::AbstractFixedArray, i::AbstractFixedArray) = map(IndexFunctor(A), i)



#Constructor:

tuple_to_string(t::(), sep) = ""
tuple_to_string(t::(Any,), sep) = "$(t[1])"
tuple_to_string(t::(Any, Any...), sep) = "$(t[1])$sep" * tuple_to_string(Base.tail(t), sep)
vec_name(sz::(Integer, Integer...)) = symbol("NVec" * tuple_to_string(sz, 'x'))
vec_name(sz::(Integer,)) = symbol("NVec" * string(first(sz)))

gen_fixedsizevector_type(name::DataType, T::Symbol, N::Int) = gen_fixedsizevector_type(symbol(string(name.name.name)), T, N)
function gen_fixedsizevector_type(SIZE::(Integer...))
    fields      = [Expr(:(::), symbol("I_$i"), :T) for i = 1:prod(SIZE)]
    NDim        = length(SIZE)
    typename    = vec_name(SIZE)
    eval(quote
        immutable $(typename){T} <: AbstractFixedArray{T, $NDim, $SIZE}
            $(fields...)
        end
    end)
    typename
end

stagedfunction call{T, NDim, SIZE}(t::Type{AbstractFixedArray{T, NDim, SIZE}}, data::T...)
    N = length(data)
    @assert prod(SIZE) == N "not the right dimension"
    typename = gen_fixedsizevector_type(SIZE)
    :($typename(data...))
end
stagedfunction call{T, NDim, SIZE}(t::Type{AbstractFixedArray{T, NDim, SIZE}}, data::(T...,))
    N = length(data)
    @assert prod(SIZE) == N "not the right dimension"
    typename = gen_fixedsizevector_type(SIZE)
    :($typename(data...))
end

Base.call{FS <: AbstractFixedArray, T, N}(::Type{FS}, a::Array{T, N}) = AbstractFixedArray{T, N, size(a)}(a...) 
nvec{T, N}(x::Array{T,N})             = AbstractFixedArray(x)
nvec{T}(x::T...)                      = AbstractFixedArray{T, 1, (length(x),)}(x)
nvec{T}(SIZE::(Integer...,), x::T...) = AbstractFixedArray{T, length(SIZE), SIZE}(x)

@show nvec((2,2,2), 1f0,1f0,1f0,1f0, 1f0,1f0,1f0,1f0)
