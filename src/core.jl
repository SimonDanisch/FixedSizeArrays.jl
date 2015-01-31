importall Base
import Base.Func

# Alot of workarounds for not having triangular dispatch
const TYPE_PARAM_POSITION = 1
const NDIM_PARAM_POSITION = 2
const SIZE_PARAM_POSITION = 3

abstract AbstractFixedArray{T, NDIM, SIZE}

abstract AbstractFixedVector{T, CARDINALITY} <: AbstractFixedArray{T, 1, (CARDINALITY,)}
abstract AbstractFixedMatrix{T, M, N} 		 <: AbstractFixedArray{T, 2, (M, N)}

# Wrapper type, for types that have an immutable as an
abstract WrappedFixedArray{T <: AbstractFixedArray} 			 


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


getindex        (A::AbstractFixedArray,           i::Real)   			                    = A.(i)
getindex{T,SZ}  (A::AbstractFixedArray{T, 2, SZ}, i::Real, j::Real)                         = A.(sub2ind(SZ, i, j))
getindex{T,SZ}  (A::AbstractFixedArray{T, 3, SZ}, i::Real, j::Real, k::Real)                = A.(sub2ind(SZ, i, j, k))
getindex{T,N,SZ}(A::AbstractFixedArray{T, N, SZ}, i::Real, j::Real, k::Real, inds::Real...) = A.(sub2ind(SZ, i, j, k, inds...))

getindex        (A::WrappedFixedArray,            inds::Real...)                            = A.(1)[inds...]

immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = getfield(f.args1, i) 
getindex(A::AbstractFixedArray, i::AbstractFixedArray) = map(IndexFunctor(A), i)





function show{T <: AbstractFixedVector}(io::IO, F::T)
    print(io, T, "[")
    print(io, foldl(F) do v0, elem
          string(v0)*" "*string(elem)
    end)
    println(io, "]")
end
#=function show{T <: AbstractFixedMatrix}(io::IO, F::T)
    println(io, T, "[")
    for i=1:size(F, 1)
        tmp = row(F, i)
        print(io, foldl(tmp) do v0, elem
          string(v0)*" "*string(elem)
        end)
        println(io, "")
    end
    println(io, "]")
end=#

#=

stagedfunction call{T <: FixedSizeVector, N, ET}(t::Type{T}, data::ET...)
    Tsuper, Nsuper = super(T).parameters
    @assert Nsuper == N "not the right dimension"
    typename = gen_fixedsizevector_type(T, Tsuper.name, N)
    :($typename(data...))
end
gen_fixedsizevector_type(name::DataType, T::Symbol, N::Int) = gen_fixedsizevector_type(symbol(string(name.name.name)), T, N)
function gen_fixedsizevector_type(name::Symbol, T::Symbol, N::Int)
    fields = [Expr(:(::), symbol("I_$i"), T) for i = 1:N]
    typename = symbol("FS$name")
    eval(quote
        immutable $(typename){$T} <: $(name){$T}
            $(fields...)
        end
    end)
    typename
end

Base.call{FS <: AbstractFixedSizeArray, T, N}(::Type{FS}, a::Array{T, N}) = AbstractFixedSizeArray{T, N, size(a)}(a...) 
Base.call{FS <: AbstractFixedSizeArray, T, N}(::Type{FS}, a::Array{T, N}) = AbstractFixedSizeArray{T, N, size(a)}(a...) 



nvec{T, N}(x::Array{T,N}) = AbstractFixedSizeArray(x)
=#