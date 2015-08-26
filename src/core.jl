abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, Tuple{CARDINALITY}}
typealias MutableFixedMatrix{T, M, N} 		 MutableFixedArray{T, 2, Tuple{M,N}}

typealias FixedVector{CARDINALITY, T} FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{Row, Column, T}        FixedArray{T, 2, Tuple{Row, Column}}

abstract FixedVectorNoTuple{CARDINALITY, T} <: FixedVector{CARDINALITY, T}
export FixedVectorNoTuple


# Get the abstract FixedSizeArray type, even for complex type hirarchies
function fixedsizearray_type{FSA <: FixedArray}(::Type{FSA})
    ff = FSA
    while ff.name.name != :FixedArray
        ff = super(ff)
        if ff == Any
            error("Uncommon type hierarchy encountered. Please report issue on Github")
        end
    end
    ff
end
isfullyparametrized{T}(::Type{T}) = !any(x-> isa(x, TypeVar), T.parameters)


eltype{T,N,SZ}(A::FixedArray{T,N,SZ}) 				= T
eltype{T,N,SZ}(A::Type{FixedArray{T,N,SZ}}) 		= T
eltype{T <: FixedArray}(A::Type{T})                 = eltype(super(T))

length{T, L}(A::FixedVector{L, T})           		= L

length{T, M, N}(A::FixedMatrix{M, N, T})           	= M*N
length{T,N,SZ}(A::FixedArray{T,N,SZ})           	= prod(SZ.parameters)::Int
length{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})         = prod(SZ.parameters)::Int
#This is soo bad. But a non fully parametrized abstract type doesn't get catched by the above function
length{T <: FixedArray}(A::Type{T})                 = isfullyparametrized(T) ? length(super(T)) : prod(super(A).parameters[3].parameters)


endof{T,N,SZ}(A::FixedArray{T,N,SZ})                = length(A)


ndims{T,N,SZ}(A::FixedArray{T,N,SZ})            	= N
ndims{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})          = N
ndims{T <: FixedArray}(A::Type{T})            		= ndims(super(T))

size{T,N,SZ}(A::FixedArray{T,N,SZ})             	= (SZ.parameters...)::NTuple{N, Int}
size{T,N,SZ}(A::FixedArray{T,N,SZ}, d::Integer) 	= SZ.parameters[d]::Int

size{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})           = (SZ.parameters...)::NTuple{N, Int}
size{T <: FixedArray}(A::Type{T})            		= size(super(T))

size{T <: FixedArray}(A::Type{T}, d::Integer) 		= size(T)[d]::Int

# Iterator
start(A::FixedArray)            					= 1
next(A::FixedArray, state::Integer) 				= (A[state], state+1)
done(A::FixedArray, state::Integer) 				= length(A) < state


immutable Mat{Row, Column, T} <: FixedMatrix{Row, Column, T}
    _::NTuple{Column, NTuple{Row, T}}
end

function show{R,C,T}(io::IO, m::Mat{R,C,T})
	println(io, typeof(m), "(")
	for i=1:R
		println(io, "    ", join(row(m, i), " "))
	end
	println(io, ")")
end
