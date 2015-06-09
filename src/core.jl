importall Base
import Base.Func

# Alot of workarounds for not having triangular dispatch
const TYPE_PARAM_POSITION = 1
const NDIM_PARAM_POSITION = 2
const SIZE_PARAM_POSITION = 3

abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, @compat(Tuple{CARDINALITY})}
typealias MutableFixedMatrix{T, M, N} 		 MutableFixedArray{T, 2, @compat(Tuple{M,N})}

typealias FixedVector{T, CARDINALITY} FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{T, M, N}        FixedArray{T, 2, Tuple{M, N}}

abstract FixedArrayWrapper{T <: FixedArray} <: FixedArray

isfullyparametrized{T}(::Type{T}) = !any(x-> isa(x, TypeVar), T.parameters)


eltype{T,N,SZ}(A::FixedArray{T,N,SZ}) 				= T
eltype{T <: FixedArray}(A::Type{T})                 = T.types[1]

length{T, L}(A::FixedVector{T,L})           		= L
length{T, M, N}(A::FixedMatrix{T,M, N})           		= M*N
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
next (A::FixedArray, state::Integer) 				= (A[state], state+1)
done (A::FixedArray, state::Integer) 				= length(A) < state


#Utilities:
name(typ::DataType) = string(typ.name.name)
fieldname(i) = symbol("i_$i")
# Function to strip of parameters from a type definition, to avoid conversion.
# eg: Point{Float32}(1) would end up as Point{Float32}(convert(Float32, 1))
# whereas Point(1) -> Point{Int}(1)
# Main is needed, as the resulting type is mostly defined in Main, so it wouldn't be found otherwise
function without_params{T}(::Type{T})
    eval(:(Main.$(symbol(name(T)))))
end
