immutable CArray{T, NDim, ElType} <: DenseArray{T, NDim}
	data::Array{ElType, NDim}
end

length{T,N}(A::CArray{T,N}) 			= div(length(A.data), length(T)) 
endof{T,N}(A::CArray{T,N})  			= length(A)
size{T,N}(A::CArray{T,N})   			= ntuple(i->size(A,i), N)
size{T,N}(A::CArray{T,N}, d::Integer) 	= div(size(A.data, d), length(T))
pointer{T,N}(A::CArray{T,N}) 			= Ptr{T}(pointer(A.data))

getindex{T,N}(A::CArray{T,N}, i::Integer) = unsafe_load(Ptr{T}(pointer(A.data)), i)

reinterpret{El, N, T}(::Type{T}, A::CArray{El, N, T}) 	= A.data
reinterpret{T}(::Type{T}, A::CArray) 					= reinterpret(T, A.data)

function Base.vect{T <: MutableFixedArray}(V::T...)
	ElType = eltype(T)
	isempty(V) && return CArray{T,1, ElType}()
	dims = (length(V),length(T))
	result = Array(ElType, prod(dims))
	i = 1
    for v in V, j=1:length(T)
    	result[i] = v[j]
    	i+=1
    end
    CArray{T, 1, ElType}(result)
end

#=
function Base.typed_vcat{FSA <: FixedArray}(::Type{FSA}, X::Real...)
	ET = eltype(FSA)
	l  = length(FSA)
	length(X) % l && throw(DimensionMismatch("Number of elements in array doesn't match FixedArray type. FixedArray: $l, elements: $(length(X))"))
	[FSA(X[i:i+l]...) for i=1:l:length(X)]
end
=#