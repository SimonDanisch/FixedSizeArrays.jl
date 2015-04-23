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
function Base.show{FSA <: FixedVector}(io::IO, a::Vector{FSA})
	print(io, FSA, "[")
	for elem in a
		print(io, "[")
		for i=1:length(elem)
			print(io, elem[i], i < length(elem) ? ", " : "")
		end
		print(io, "]")
	end
	println(io, "]")
end
