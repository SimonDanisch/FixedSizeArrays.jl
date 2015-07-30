immutable CArray{T, NDim, ElType} <: DenseArray{T, NDim}
	data::Array{ElType, NDim}
end

length{T,N}(A::CArray{T,N})           = div(length(A.data), length(T)) 
endof{T,N}(A::CArray{T,N})            = length(A)
size{T,N}(A::CArray{T,N})             = ntuple(i->size(A,i), N)
size{T,N}(A::CArray{T,N}, d::Integer) = div(size(A.data, d), length(T))
pointer{T,N}(A::CArray{T,N})          = Ptr{T}(pointer(A.data))

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
# This still doesn't work -.-
function Base.typed_vcat{T <: FixedArray}(:Type{T}, V1::Real,  Rest::Real...)
	println(T)
	@assert length(Rest+1) % length(T) == 0 "cannot create array, as the elements can't be fitted into a FixedSizeVector. Length: $(length(Rest)+1), length FixedArray = $(length(T))"
	@assert length(Rest+1) >= length(T) "Not enough elements given. Length: $(length(Rest)+1), needs at least = $(length(T))"
	len   	 	= length(T)
	result 		= Array(T, div(length(Rest+1), len))
	result[1] 	= T(V1, Rest[1:len-1])

	for i=len:len:length(Rest)
		result[div(i, len)] = T(Rest[i:i+len]...)
	end
	result
end
=#
function show{FSA <: FixedVector}(io::IO, a::Vector{FSA})
	print(io, "Vector with: ", length(a), "x", FSA, "[")
	for elem in a
		print(io, "[")
		print(io, join(elem, ", "))
		print(io, "]")
	end
	println(io, "]")
end

function (.+){T<:FixedArray, ND}(a::Array{T, ND}, x::T)
	T[elem + x for elem in a]
end
