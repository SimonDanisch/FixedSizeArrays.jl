eltype{T}(A::FixedArrayWrapper{T}) 					= eltype(T)

length{T}(A::FixedArrayWrapper{T})            		= length(T)

endof{T}(A::FixedArrayWrapper{T})                	= length(A)

ndims{T}(A::FixedArrayWrapper{T})            		= ndims(T)

size{T}(A::FixedArrayWrapper{T})             		= size(T)
size{T}(A::FixedArrayWrapper{T}, d::Integer) 		= size(T, d)


