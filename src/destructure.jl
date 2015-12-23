import Base.linearindexing
import Base.similar

"Destructuring view around an arbitrary abstract array.  See destructure()"
immutable DestructuredArray{T,N,ArrayT} <: AbstractArray{T,N}
    base::ArrayT
end

function DestructuredArray{ArrayT<:AbstractArray}(A::ArrayT)
    DestructuredArray{eltype(eltype(A)), ndims(A)+ndims(eltype(A)), ArrayT}(A)
end

linearindexing(::Type{DestructuredArray}) = Base.LinearSlow()
size{T,N,ArrayT}(A::DestructuredArray{T,N,ArrayT}) = (size(eltype(ArrayT))..., size(A.base)...)

# TODO: Should define similar() but it doesn't work for sparse matrices in 0.4
# since they don't support more than two dimensions.
# similar{S,D}(A::DestructuredArray, ::Type{S}, dims::NTuple{D,Int}) = similar(A.base, S, dims...)

@generated function getindex{T,N,ArrayT}(A::DestructuredArray{T,N,ArrayT}, inds::Int...)
    eldims = ndims(eltype(ArrayT))
    # We'd like to just do
    # A.base[inds[eldims+1:end]...][inds[1:eldims]...]
    # but splatting is currently (2015-12) a performance disaster
    fixinds = [:(inds[$i]) for i in 1:eldims]
    baseinds = [:(inds[$i]) for i in eldims+1:length(inds)]
    quote
        A.base[$(baseinds...)][$(fixinds...)]
    end
end

@generated function setindex!{T,N,ArrayT}(A::DestructuredArray{T,N,ArrayT}, value, inds::Int...)
    eldims = ndims(eltype(ArrayT))
    fixinds = [:(inds[$i]) for i in 1:eldims]
    baseinds = [:(inds[$i]) for i in eldims+1:length(inds)]
    quote
        # TODO: FixedSizeArrays.setindex gives a rather severe performance
        # penalty here.
        A.base[$(baseinds...)] = setindex(A.base[$(baseinds...)], value, $(fixinds...))
    end
end


"""
Destructure the elements of an array `A` to create `M=ndims(eltype(A))`
additional dimensions prepended to the dimensions of `A`.  The returned array
is a view onto the original elements; additional dimensions occur first
for consistency with the natural memory ordering.

For example, `AbstractArray{F<:FixedArray{T,M,SIZE},N}` appears as an
`AbstractArray{T,M+N}` after destructuring.
"""
destructure(A::AbstractArray) = DestructuredArray(A)

# Faster reinterpret-based version (as of 2015-12) for plain Array
destructure{T,N}(A::Array{T,N}) = reinterpret(eltype(T), A, (size(T)..., size(A)...))

