"Destructuring view around an arbitrary abstract array.  See destructure()"
immutable DestructuredArray{T,N,ArrayT} <: AbstractArray{T,N}
    base::ArrayT
end

import Base.linearindexing

linearindexing(::Type{DestructuredArray}) = Base.LinearSlow()
size{T,N,ArrayT}(a::DestructuredArray{T,N,ArrayT}) = (size(eltype(ArrayT))..., size(a.base)...)

@generated function getindex{T,N,ArrayT}(a::DestructuredArray{T,N,ArrayT}, inds::Int...)
    eldims = ndims(eltype(ArrayT))
    # We'd like to just do
    # a.base[inds[eldims+1:end]...][inds[1:eldims]...]
    # but splatting is currently (2015-12) a performance disaster
    fixinds = [:(inds[$i]) for i in 1:eldims]
    baseinds = [:(inds[$i]) for i in eldims+1:length(inds)]
    quote
        a.base[$(baseinds...)][$(fixinds...)]
    end
end

@generated function setindex!{T,N,ArrayT}(a::DestructuredArray{T,N,ArrayT}, value, inds::Int...)
    eldims = ndims(eltype(ArrayT))
    fixinds = [:(inds[$i]) for i in 1:eldims]
    baseinds = [:(inds[$i]) for i in eldims+1:length(inds)]
    quote
        # TODO: FixedSizeArrays.setindex gives a rather severe performance
        # penalty here.
        a.base[$(baseinds...)] = setindex(a.base[$(baseinds...)], value, $(fixinds...))
    end
end


"""
Destructure the elements of an array `a` to create `M=ndims(eltype(a))`
additional dimensions prefixing the dimensions of `a`.  The returned array is a
view onto the original elements, with fixed dimensions first according to the
natural memory ordering.

For example, `AbstractArray{F<:FixedArray{T,M,SIZE},N}` appears as as an
`AbstractArray{T,M+N}`.
"""
function destructure{ArrayT<:AbstractArray}(a::ArrayT)
    DestructuredArray{eltype(eltype(a)), ndims(a)+ndims(eltype(a)), ArrayT}(a)
end

# Faster (as of 2015-12) reinterpret-based version for plain Array
destructure{T,N}(a::Array{T,N}) = reinterpret(eltype(T), a, (size(T)..., size(a)...))

