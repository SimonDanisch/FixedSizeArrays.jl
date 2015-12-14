@inline getindex{T <: FixedVector}(x::T, i::Union{Range, Integer}) = x.(1)[i]
@inline getindex{T <: FixedVectorNoTuple}(x::T, i::Integer) = x.(i)
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Range, j::Int) = ntuple(IndexFunc(a, j), Val{length(i)})::NTuple{length(i), T}
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Int, j::Union{Range, Int}) = a.(1)[j][i]
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Int) = a[ind2sub((N,M), i)...]
@inline getindex(A::FixedArray, I::Tuple) = map(IndexFunctor(A), I)

@inline setindex(a::FixedArray, value, index::Int...) = map(SetindexFunctor(a, value, index), typeof(a))

@inline column{N, T}(v::FixedVector{N, T}) = v
@inline column{R, C, T}(a::Mat{R, C, T}, i::Union{Range, Int}) = a.(1)[i]

@inline row{N, T}(v::FixedVector{N, T}) = Mat{1, N, T}(v...)
@inline row{N, T}(v::FixedVector{N, T}, i::Int) = (v[i],)
@inline row{R, C, T}(a::Mat{R, C, T}, j::Int) = ntuple(IndexFunc(a, j), Val{C})::NTuple{C, T}
@inline row{R, T}(a::Mat{R, 1, T}, j::Int) = (a.(1)[1][j],)
@inline row{R, T}(a::Mat{R, 2, T}, j::Int) = (a.(1)[1][j], a.(1)[2][j])
@inline row{R, T}(a::Mat{R, 3, T}, j::Int) = (a.(1)[1][j], a.(1)[2][j], a.(1)[3][j])
@inline row{R, T}(a::Mat{R, 4, T}, j::Int) = (a.(1)[1][j], a.(1)[2][j], a.(1)[3][j], a.(1)[4][j],)


# the columns of the ctranspose are the complex conjugate rows
@inline crow{R, C, T}(a::Mat{R, C, T}, j::Int) = ntuple(IndexFunc(a, j)', Val{C})::NTuple{C, T}
@inline crow{R, T}(a::Mat{R, 1, T}, j::Int) = (a.(1)[1][j]',)
@inline crow{R, T}(a::Mat{R, 2, T}, j::Int) = (a.(1)[1][j]', a.(1)[2][j]')
@inline crow{R, T}(a::Mat{R, 3, T}, j::Int) = (a.(1)[1][j]', a.(1)[2][j]', a.(1)[3][j]')
@inline crow{R, T}(a::Mat{R, 4, T}, j::Int) = (a.(1)[1][j]', a.(1)[2][j]', a.(1)[3][j]', a.(1)[4][j]',)


# Placeholder type for dispatch to slice along fixed dimensions of
# Array{F<:FixedArray}
immutable FSlice; end

"""
Remove the fixed dimensional structure from an
`Array{F<:FixedArray{T,M,SIZE},N}`, reinterpreting it as an `Array{T,M+N}`, in
natural memory ordering (ie, fixed dimensions first).
"""
function destructure{F<:FixedArray}(a::Array{F})
    SIZE = size(F)
    T = eltype(F)
    reinterpret(T, a, (SIZE..., size(a)...))
end

function getindex{F<:FixedArray}(a::Array{F}, ::Type{FSlice}, inds...)
    destructure(a)[inds...]
end

function setindex!{F<:FixedArray}(a::Array{F}, val, ::Type{FSlice}, inds...)
    destructure(a)[inds...] = val
end

function getindex{F<:FixedVectorNoTuple}(a::Array{F}, fieldname::Symbol, inds...)
    # Not terribly efficient; can we somehow force the field index to be
    # resolved at compile time?
    fieldindex = findfirst(fieldnames(F), fieldname)
    fieldindex > 0 || error("No field named $fieldname in $(typeof(a))")
    destructure(a)[fieldindex, inds...]
end

