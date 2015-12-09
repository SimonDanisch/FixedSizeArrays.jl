@inline getindex{T <: FixedVector}(x::T, i::Union{Range, Integer}) = x.(1)[i]
@inline getindex{T <: FixedVectorNoTuple}(x::T, i::Integer) = x.(i)
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Range, j::Int) = ntuple(IndexFunc(a, j), Val{length(i)})::NTuple{length(i), T}
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Int, j::Union{Range, Int}) = a.(1)[j][i]
@inline getindex{N, M, T}(a::Mat{N, M, T}, i::Int) = a[ind2sub((N,M), i)...]
@inline getindex(A::FixedArray, I::Tuple) = map(IndexFunctor(A), I)

@inline setindex(a::FixedArray, value, index::Int...) = map(SetindexFunctor(a, value, index), typeof(a))

@inline column{N, T}(v::FixedVector{N, T}) = v
@inline column{R, C, T}(a::Mat{R, C, T}, i::Union{Range, Int}) = a.(1)[i]

@inline row{N, T}(v::FixedVector{N, T}) = Mat{1, N, T}(v.(1))
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

