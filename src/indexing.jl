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


# Unified slicing along combinations of fixed and variable dimensions

"Get index of a field with `name` in type `T`.  Should this be in Base?"
@generated function fieldindex{T}(::Type{T}, name::Symbol)
    # Expand all field names inline to allow this to become a constant
    # expression.  Simple alternative is `findfirst(fieldnames(T), name)`
    exprs = [:($(Expr(:quote, n)) == name && return $i) for (i,n) in enumerate(fieldnames(T))]
    quote
        $(Expr(:meta, :inline))
        $(exprs...)
        error("No field \"$name\" in type $T")
    end
end


# Turn `A[1,2,...]` into `destructure(A)[1,2,3,4]`
#
# Turn `A[:x,:y,1,2,...]` into
# `destructure(A)[fieldindex(eltype(A),:x), fieldindex(eltype(A),:y), 1,2,...]`
#
# Also works on the left side of an assignment
function fixed_slice_expr(expr)
    assignrhs = nothing
    if expr.head == :(=)
        @assert length(expr.args) == 2
        assignrhs = expr.args[2]
        expr = expr.args[1]
    end
    expr.head == :ref || error("Array reference not found in expression $expr")
    inds = Any[]
    name = expr.args[1]
    for i = 2:length(expr.args)
        ind = expr.args[i]
        if isa(ind,Expr) && ind.head == :quote
            push!(inds, :(fieldindex(eltype(tmp), $(esc(ind)))))
        else
            push!(inds, esc(ind))
        end
    end
    if assignrhs === nothing
        quote
            tmp = $(esc(name))
            destructure(tmp)[$(inds...)]
        end
    else
        quote
            tmp = $(esc(name))
            destructure(tmp)[$(inds...)] = $(esc(assignrhs))
        end
    end
end

"""
Slice across both fixed and variable size dimensions of
`Array{F<:FixedArray,M}`.  Before slicing, the array is reshaped to the natural
memory ordering as in `destructure()`: the `N=ndims(F)` fixed dimensions come
first, followed by the `M` dimensions of variable length.

Examples:

    # Array of fixed size vectors
    a = [Vec(i,j) for i=1:5, j=1:5]
    @fslice a[2,:,:] = 10
    xcomps = @fslice a[1,:,:]

    # Vector of fixed size matrices
    m = [@fsa([i 0; 0 i^2]) for i=1:4]
    @fslice m[1,1,:]

    # Slice immutables by field name
    immutable MyVec <: FixedVectorNoTuple{2,Float64}
        x::Float64
        y::Float64
    end
    v = MyVec[MyVec(i,-i) for i=1:5]
    @fslice v[:x, :]
"""
macro fslice(expr)
    fixed_slice_expr(expr)
end
