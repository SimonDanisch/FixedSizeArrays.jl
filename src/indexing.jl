@inline getindex{T <: FixedVector}(x::T, i::Integer) = x.(1)[i]
@inline getindex(x::FixedMatrix, i::Integer, j::Union(Integer, Range)) = x.(1)[i][j]
@inline getindex(x::FixedMatrix, i::Integer) = x[ind2sub(size(x), i)...]
    

@inline row(x::FixedMatrix, i::Integer) = x[i]

@generated function column{T, Row, Column}(A::FixedMatrix{T, Row, Column}, j::Integer)
    fields = ntuple(i-> :(A[$i][j]), Row)
    :(tuple($(fields...)))
end

immutable IndexFunctorTuple{T, T2} <: Func{1}
    indexes::T
    target::T2
end
call(f::IndexFunctorTuple, i) = f.target[f.indexes[i]]

immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = f.args1[i] 

function getindex{T <: FixedArray, TuPl <: Tuple}(A::T, I::Type{TuPl})
    range = TuPl.parameters[1] #Of form Tuple{1:3}
    map(IndexFunctorTuple(range, A), FixedVector{eltype(A), length(range)})
end

getindex(A::FixedArray, I::Union(NTuple, FixedArray)) = map(IndexFunctor(A), I)


row{T, N}(v::FixedVector{T, N}) = convert(FixedMatrix{T, 1, N}, v)
column{T, N}(v::FixedVector{T, N}) = v