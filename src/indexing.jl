@inline getindex{T,SZ}(A::FixedArray{T, 1, SZ}, i::Int) = A.(i)
@inline getindex{T,N,SZ}(A::FixedArray{T, N, SZ}, inds::Real...) = A.(sub2ind(size(A)::NTuple{N, Int}, inds...))

function getindex{T, SZ}(A::FixedArray{T, 2, SZ}, i::Real, j::UnitRange)
    FixedVector{T, length(j)}(
        [A[i, k] for k in j]...
    )
end
function getindex{T, SZ}(A::FixedArray{T, 2, SZ}, j::UnitRange, i::Real)
    FixedVector{T, length(j)}(
        [A[k, i] for k in j]...
    )
end

immutable IndexFunctorTuple{T, T2} <: Func{1}
    indexes::T
    target::T2
end
call(f::IndexFunctorTuple, i) = f.target[f.indexes[i]]
function getindex{T <: FixedArray, TuPl <: Tuple}(A::T, I::Type{TuPl})
    range = TuPl.parameters[1] #Of form Tuple{1:3}
    map(IndexFunctorTuple(range, A), FixedVector{eltype(A), length(range)})
end
immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = f.args1[i] 
getindex(A::FixedArray, I::FixedArray) = map(IndexFunctor(A), I)

#Wrapper 

@generated function row{T, Column, Row}(A::FixedMatrix{T, Row, Column}, i::Integer)
    fields = [:(A[i,$j]) for j=1:Column]
    returntype = gen_fixedsizevector_type((Column,), A.mutable)
    :($returntype($(fields...)))
end
@generated function column{T, Column, Row}(A::FixedMatrix{T, Row, Column}, j::Integer)
    fields = [:(A[$i,j]) for i=1:Row]
    returntype = gen_fixedsizevector_type((Row,), A.mutable)
    :($returntype($(fields...)))
end

@generated function row{T, Column, Row, VT}(A::FixedMatrix{T, Row, Column}, i::Integer, t::Type{VT})
    fields = [:(A[i,$j]) for j=1:Column]
    :(VT($(fields...)))
end
@generated function column{T, Column, Row, VT}(A::FixedMatrix{T, Row, Column}, j::Integer, t::Type{VT})
    fields = [:(A[$i,j]) for i=1:Row]
    :(VT($(fields...)))
end


function row{T, N}(v::FixedVector{T, N})
    convert(FixedMatrix{T, 1, N}, v)
end
column{T, N}(v::FixedVector{T, N}) = v