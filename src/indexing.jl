getindex{T,N,SZ}(A::FixedArray{T, N, SZ}, inds::Real...) = A.(sub2ind(SZ, inds...))

function getindex{T, SZ}(A::FixedArray{T, 2, SZ}, i::Real, j::UnitRange)
    FixedArray{T, 1, (1, length(j))}(
        [A[i, k] for k in j]...
    )
end
function getindex{T, SZ}(A::FixedArray{T, 2, SZ}, j::UnitRange, i::Real)
    FixedArray{T, 1, (length(j), 1)}(
        [A[k, i] for k in j]...
    )
end
immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = getfield(f.args1, i) 
getindex(A::FixedArray, I::FixedArray) = map(IndexFunctor(A), I)

#Wrapper 
getindex{T,N,SZ}(A::FixedArray{T, N, SZ}, inds...)  = A.(1)[inds...]

