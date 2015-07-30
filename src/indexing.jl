immutable IndexFunc{T} <: Base.Func{1}
    arg::T
    i::Int
end
Base.call{T}(a::IndexFunc{T}, j) = a.arg[a.i,j] 

getindex{T <: FixedVector}(x::T, i::Integer) = x.(1)[i]


getindex{N, M, T}(a::Mat{N, M, T}, i::Int, j::Union(Range, Int)) = a.(1)[j][i]
getindex{N, M, T}(a::Mat{N, M, T}, i::Int) 	 = a[ind2sub((N,M), i)...]

column{R, C, T}(a::Mat{R, C, T}, i::Union(Range, Int)) = a.(1)[i]
row{R, C, T}(a::Mat{R, C, T}, j::Int) = ntuple(IndexFunc(a, j), Val{C})::NTuple{C, T}

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


row{N, T}(v::FixedVector{N, T}) = convert(FixedMatrix{1, N, T}, v)

row{N, T}(v::FixedVector{N, T}, i::Int) = (v[i],)

column{N, T}(v::FixedVector{N, T}) = v