immutable RandFunctor{T} <: Func{1}
    range::Range{T}
end
call{T}(rf::RandFunctor{T}, i...) = rand(rf.range)

immutable ConstFunctor{T} <: Base.Func{1}
    args::T
end
call(f::ConstFunctor, i...) = f.args


immutable EyeFunc{N} <: Func{1}
    size::NTuple{N, Int}
    eltype::DataType
end
function call{T}(ef::EyeFunc{T}, i::Int, j::Int)
    i==j ? one(ef.eltype) : zero(ef.eltype)
end

immutable UnitFunctor <: Func{1}
    i::Int
    eltype::DataType
end
call(ef::UnitFunctor, i) = ef.i==i ? one(ef.eltype) : zero(ef.eltype)


immutable ConversionIndexFunctor{T, T2} <: Func{1}
    args1::T
    target::Type{T2}
end
call(f::ConversionIndexFunctor, i...) = f.target(f.args1[i...])

immutable IndexFunctorTuple{T, T2} <: Func{1}
    indexes::T
    target::T2
end
call(f::IndexFunctorTuple, i) = f.target[f.indexes[i]]

immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = f.args1[i]

immutable IndexFunc{T} <: Base.Func{1}
    arg::T
    i::Int
end
Base.call{T}(a::IndexFunc{T}, j) = a.arg[a.i,j]


immutable RowFunctor{M}
    mat::M
end
call(r::RowFunctor, i::Int) = row(r.mat, i)



immutable SetindexFunctor{T <: FixedArray, V, N} <: Func{1}
    target::T
    value::V
    index::NTuple{N, Int}
end

function call(sf::SetindexFunctor, i::Int...)
    sf.index == i && return sf.value
    sf.target[i...]
end