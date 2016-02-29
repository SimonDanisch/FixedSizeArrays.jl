immutable MersenneFunctor{T} <: Func{1}
    mt::MersenneTwister
end
@inline call{T}(rf::MersenneFunctor{T}, i...) = rand(rf.mt, T)

immutable RandFunctor{T} <: Func{1}
    valuerange::T
end
@inline call(rf::RandFunctor, i...) = rand(rf.valuerange)
immutable RandnFunctor{T} <: Func{1}
    mt::MersenneTwister
end
@inline call(rf::RandnFunctor{Float64}, i...) = randn(rf.mt)
@inline call(rf::RandnFunctor{Complex{Float64}}, i...) = randn(rf.mt) + im*randn(rf.mt)

immutable ConstFunctor{T} <: Base.Func{1}
    args::T
end
call(f::ConstFunctor, i...) = f.args

immutable EyeFunc{T} <: Func{1} end
@inline call{T}(ef::Type{EyeFunc{T}}, i::Int, j::Int) = (i==j ? one(T) : zero(T))

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

immutable IndexFunctor{T} <: Func{1}
    args1::T
end
call(f::IndexFunctor, i) = f.args1[i]

immutable IndexFunc{T} <: Base.Func{1}
    arg::T
    i::Int
end
Base.call{T}(a::IndexFunc{T}, j) = a.arg[a.i,j]


immutable CRowFunctor{M}
    mat::M
end
call(r::CRowFunctor, i::Int) = crow(r.mat, i)

immutable SetindexFunctor{T <: FixedArray, V, N} <: Func{1}
    target::T
    value::V
    index::NTuple{N, Int}
end

function call(sf::SetindexFunctor, i::Int...)
    sf.index == i && return eltype(sf.target)(sf.value)
    sf.target[i...]
end

immutable RowFunctor{M}
    mat::M
end
call(r::RowFunctor, i::Int) = row(r.mat, i)
