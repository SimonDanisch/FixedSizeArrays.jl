
immutable MersenneFunctor{T} <: Functor{1}
    mt::MersenneTwister
end
@compat( @inline (rf::MersenneFunctor{T}){T}(i...) = rand(rf.mt, T))

immutable RandFunctor{T} <: Functor{1}
    valuerange::T
end
@compat( @inline (rf::RandFunctor)(i...) = rand(rf.valuerange) )
immutable RandnFunctor{T} <: Functor{1}
    mt::MersenneTwister
end
@compat( @inline (rf::RandnFunctor{Float64})(i...) = randn(rf.mt) )
@compat( @inline (rf::RandnFunctor{Complex{Float64}})(i...) = randn(rf.mt) + im*randn(rf.mt) )

immutable ConstFunctor{T} <: Functor{1}
    args::T
end
@compat (f::ConstFunctor)(i...) = f.args

immutable EyeFunc{T} <: Functor{1} end
@compat( @inline (::Type{EyeFunc{T}}){T}(i::Int, j::Int) = (i==j ? one(T) : zero(T)) )

immutable UnitFunctor <: Functor{1}
    i::Int
    eltype::DataType
end
@compat (ef::UnitFunctor)(i) = ef.i==i ? one(ef.eltype) : zero(ef.eltype)

immutable ConversionIndexFunctor{T, T2} <: Functor{1}
    args1::T
    target::Type{T2}
end
@compat (f::ConversionIndexFunctor)(i...) = f.target(f.args1[i...])

immutable IndexFunctor{T} <: Functor{1}
    args1::T
end
@compat (f::IndexFunctor)(i) = f.args1[i]

immutable IndexFunc{T} <: Functor{1}
    arg::T
    i::Int
end
@compat (a::IndexFunc{T}){T}(j) = a.arg[a.i,j]

immutable CIndexFunc{T} <: Functor{1}
    arg::T
    i::Int
end
@compat (a::CIndexFunc{T}){T}(j) = a.arg[a.i,j]'


immutable CRowFunctor{M}
    mat::M
end
@compat (r::CRowFunctor)(i::Int) = crow(r.mat, i)

immutable SetindexFunctor{T <: FixedArray, V, N} <: Functor{1}
    target::T
    value::V
    index::NTuple{N, Int}
end

@compat function (sf::SetindexFunctor)(i::Int...)
    sf.index == i && return eltype(sf.target)(sf.value)
    sf.target[i...]
end

immutable RowFunctor{M}
    mat::M
end
@compat (r::RowFunctor)(i::Int) = row(r.mat, i)
