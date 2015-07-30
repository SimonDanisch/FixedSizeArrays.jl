
immutable RandFunc{T} <: Func{1} 
    range::Range{T}
end
call{T}(rf::RandFunc{T}, x) = rand(rf.range)

immutable ConstFunctor{T} <: Func{1}
    args::T
end
call(f::ConstFunctor, i) = f.args

immutable EyeFunc{N} <: Func{1}
    size::NTuple{N, Int}
    eltype::DataType
end
function call{T}(ef::EyeFunc{T}, x)
    i,j = ind2sub(ef.size, x)
    i==j ? one(ef.eltype) : zero(ef.eltype)
end

immutable UnitFunctor <: Func{1}
    i::Int
    eltype::DataType
end
call(ef::UnitFunctor, x) = ef.i==x ? one(ef.eltype) : zero(ef.eltype)


zero{FSA <: FixedArray}(::Type{FSA}) = map(ConstFunctor(zero(eltype(FSA))), FSA)
one{FSA <: FixedArray}(::Type{FSA})  = map(ConstFunctor(one(eltype(FSA))), FSA)
eye{FSA <: FixedArray}(::Type{FSA})  = map(EyeFunc(size(FSA), eltype(FSA)), FSA)
unit{FSA <: FixedVector}(::Type{FSA}, i::Integer) = map(UnitFunctor(i, eltype(FSA)), FSA)
function rand{FSA <: FixedArray}(x::Type{FSA})
    T = eltype(FSA)
    T <: applicable(eps, T) && return map(RandFunc(zero(T) : eps(T) : one(T)), FSA) # this case is basically for FixedPointNumbers
    map(RandFunc(typemin(T) : typemax(T)), FSA)
end
rand{FSA <: FixedArray}(x::Type{FSA}, range::Range) = map(RandFunc(range), FSA)


call{FSA <: FixedArray, T}(::Type{FSA}, a::T, b::T, c::T...) = FSA(tuple(a, b, c...))
call{FSA <: FixedVector}(::Type{FSA}, a::Real) = FSA(ntuple(x->a, length(FSA)))

#call{Row, Column, T}(::Type{FixedMatrix{Row, Column, T}}, a::Real) = FSA(ntuple(x->ntuple(y->a, Column), Row))

