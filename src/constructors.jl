
immutable RandFunc{T} <: Func{1} 
    range::Range{T}
end
call{T}(rf::RandFunc{T}, x) = rand(rf.range)


immutable ConstFunctor{T} <: Base.Func{1}
    args::T
end
Base.call(f::ConstFunctor, i) = f.args
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


#TODO: making this look not like a total hack would be nice!

call{N, T}(::Type{NTuple{N, T}}, a::Vararg{T}) = a

convert{T <: Tuple}(::Type{T}, x::Real) =  ntuple(ConstFunctor(eltype(T)(x)), Val{length(T.parameters),})



call{FSA <: FixedVectorNoTuple}(::Type{FSA}, a...) = FSA(a...)

call{N, T}(::Type{FixedVector{N, T}}, a::Real, b::Real, c::Real...) = FSA(tuple(T(a),T(b), map(T, c)...))

call{FSA <: FixedVector}(::Type{FSA}, a::Real, b::Real, c::Real...) = FSA(promote(a,b, c...))
call{FSA <: Mat}(::Type{FSA}, a::Tuple, b::Tuple, c::Tuple...) = FSA((a,b,c...))

convert{FSA<: FixedVector}(T::Type{FSA}, a::AbstractVector) = FSA(tuple(a...))


call(::Type{Mat}, a::AbstractMatrix) = Mat(ntuple(x->ntuple(y->a[y,x], size(a,1)), size(a,2)))


call{R,C,T}(::Type{Mat{R,C,T}}, a::AbstractArray) = Mat(ntuple(x->ntuple(y->a[sub2ind((R,C), y, x)], C), R))


immutable ConversionIndexFunctor{T, T2} <: Func{1}
    args1::T
    target::Type{T2}
end
call(f::ConversionIndexFunctor, i) = f.target(f.args1[i])

convert{T <: Tuple}(a::Type{T}, b::FixedArray) =
    ntuple(ConversionIndexFunctor(b, eltype(T)), Val{length(T.parameters)})

convert{N, T}(a::Type{NTuple{N, T}}, b::FixedArray) =
    ntuple(IndexFunctor(b), Val{length(T.parameters)})

convert{T <: FixedVector}(a::Type{T}, b::FixedVector, x::Real...) =
    T(ntuple(ConversionIndexFunctor((b..., x...), eltype(b)), Val{length(b)+length(x)}))


convert{T <: FixedVectorNoTuple}(a::Type{T}, b::T) = b

convert{R, C, T}(a::Type{Mat{R, C, T}}, b::FixedVector) =
    map(IndexFunctor(b), a)

function convert{N, T}(a::Type{FixedVector{N, T}}, b::FixedVector)
    map(IndexFunctor(b), a)
end
#call{Row, Column, T}(::Type{FixedMatrix{Row, Column, T}}, a::Real) = FSA(ntuple(x->ntuple(y->a, Column), Row))

