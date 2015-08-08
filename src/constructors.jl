function size_or{FSA <: FixedArray}(::Type{FSA}, SZ)
    fsa = fixedsizearray_type(FSA)
    s  = tuple(fsa.parameters[3].parameters...)
    any(sz -> isa(sz, TypeVar), s) ? SZ : s
end
function eltype_or{FSA <: FixedArray}(::Type{FSA}, ElType)
    fsa = fixedsizearray_type(FSA)
    s  = fsa.parameters[1]
    isa(s, TypeVar) ? ElType : s
end


@generated function Base.call{FSA <: FixedArray}(::Type{FSA}, a)
    expr = fsa_constructor(FSA, a)
    expr
end
@generated function Base.call{FSA <: FixedArray}(::Type{FSA}, a...)
    fsa_constructor(FSA, a...)
end


function fsa_constructor{FSA <: FixedArray, T1 <: FixedArray}(::Type{FSA}, a::Type{T1}, b...)
    SZ = size_or(FSA, nothing)
    ElType = eltype_or(FSA, eltype(T1))
    if SZ != nothing
       (prod(SZ) != (length(T1) + length(b))) && throw(DimensionMismatch("size of $FSA is not fitting array $a + $b arguments"))
    end
    if isempty(b)
        if ElType != eltype(T1)
            return :( $FSA(map($ElType, a)...) )
        else
            return :( $FSA(a...) )
        end
    else
        return :( $FSA(a[1]..., a[2]...) )
    end
end


function fsa_constructor{FSA <: FixedArray, T <: Array}(::Type{FSA}, a::Type{T})
    SZ     = size_or(FSA, :(size(a)))
    ElType = eltype_or(FSA, eltype(a))
    expr = :($FSA(fill_tuples((sz, i...)->$ElType(a[i...]), $SZ)))
    if FSA <: FixedVectorNoTuple
        expr = :($FSA(a...))
    end
    quote
        $SZ != size(a) && throw(DimensionMismatch("size of $FSA is not fitting array $(typeof(a)), with size: $(size(a))"))
        $expr
    end
end


function fsa_constructor{FSA <: FixedArray, X}(::Type{FSA}, a::X)
    SZ      = size_or(FSA, (1,))
    ElType  = eltype_or(FSA, a)
    Len     = prod(SZ)
    if FSA <: FixedVectorNoTuple
        return :($FSA($(ntuple(i-> :($ElType(a)), Len)...)))
    else
        expr = fill_tuples_expr(const_fill, ElType, :a, SZ)
        return :($FSA($expr))
    end
end

function fsa_constructor{FSA <: FixedArray}(::Type{FSA}, a...)
    SZ     = size_or(FSA, (length(a),))
    ElType = eltype_or(FSA, promote_type(a...))
    if FSA <: FixedVectorNoTuple
        return :($FSA(map($ElType, a)...))
    else
        a[1] <: Tuple && return :($FSA(a))
        any(x-> x != ElType, a) && return :($FSA(map($ElType, a))) 
        return :($FSA(a))
    end
end





const_fill(T, sym, SZ, inds...)        = :($T($sym[1]))
flat_fill(T, sym, SZ, inds...)         = :($T($sym[$(sub2ind(SZ, inds...))]))
array_fill(T, sym, SZ, inds...)        = :($T($sym[1][$(inds...)]))
hierarchical_fill(T, sym, SZ, inds...) = :($T($sym[$(inds...)]))


_fill_tuples_expr(inner::Function, T, sym, originalSZ, SZ::Tuple{Int}, inds...) = 
    :(tuple($(ntuple(i->inner(T, sym, originalSZ, i, inds...), SZ[1])...)))
_fill_tuples_expr{N}(inner::Function, T, sym, originalSZ, SZ::NTuple{N, Int}, inds...) = 
    :(tuple($(ntuple(i->_fill_tuples_expr(inner, T, sym, originalSZ, Base.tail(SZ), i, inds...), SZ[1]))))
fill_tuples_expr(inner::Function, T, sym, SZ) = _fill_tuples_expr(inner, T, sym, SZ, SZ)


_fill_tuples(inner, originalSZ, SZ::Tuple{Int}, inds::Int...) = 
    ntuple(i->inner(SZ, i, inds...), Val{SZ[1]})
_fill_tuples{N}(inner, originalSZ, SZ::NTuple{N, Int}, inds::Int...) = 
    ntuple(i->_fill_tuples(inner, originalSZ, SZ[1:end-1], i, inds...), Val{SZ[end]})
fill_tuples{N}(inner, SZ::NTuple{N, Int}) = _fill_tuples(inner, SZ, SZ)


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


immutable ConversionIndexFunctor{T, T2} <: Func{1}
    args1::T
    target::Type{T2}
end
call(f::ConversionIndexFunctor, i) = f.target(f.args1[i])


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


#conversion
convert{T <: Tuple}(::Type{T}, x::Real) =  ntuple(ConstFunctor(eltype(T)(x)), Val{length(T.parameters),})
convert{T <: Tuple, FSA <: FixedArray}(::Type{T}, x::FSA) = map(eltype(T), x.(1)[1:length(T.parameters)])
convert{T <: FixedArray}(t::Type{T}, f::T) = f
convert{FSA1 <: FixedArray}(t::Type{FSA1}, f::FixedArray) =
    map(ConversionIndexFunctor(f, eltype_or(FSA1, eltype(typeof(f)))), FSA1)
