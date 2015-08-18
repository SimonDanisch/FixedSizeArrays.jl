type_name(fsa) = :($(symbol(fsa.name)))

function size_or{FSA <: FixedArray}(::Type{FSA}, SZ)
    fsa = fixedsizearray_type(FSA)
    s  = tuple(fsa.parameters[3].parameters...)
    any(sz -> isa(sz, TypeVar), s) ? SZ : s
end
function eltype_or{FSA <: FixedArray}(::Type{FSA}, ElType)
    fsa = fixedsizearray_type(FSA)
    s   = fsa.parameters[1]
    isa(s, TypeVar) ? ElType : s
end


@generated function call{FSA <: FixedArray, T1 <: FixedArray}(::Type{FSA}, a::T1, b...)
    SZ      = size_or(FSA, nothing)
    ElType  = eltype_or(FSA, eltype(T1))
    a_expr  = :( a )
    if SZ != nothing
        if (prod(SZ) > (length(T1) + length(b)))
            throw(DimensionMismatch("$FSA is too small, can not be constructed from array $a and $b arguments"))
        elseif prod(SZ) < length(T1) && isempty(b) #shrinking constructor, e.g. Vec3(Vec4(...))
            a_expr = :( a[1:$(prod(SZ))] )
        end
    end
    if isempty(b)
        if ElType != eltype(T1)
            return :( $FSA(map($ElType, $a_expr)...) )
        else
            return :( $FSA($a_expr...) )
        end
    else
        return :( $FSA(a..., b...) )
    end
end


@generated function call{FSA <: FixedArray, T <: Array}(::Type{FSA}, a::T)
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

call{FSA <: FixedVectorNoTuple}(::Type{FSA}, a::Tuple, b::Tuple...) = error("$FSA can't be constructed from $a")
call{FSA <: FixedVectorNoTuple}(::Type{FSA}, a::Tuple) = FSA(a...)

call{FSA <: FixedArray, T}(::Type{FSA}, a::T..., ) = FSA(a)



@generated function call{FSA <: FixedArray, X}(::Type{FSA}, a::X)
    SZ      = size_or(FSA, (1,))
    ElType  = eltype_or(FSA, a)
    Len     = prod(SZ)
    T_N     = FSA
    if a <: Tuple # a::Tuple is ambigous to default constructor, so need to do it here
        tuple_expr = any(x-> x!=ElType, a.parameters) ? :( map($ElType, a) ) : :(a)
        if FSA <: FixedVectorNoTuple
            return :( $T_N($tuple_expr...) )
        else
            return :($T_N($tuple_expr))
        end
    end
    if FSA <: FixedVectorNoTuple
        return :($T_N($(ntuple(i-> :($ElType(a)), Len)...)))
    else
        expr = fill_tuples_expr((inds...)->:($ElType(a[1])), SZ)
        return :($T_N($expr))
    end
end

@generated function call{FSA <: FixedArray}(::Type{FSA}, a...)
    SZ     = size_or(FSA, (length(a),))
    ElType = eltype_or(FSA, promote_type(a...))
    if FSA <: FixedVectorNoTuple
        length(a) != prod(SZ) && throw(DimensionMismatch("can't construct $FSA with $(length(a)) arguments. Args: $a"))
        return :($FSA(map($ElType, a)...))
    else
        all(x-> x <: Tuple, a) && return :( $FSA(a) ) # TODO be smarter about this
        any(x-> x != ElType, a) && return :($FSA(map($ElType, a)))
        return :($FSA(a))
    end
end


const_fill(T, sym, SZ, inds...)        = :($T($sym[1]))
flat_fill(T, sym, SZ, inds...)         = :($T($sym[$(sub2ind(SZ, inds...))]))
array_fill(T, sym, SZ, inds...)        = :($T($sym[1][$(inds...)]))
hierarchical_fill(T, sym, SZ, inds...) = :($T($sym[$(inds...)]))


_fill_tuples_expr(inner::Function, SZ::Tuple{Int}, inds...) =
    :(tuple($(ntuple(i->inner(i, inds...), SZ[1])...)))
_fill_tuples_expr{N}(inner::Function, SZ::NTuple{N, Int}, inds...) =
    :(tuple($(ntuple(i->_fill_tuples_expr(inner, SZ[1:end-1], i, inds...),SZ[end])...)))
fill_tuples_expr(inner::Function, SZ::Tuple) = _fill_tuples_expr(inner, SZ)


_fill_tuples(inner, originalSZ, SZ::Tuple{Int}, inds::Int...) =
    ntuple(i->inner(SZ, i, inds...), Val{SZ[1]})
_fill_tuples{N}(inner, originalSZ, SZ::NTuple{N, Int}, inds::Int...) =
    ntuple(i->_fill_tuples(inner, originalSZ, SZ[1:end-1], i, inds...), Val{SZ[end]})
fill_tuples{N}(inner, SZ::NTuple{N, Int}) = _fill_tuples(inner, SZ, SZ)




zero{FSA <: FixedArray}(::Type{FSA}) = map(ConstFunctor(zero(eltype(FSA))), FSA)
one{FSA <: FixedArray}(::Type{FSA})  = map(ConstFunctor(one(eltype(FSA))), FSA)
eye{FSA <: FixedArray}(::Type{FSA})  = map(EyeFunc(size(FSA), eltype(FSA)), FSA)
unit{FSA <: FixedVector}(::Type{FSA}, i::Integer) = map(UnitFunctor(i, eltype(FSA)), FSA)

function rand{FSA <: FixedArray}(x::Type{FSA})
    T = eltype(FSA)
    T <: applicable(eps, T) && return map(RandFunctor(zero(T) : eps(T) : one(T)), FSA) # this case is basically for FixedPointNumbers
    map(RandFunctor(typemin(T) : typemax(T)), FSA)
end
rand{FSA <: FixedArray}(x::Type{FSA}, range::Range) = map(RandFunctor(range), FSA)


#conversion
convert{T <: Tuple}(::Type{T}, x::Real)                     =  ntuple(ConstFunctor(eltype(T)(x)), Val{length(T.parameters),})
convert{T <: Tuple, FSA <: FixedArray}(::Type{T}, x::FSA)   = map(eltype(T), x.(1)[1:length(T.parameters)])
convert{T <: FixedArray}(t::Type{T}, f::T)                  = f
convert{FSA1 <: FixedArray}(t::Type{FSA1}, f::FixedArray) =
    map(ConversionIndexFunctor(f, eltype_or(FSA1, eltype(typeof(f)))), FSA1)
