@generated function fsa_abstract{FSA <: FixedArray}(::Type{FSA})
    ff = FSA
    while ff.name.name != :FixedArray
       ff = super(ff)
    end
    :($ff)
end
@generated function size_or{FSA <: FixedArray}(::Type{FSA}, OR)
    fsatype = fsa_abstract(FSA)
    sz = fsatype.parameters[3]
    any(x->isa(x, TypeVar), sz.parameters) && return :(OR)
    :($(_size(sz)))
end
@generated function eltype_or{FSA <: FixedArray}(::Type{FSA}, OR)
    fsatype = fsa_abstract(FSA)
    T = fsatype.parameters[1]
    isa(T, TypeVar) && return :(OR)
    :($T)
end

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
        if ElType != eltype(T1) && FSA <: T1
            return :( map($ElType, $a_expr) )
        elseif ElType == eltype(T1) && !(FSA <: T1)
            return :( $FSA($a_expr...) )
        else
            return :( $FSA($a_expr...) )
        end
    else
        return :( $FSA(a..., b...) )
    end
end

immutable ParseFunctor{T, S <: AbstractString} <: Func{1}
    t::Type{T}
    a::Vector{S}
end
call{T}(pf::ParseFunctor{T}, i::Int) = parse(T, pf.a[i])
call(pf::ParseFunctor{Void}, i::Int) = parse(pf.a[i])

@generated function call{FSA <: FixedArray, T <: Array}(::Type{FSA}, a::T)
    if eltype(a) <: AbstractString
        ElType = eltype_or(FSA, Void)
        return :(map(ParseFunctor($ElType, a), FSA)) # can't be defined in another method as it leads to lots of ambigouity with the default constructor
    end
    SZ     = size_or(FSA, :(size(a)))
    ElType = eltype_or(FSA, eltype(a))
    if isa(SZ, Expr)
        expr = :($FSA(fill_tuples((sz, i...)->$ElType(a[i...]), $SZ)))
    else
        tupexpr = fill_tuples_expr((i,inds...) -> :($ElType(a[$i, $(inds...)])), SZ)
        expr = :($FSA($tupexpr))
    end
    if FSA <: FixedVectorNoTuple
        expr = :($FSA(a...))
    end
    quote
        $SZ != size(a) && throw(DimensionMismatch("size of $FSA is not fitting array $(typeof(a)), with size: $(size(a))"))
        $expr
    end
end

call{FSA <: FixedVectorNoTuple}(::Type{FSA}, a::Tuple, b::Tuple...) = throw(DimensionMismatch("$FSA can't be constructed from $a"))
call{FSA <: FixedVectorNoTuple}(::Type{FSA}, a::Tuple) = FSA(a...)
#call{FSA <: FixedArray, T}(::Type{FSA}, a::T...) = FSA(a) # all a's have the same type, so we can just use that tuple

@generated function call{FSA <: FixedArray, X}(::Type{FSA}, a::X)
    SZ      = size_or(FSA, (1,))
    ElType  = eltype_or(FSA, a)
    Len     = prod(SZ)
    if FSA <: FixedVectorNoTuple
        return :($FSA($(ntuple(i-> :($ElType(a)), Len)...)))
    else
        name = similar(FSA, ElType, SZ)
        expr = fill_tuples_expr((inds...)->:($ElType(a[1])), SZ)
        return :($name($expr))
    end
end

@generated function call{FSA <: FixedArray}(::Type{FSA}, a...)
    SZ     = size_or(FSA, (length(a),))
    ElType = eltype_or(FSA, promote_type(a...))
    if FSA <: FixedVectorNoTuple
        length(a) != prod(SZ) && throw(DimensionMismatch("can't construct $FSA with $(length(a)) arguments. Args: $a"))
        return :($FSA(map($ElType, a)...))
    else
        FSA <: Mat && return :($FSA(a))
        inner = any(x-> x != ElType, a) ? :(map($ElType, a)) : :(a)
        name = similar(FSA, ElType, SZ)
        return :($name($inner))
    end
end



@inline zero{FSA <: FixedArray}(::Type{FSA}) = map(ConstFunctor(zero(eltype(FSA))), FSA)
zero(fsa::FixedArray) = zero(typeof(fsa))
@inline one{FSA <: FixedArray}(::Type{FSA})  = map(ConstFunctor(one(eltype(FSA))), FSA)
@inline eye{FSA <: FixedArray}(::Type{FSA})  = map(EyeFunc{eltype(FSA)}, FSA)
@inline unit{FSA <: FixedVector}(::Type{FSA}, i::Integer) = map(UnitFunctor(i, eltype(FSA)), FSA)

@inline rand{FSA <: FixedArray}(x::Type{FSA}) = map(RandFunctor{eltype(FSA)}, FSA)
@inline rand{FSA <: FixedArray}(x::Type{FSA}, range::Range) = (T = eltype(FSA) ; map(RandFunctor{T(first(range)):T(step(range)):T(last(range))}, FSA)) # there's no easy way to convert eltypes of ranges (I think)

"""
Macro `fsa` helps to create fixed size arrays like Julia arrays.
E.g.
```
@fsa([1 2 3;4 5 6])
@fsa([a,2,3]) # you can also use variables
@fsa([a 2 3])
```
"""
macro fsa(expr)
    if expr.head == :vect
        result = Expr(:call, :Vec, expr.args...)
    elseif expr.head == :hcat
        result = Expr(:call, :Mat, [Expr(:tuple, a) for a in expr.args]...)
    elseif expr.head == :vcat
        if isa(expr.args[1], Expr) && expr.args[1].head == :row
            cols = [Any[] for _=1:length(expr.args[1].args)]
            for row in expr.args
                for (j, a) in enumerate(row.args)
                    push!(cols[j], a)
                end
            end
            result = Expr(:call, :Mat, Expr(:tuple, [Expr(:tuple, col...) for col in cols]...))
        else
            result = Expr(:call, :Vec, expr.args...)
        end
    end
    esc(result)
end
