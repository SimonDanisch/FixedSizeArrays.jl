
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


#=
Constructors for homogenous non tuple arguments
Is unrolled for the first 4, since a::T... leads to slow performance
=#
call{FSA <: FixedArray, T}(::Type{FSA}, a::T) = FSA(NTuple{1,T}((a,)))
call{FSA <: FixedArray, T}(::Type{FSA}, a::T, b::T) = FSA(NTuple{2,T}((a,b)))
call{FSA <: FixedArray, T}(::Type{FSA}, a::T, b::T, c::T) = FSA(NTuple{3,T}((a,b,c)))
call{FSA <: FixedArray, T}(::Type{FSA}, a::T, b::T, c::T, d::T) = FSA(NTuple{4,T}((a,b,c,d)))
call{FSA <: FixedArray, T}(::Type{FSA}, a::T...) = FSA(a)

immutable ParseFunctor{T, S <: AbstractString} <: Func{1}
    t::Type{T}
    a::Vector{S}
end
call{T}(pf::ParseFunctor{T}, i::Int) = parse(T, pf.a[i])
call(pf::ParseFunctor{Void}, i::Int) = parse(pf.a[i])


"""
Constructs a fixedsize array from a Base.Array
"""
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
    quote
        $SZ != size(a) && throw(DimensionMismatch("size of $FSA is not fitting array $(typeof(a)), with size: $(size(a))"))
        $expr
    end
end


"""
Constructor for singular arguments.
Can be a tuple, is not declared as that, because it will generate ambigouities
and overwrites the default constructor.
"""
@generated function call{FSA <: FixedArray, X}(::Type{FSA}, a::X)
    ND = ndims(FSA)
    if X <: Tuple
        types_svec = a.parameters
        if all(x-> x <: Tuple, types_svec) && ND == 1
            return :(throw(
                DimensionMismatch("tried to construct $FSA from $a. I can't allow that!")
            ))
        end
        orlen = length(types_svec)
        ortyp = promote_type(types_svec...)
    else
        orlen = 1
        ortyp = a
    end
    SZ      = size_or(FSA, (orlen, ntuple(x->1, ND-1)...))
    T       = eltype_or(FSA, ortyp)
    FSAT    = similar(FSA, T, SZ)
    if X <: Tuple
        expr = fill_tuples_expr((inds...)->:($T(a[$(inds[1])])), SZ)
    else
        expr = fill_tuples_expr((inds...)->:($T(a)), SZ)
    end
    return :($FSAT($expr))
end
"""
Constructors for heterogenous multiple arguments.
E.g. 1, 2f0, 4.0
"""
@generated function call{FSA <: FixedArray}(::Type{FSA}, a...)
    SZ     = size_or(FSA, (length(a),))
    ElType = eltype_or(FSA, promote_type(a...))
    all(x-> x <: Tuple, a) && return :( $FSA(a) ) # TODO be smarter about this
    any(x-> x != ElType, a) && return :($FSA(map($ElType, a)))
    return :($FSA(a))
end

"""
Construction from other FixedSizeArrays + X arguments
E.g. Vec4f0(Vec3f0(1), 0)
"""
@generated function call{FSA <: FixedArray, T1 <: FixedArray}(::Type{FSA}, a::T1, b...)
    if isempty(b) # this is the conversion constructor for 2 FSA's
        #easiest way is to just construct from the tuple
        expr = :(FSA(get_tuple(a)))
        if size_or(FSA, nothing) == nothing # no concrete size
            return expr
        else #has a size
            len1 = size(FSA, 1)
            len2 = size(T1, 1)
            if len1 < len2 # we need to shrink
                return :(FSA(get_tuple(a)[1:$len1]))
            elseif len1==len2
                return expr
            else
                return :(throw(DimensionMismatch(
                    "tried to create $FSA from $T1. The latter has too many elements"
                )))
            end
        end
    end
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
    return :( $FSA(a..., b...) )
end

@inline zero{FSA <: FixedArray}(::Type{FSA}) = map(ConstFunctor(zero(eltype(FSA))), FSA)
zero(fsa::FixedArray) = zero(typeof(fsa))
@inline one{FSA <: FixedArray}(::Type{FSA})  = map(ConstFunctor(one(eltype(FSA))), FSA)
@inline eye{FSA <: FixedArray}(::Type{FSA})  = map(EyeFunc{eltype(FSA)}, FSA)
@inline unit{FSA <: FixedVector}(::Type{FSA}, i::Integer) = map(UnitFunctor(i, eltype(FSA)), FSA)

@inline rand{FSA <: FixedArray}(m::MersenneTwister, x::Type{FSA}) = map(MersenneFunctor{eltype(FSA)}(m), FSA)
@inline function rand{FSA <: FixedArray}(x::Type{FSA}, range::Range)
    T = eltype(FSA)
    map(RandFunctor(T(first(range)):T(step(range)):T(last(range))), FSA) # there's no easy way to convert eltypes of ranges (I think)
end
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
