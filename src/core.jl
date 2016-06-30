abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, Tuple{CARDINALITY}}
typealias MutableFixedMatrix{T, M, N}        MutableFixedArray{T, 2, Tuple{M,N}}

typealias FixedVector{CARDINALITY, T}        FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{Row, Column, T}        FixedArray{T, 2, Tuple{Row, Column}}

abstract FixedVectorNoTuple{CARDINALITY, T} <: FixedVector{CARDINALITY, T}
export FixedVectorNoTuple


_size{T <: Tuple}(::Type{T}) = (T.parameters...)
_size{N, N2}(::Type{Tuple{N, N2}}) = (N,N2)
_size{N}(::Type{Tuple{N}}) = (N,)

eltype{T <: FixedArray}(A::Type{T}) = eltype_or(T, Any)
eltype{T <: FixedArray,N,SZ}(A::FixedArray{T,N,SZ}) = T


length{T <: FixedArray}(A::Type{T}) = prod(size(T))
length{T <: FixedArray}(A::T) = length(T)

endof{T <: FixedArray}(A::Type{T}) = length(T)
endof{T <: FixedArray}(A::T) = endof(T)


@generated function ndims{T <: FixedArray}(A::Type{T})
    :($(fsa_abstract(T).parameters[2]))
end
ndims{T <: FixedArray}(A::T) = ndims(T)


size{T,N,SZ}(A::Type{FixedArray{T,N,SZ}}) = _size(SZ)
size{T <: FixedArray}(A::Type{T}) = size_or(T, ())
size{T <: FixedArray}(A::T) = size(T)

size{T <: FixedArray}(A::Type{T}, d::Integer) = size(T)[d]
size{T <: FixedArray}(A::T, d::Integer) = size(T, d)

@generated function fsa_abstract{FSA <: FixedArray}(::Type{FSA})
    ff = FSA
    while ff.name.name != :FixedArray
       ff = supertype(ff)
    end
    :($ff)
end
@generated function size_or{FSA <: FixedArray}(::Type{FSA}, OR)
    fsatype = fsa_abstract(FSA)
    sz = fsatype.parameters[3]
    isa(sz, TypeVar) && return :(OR)
    any(x->isa(x, TypeVar), sz.parameters) && return :(OR)
    :($(_size(sz)))
end
@generated function eltype_or{FSA <: FixedArray}(::Type{FSA}, OR)
    fsatype = fsa_abstract(FSA)
    T = fsatype.parameters[1]
    isa(T, TypeVar) && return :(OR)
    :($T)
end
@generated function ndims_or{FSA <: FixedArray}(::Type{FSA}, OR)
    fsatype = fsa_abstract(FSA)
    N = fsatype.parameters[2]
    isa(N, TypeVar) && return :(OR)
    :($N)
end

# Iterator
start(A::FixedArray) = 1
function next(A::FixedArray, state::Integer)
    @inbounds x = A[state]
    (x, state+1)
end
done(A::FixedArray, state::Integer) = length(A) < state

@generated function basetype{T<:FixedArray}(::Type{T})
    :($(T.name.primary))
end

"""
    similar_type(::Type{FSA}, [::Type{T}=eltype(FSA)], [sz=size(FSA)])

Given an array type `FSA`, element type `T` and size `sz`, return a `FixedArray`
subtype which is as similar as possible.  `similar_type` is used in the same
spirit as `Base.similar` to store the results of `map()` operations, etc.
(`similar` cannot work here, because the types are generally immutable.)

By default, `similar_type` introspects `FSA` to determine whether `T` and `sz`
can be used; if not a canonical FixedArray container is returned instead.
"""
@pure function similar_type{T}(::Type{FixedArray}, ::Type{T}, sz::Tuple)
    if length(sz) == 1
        return Vec{sz[1],T}
    elseif length(sz) == 2
        return Mat{sz[1],sz[2],T}
    else
        throw(ArgumentError("No built in FixedArray type is implemented for eltype $T and size $sz"))
    end
end

@pure function similar_type{FSA <: FixedArray, T}(::Type{FSA}, ::Type{T}, sz::Tuple)
    fsa_size = fsa_abstract(FSA).parameters[3].parameters
    if eltype(FSA) == T && fsa_size == sz
        return FSA # Common case optimization
    end

    # The default implementation for similar_type is follows.  It involves a
    # fair bit of crazy type introspection: We check whether the type `FSA` has
    # the necessary type parameters to replace with `T` and `sz`, and if so
    # figure out how to do the replacement.  It's complicated because users may
    # arbitrarily rearrange type parameters in their subtypes, and possibly
    # even add new type parameters which aren't related to the abstract
    # FixedArray but should be preserved.

    # Propagate the available type parameters of FSA down to the abstract base
    # FixedArray as `TypeVar`s.
    pritype = FSA.name.primary
    fsatype = fsa_abstract(pritype)
    T_parameter    = fsatype.parameters[1]
    ndim_parameter = fsatype.parameters[2]
    sz_parameter   = fsatype.parameters[3]
    sz_parameters  = fsatype.parameters[3].parameters

    # Figure out whether FSA can accommodate the new eltype `T` and size `sz`.
    # If not, delegate to the fallback by default.
    if !((eltype(FSA) == T          || isa(T_parameter,    TypeVar)) &&
         (ndims(FSA)  == length(sz) || isa(ndim_parameter, TypeVar)) &&
         (fsa_size    == sz         || all(i -> (sz[i] == fsa_size[i] || isa(sz_parameters[i],TypeVar)), 1:length(sz))))
        return similar_type(FixedArray, T, sz)
    end

    # Iterate type parameters, replacing as necessary with T and sz
    params = collect(FSA.parameters)
    priparams = pritype.parameters
    for i=1:length(params)
        if priparams[i] === T_parameter
            params[i] = T
        elseif priparams[i] === ndim_parameter
            params[i] = length(sz)
        elseif priparams[i] === sz_parameter
            params[i] = Tuple{sz...}
        else
            for j = 1:length(sz_parameters)
                if priparams[i] === sz_parameters[j]
                    params[i] = sz[j]
                end
            end
        end
    end
    pritype{params...}
end

# similar_type versions with defaulted eltype and size
@pure similar_type{FSA <: FixedArray, T}(::Type{FSA}, ::Type{T}) = similar_type(FSA, T, size(FSA))
@pure similar_type{FSA <: FixedArray}(::Type{FSA}, sz::Tuple) = similar_type(FSA, eltype(FSA), sz)

# Deprecated similar() -> similar_type()
function get_tuple(f::FixedArray)
    Base.depwarn("get_tuple(f::FixedArray) is deprecated, use Tuple(f) instead", :get_tuple)
    Tuple(f)
end
function similar{FSA <: FixedArray}(::Type{FSA}, args...)
    Base.depwarn("similar{FSA<:FixedArray}(::Type{FSA}, ...) is deprecated, use similar_type instead", :similar)
    similar_type(FSA, args...)
end
similar{FSA <: FixedArray}(::Type{FSA}, sz::Int...) = similar(FSA, eltype(FSA), sz)
similar{FSA <: FixedArray,T}(::Type{FSA}, ::Type{T}, sz::Int...) = similar(FSA, T, sz)

@compat @generated function (::Type{T}){T<:Tuple, N, T1}(f::FixedVectorNoTuple{N, T1})
    return Expr(:tuple, ntuple(i->:(f[$i]), N)...)
end
@compat function (::Type{T}){T<:Tuple}(f::FixedArray)
    getfield(f, 1)
end


# Infrastructure for construct_similar
"Compute promoted element type of the potentially nested Tuple type `ty`"
function promote_type_nested(ty)
    if ty.name.primary === Tuple
        promote_type([promote_type_nested(t) for t in ty.parameters]...)
    else
        ty
    end
end

"""
Construct tuple expression converting inner elements of the nested Tuple type
`ty` with name `varexpr` to the scalar `outtype`
"""
function convert_nested_tuple_expr(outtype, varexpr, ty)
    if ty.name.primary === Tuple
        Expr(:tuple, [convert_nested_tuple_expr(outtype, :($varexpr[$i]), t)
                      for (i,t) in enumerate(ty.parameters)]...)
    else
        :(convert($outtype, $varexpr))
    end
end

"""
Compute the N-dimensional array shape of a nested Tuple if it were used as
column-major storage for a FixedArray.
"""
function nested_Tuple_shape(ty)
    if ty.name.primary !== Tuple
        return 0
    end
    subshapes = [nested_Tuple_shape(t) for t in ty.parameters]
    if isempty(subshapes)
        return ()
    end
    if any(subshapes .!= subshapes[1])
        throw(DimensionMismatch("Nested tuples must have equal length to form a FixedSizeArray"))
    end
    if subshapes[1] == 0
        return (length(subshapes),) # Scalar elements
    end
    return (subshapes[1]..., length(subshapes))
end

"""
    construct_similar(::Type{FSA}, elements::Tuple)

Construct FixedArray as similar as possible to `FSA`, but with shape given by
the shape of the nested tuple `elements` and with `eltype` equal to the
promoted element type of the nested tuple `elements`.
"""
@generated function construct_similar{FSA <: FixedArray}(::Type{FSA}, elements::Tuple)
    etype = promote_type_nested(elements)
    shape = nested_Tuple_shape(elements)
    outtype = similar_type(FSA, etype, shape)
    converted_elements = convert_nested_tuple_expr(etype, :elements, elements)
    constructor_expr(outtype, converted_elements)
end
