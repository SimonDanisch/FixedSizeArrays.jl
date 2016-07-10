#-------------------------------------------------------------------------------
# Abstract type hierarchy

"""
    FixedArray{T,D}

Super type for all fixed size arrays of eltype `T` and dimension `D`.

The extents of the fixed dimensions aren't specified here due to technical
limitations.  (It'd have to be an `NTuple{Int,D}`, but `TypeVar`s representing
the extent of the dimensions in `FixedArray` subtypes can't be placed into this
tuple when subtyping.)
"""
abstract FixedArray{T,D}

# Define a subtype of FixedArray for each dimension.  Ideally we'd have the
# fixed size type parameters in FixedArray itself, but it's not clear how we
# can actually do this in a good way since type constructors are fairly
# restricted.
abstract FixedArray1{N,T}       <: FixedArray{T,1}
abstract FixedArray2{N,M,T}     <: FixedArray{T,2}
abstract FixedArray3{N,M,P,T}   <: FixedArray{T,3}
abstract FixedArray4{N,M,P,Q,T} <: FixedArray{T,4}
# Is there any point going to more dimensions than MAX_TUPLE_DEPTH?
abstract FixedArray5{N,M,P,Q,R,T} <: FixedArray{T,5}

typealias FixedVector FixedArray1
typealias FixedMatrix FixedArray2

abstract FixedVectorNoTuple{N,T} <: FixedVector{N,T}

#-------------------------------------------------------------------------------
# Helper functions

# Helper function: walk up the type tree until FixedArray is reached,
# propagating any TypeVars in the original type parameters.
@generated function fsa_abstract{FSA <: FixedArray}(::Type{FSA})
    ff = FSA
    while ff.name.name != :FixedArray
       ff = supertype(ff)
    end
    :($ff)
end
# Walk up the type tree until a subtype of FixedArray with known size is
# reached.  Propagates TypeVars correctly.
@generated function fsa_abstract_sized{FSA <: FixedArray}(::Type{FSA})
    if FSA.name.name == :FixedArray
        return :(error("FixedArray has no concretely sized supertype"))
    end
    T = FSA
    while supertype(T).name.name != :FixedArray
       T = supertype(T)
    end
    :($T)
end

# TODO: Deprecate the following?
@generated function size_or{FSA <: FixedArray}(::Type{FSA}, OR)
    FSA.name.name != :FixedArray || return :(OR)
    fsatype = fsa_abstract_sized(FSA)
    sz = (fsatype.parameters[1:end-1]...)
    any(x->isa(x, TypeVar), sz) && return :(OR)
    :($sz)
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


#-------------------------------------------------------------------------------
# General AbstractArray-like core interface
Base.eltype{T <: FixedArray}(::Type{T}) = eltype_or(T, Any)
Base.eltype{T,N}(::FixedArray{T,N}) = T

Base.length{T <: FixedArray}(::Type{T}) = prod(size(T))
Base.length{T <: FixedArray}(::T) = length(T)

Base.endof{T <: FixedArray}(::Type{T}) = length(T)
Base.endof{T <: FixedArray}(::T) = endof(T)

Base.ndims{T <: FixedArray}(::Type{T}) = ndims_or(T, nothing)
Base.ndims{T,N}(::FixedArray{T,N}) = N

Base.size{T <: FixedArray}(::Type{T}) = size_or(T, nothing)
Base.size{T <: FixedArray}(::T) = size(T)

Base.size{T <: FixedArray}(::Type{T}, d::Integer) = size(T)[d]
Base.size{T <: FixedArray}(::T, d::Integer) = size(T, d)


#-------------------------------------------------------------------------------
# Iterator
start(A::FixedArray) = 1
function next(A::FixedArray, state::Integer)
    @inbounds x = A[state]
    (x, state+1)
end
done(A::FixedArray, state::Integer) = length(A) < state

#-------------------------------------------------------------------------------
# Primary interface to unwrap elements: Tuple(FSA)
Base.Tuple(A::FixedArray) = getfield(A,1)

@generated function Base.Tuple{N,T}(A::FixedVectorNoTuple{N, T})
    return Expr(:tuple, ntuple(i->:(A[$i]), N)...)
end

#-------------------------------------------------------------------------------
# similar_type implementation
"""
    similar_type(::Type{FSA}, [::Type{T}=eltype(FSA)], [sz=size(FSA)])

Given an array type `FSA`, element type `T` and size `sz`, return a `FixedArray`
subtype which is as similar as possible.  `similar_type` is used in the same
spirit as `Base.similar` to store the results of `map()` operations, etc.
(`similar` cannot work here, because the types are generally immutable.)

By default, `similar_type` introspects `FSA` to determine whether `T` and `sz`
can be used; if not a canonical FixedArray container is returned instead.
"""
@pure function similar_type{FSA <: FixedArray, T}(::Type{FSA}, ::Type{T}, sz::Tuple)
    fsa_size = size_or(FSA, nothing)
    if eltype(FSA) == T && fsa_size == sz
        return FSA # Common case optimization
    end
    if fsa_size === nothing
        # Unsized
        return default_similar_type(T, sz)
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
    abstract_params = fsa_abstract_sized(pritype).parameters
    sz_parameters  = abstract_params[1:end-1]
    T_parameter    = abstract_params[end]

    # Figure out whether FSA can accommodate the new eltype `T` and size `sz`.
    # If not, delegate to the fallback by default.
    if !((eltype(FSA) == T          || isa(T_parameter, TypeVar)) &&
         (fsa_size    == sz         || (length(fsa_size) == length(sz) && all(i -> (sz[i] == fsa_size[i] || isa(sz_parameters[i],TypeVar)), 1:length(sz)))))
        return default_similar_type(T, sz)
    end

    # Iterate type parameters, replacing as necessary with T and sz
    newparams = collect(FSA.parameters)
    priparams = pritype.parameters
    for i=1:length(newparams)
        if priparams[i] === T_parameter
            newparams[i] = T
        else
            for j = 1:length(sz_parameters)
                if priparams[i] === sz_parameters[j]
                    newparams[i] = sz[j]
                end
            end
        end
    end
    pritype{newparams...}
end

# similar_type versions with defaulted eltype and size
@pure similar_type{FSA <: FixedArray, T}(::Type{FSA}, ::Type{T}) = similar_type(FSA, T, size(FSA))
@pure similar_type{FSA <: FixedArray}(::Type{FSA}, sz::Tuple) = similar_type(FSA, eltype(FSA), sz)


#-------------------------------------------------------------------------------
# construct_similar implementation

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
