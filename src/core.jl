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

similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T}, n::Tuple) = similar(FSA, T, n...)
@generated function similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T}, n::Int)
    name = basetype(FSA)
    :($name{n, T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T})
    name = basetype(FSA)
    :($name{$(FSA.parameters[1]), T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVectorNoTuple, T}(::Type{FSA}, ::Type{T})
    name = basetype(FSA)
    :($name{T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVectorNoTuple, T}(::Type{FSA}, ::Type{T}, n::Int)
    name = basetype(FSA)
    :($name{T, $(FSA.parameters[3:end]...)})
end

@generated function get_tuple{N, T}(f::FixedVectorNoTuple{N, T})
    :(tuple($(ntuple(i->:(f[$i]), N)...)))
end
function get_tuple(f::FixedArray)
    f._
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
    simtype = similar(FSA, etype, shape)
    converted_elements = convert_nested_tuple_expr(etype, :elements, elements)
    :($simtype($converted_elements))
end

