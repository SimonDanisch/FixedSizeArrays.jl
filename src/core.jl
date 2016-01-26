abstract FixedArray{T, NDim, SIZE}
abstract MutableFixedArray{T, NDim, SIZE} <: FixedArray{T, NDim, SIZE}

typealias MutableFixedVector{T, CARDINALITY} MutableFixedArray{T, 1, Tuple{CARDINALITY}}
typealias MutableFixedMatrix{T, M, N}        MutableFixedArray{T, 2, Tuple{M,N}}

typealias FixedVector{CARDINALITY, T}        FixedArray{T, 1, Tuple{CARDINALITY,}}
typealias FixedMatrix{Row, Column, T}        FixedArray{T, 2, Tuple{Row, Column}}

abstract FixedVectorNoTuple{CARDINALITY, T} <: FixedVector{CARDINALITY, T}
export FixedVectorNoTuple


_size{T <: Tuple}(::Type{T})						= (T.parameters...)
_size{N, N2}(::Type{Tuple{N, N2}})					= (N,N2)
_size{N}(::Type{Tuple{N}})							= (N,)

eltype{T <: FixedArray}(A::Type{T})                 = eltype_or(T, Any)
eltype{T <: FixedArray,N,SZ}(A::FixedArray{T,N,SZ}) = T


function length{T <: FixedArray}(A::Type{T})
    prod(size(T))
end
length{T <: FixedArray}(A::T) = length(T)

endof{T <: FixedArray}(A::Type{T}) = length(T)
endof{T <: FixedArray}(A::T) = endof(T)


@generated function ndims{T <: FixedArray}(A::Type{T})
    :($(fsa_abstract(T).parameters[2]))
end
ndims{T <: FixedArray}(A::T) = ndims(T)


size{T,N,SZ}(A::Type{FixedArray{T,N,SZ}})           = _size(SZ)
size{T <: FixedArray}(A::Type{T})                   = size_or(T, ())
size{T <: FixedArray}(A::T)                         = size(T)

size{T <: FixedArray}(A::Type{T}, d::Integer)       = size(T)[d]
size{T <: FixedArray}(A::T, d::Integer)             = size(T, d)

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


similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T}, n::Tuple) = similar(FSA, T, n...)
@generated function similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T}, n::Int)
    name = parse(string("Main.", FSA.name))
    :($name{n, T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVector, T}(::Type{FSA}, ::Type{T})
    name = parse(string("Main.", FSA.name))
    :($name{$(FSA.parameters[1]), T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVectorNoTuple, T}(::Type{FSA}, ::Type{T})
    name = parse(string("Main.", FSA.name))
    :($name{T, $(FSA.parameters[3:end]...)})
end
@generated function similar{FSA <: FixedVectorNoTuple, T}(::Type{FSA}, ::Type{T}, n::Int)
    name = parse(string("Main.", FSA.name))
    :($name{T, $(FSA.parameters[3:end]...)})
end

@generated function get_tuple{N, T}(f::FixedVectorNoTuple{N, T})
    :(tuple($(ntuple(i->:(f[$i]), N)...)))
end
function get_tuple(f::FixedArray)
    f.(1) # a little wonky, but there really no much sense if isn't the only field
end
