#-------------------------------------------------------------------------------
# Canonical fixed size array / linear algebra types
"Canonical length-`N` vector with element type `T`"
immutable Vec{N,T} <: FixedVector{N,T}
    values::NTuple{N,T}
end

"Canonical N×M matrix with element type `T`"
immutable Mat{N,M,T} <: FixedMatrix{N,M,T}
    values::NTuple{M,NTuple{N,T}}
end

"Canonical rank three tensor with element type `T`"
immutable Ar3{N,M,P,T} <: FixedArray3{N,M,P,T}
    values::NTuple{P,NTuple{M,NTuple{N,T}}}
end

"Canonical rank four tensor with element type `T`"
immutable Ar4{N,M,P,Q,T} <: FixedArray4{N,M,P,Q,T}
    values::NTuple{Q,NTuple{P,NTuple{M,NTuple{N,T}}}}
end


# Return one of the standard array / linear algebra types above.
@pure function default_similar_type{T}(::Type{T}, sz::Tuple)
    if length(sz) == 1
        return Vec{sz[1],T}
    elseif length(sz) == 2
        return Mat{sz[1],sz[2],T}
    elseif length(sz) == 3
        return Ar3{sz[1],sz[2],sz[3],T}
    elseif length(sz) == 4
        return Ar4{sz[1],sz[2],sz[3],sz[4],T}
    else
        throw(ArgumentError("No built in FixedArray type is implemented for eltype $T and size $sz"))
    end
end

#-------------------------------------------------------------------------------
# Other types
"Generic member of an `N`-dimensional Cartesian space"
immutable Point{N,T} <: FixedVector{N,T}
    values::NTuple{N,T}
end
