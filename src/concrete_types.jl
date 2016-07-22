#-------------------------------------------------------------------------------
# Canonical fixed size array / linear algebra types
"""
    Vec{N,T}

Generic length-`N` vector with element type `T`
"""
immutable Vec{N,T} <: FixedVector{N,T}
    values::NTuple{N,T}
end

"""
    Mat{N,M,T}

Generic `N×M` matrix with element type `T`
"""
immutable Mat{N,M,T} <: FixedMatrix{N,M,T}
    values::NTuple{M,NTuple{N,T}}
end

"""
    FArray3{N,M,P,T}

Generic rank three `N×M×P` tensor with element type `T`
"""
immutable FArray3{N,M,P,T} <: FixedArray3{N,M,P,T}
    values::NTuple{P,NTuple{M,NTuple{N,T}}}
end

"""
    FArray3{N,M,P,Q,T}

Generic rank four `N×M×P×Q` tensor with element type `T`
"""
immutable FArray4{N,M,P,Q,T} <: FixedArray4{N,M,P,Q,T}
    values::NTuple{Q,NTuple{P,NTuple{M,NTuple{N,T}}}}
end


# Return one of the standard array / linear algebra types above.
@pure function default_similar_type{T}(::Type{T}, sz::Tuple)
    if length(sz) == 1
        return Vec{sz[1],T}
    elseif length(sz) == 2
        return Mat{sz[1],sz[2],T}
    elseif length(sz) == 3
        return FArray3{sz[1],sz[2],sz[3],T}
    elseif length(sz) == 4
        return FArray4{sz[1],sz[2],sz[3],sz[4],T}
    else
        throw(ArgumentError("No built in FixedArray type is implemented for eltype $T and size $sz"))
    end
end

#-------------------------------------------------------------------------------
# Other types
"""
    Point{N,T}

Generic member of an `N`-dimensional Cartesian space with element type `T`
"""
immutable Point{N,T} <: FixedVector{N,T}
    values::NTuple{N,T}
end
