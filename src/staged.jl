
reduce{F <: Func{2}, T}(f::Type{F}, a::AbstractFixedVector{T, 1}) = a

stagedfunction reduce(f::Func{2}, a::AbstractFixedVector)
    alength = length(a)
    @assert alength >= 2 "type $a not long enough for reduce (needs at least two elements)"
    quote 
        s = f(a[1], a[2])
        $([:(s = f(s, a[$i])) for i=3:alength]...)
    end
end

# Needed to 
accessor{T <: AbstractFixedArray}(arg::Type{T}, i::Integer, name::Symbol) = :($name[$i])
accessor{T                      }(arg::Type{T}, i::Integer, name::Symbol) = :($name)

stagedfunction map{T <: AbstractFixedArray}(f::Func{1}, a::T)
    expr = ntuple(i -> :( f( $(accessor(a, i, :a))) ), length(T))
    type_name = a.name.name
    :($type_name($(expr...)))
end
stagedfunction map{T <: AbstractFixedArray}(f::Func{2}, a::T, b::T)
    expr = ntuple(i -> :( f( $(accessor(a, i, :a)), $(accessor(b, i, :b))) ), length(T))
    :($a($(expr...)))
end
stagedfunction map{T <: AbstractFixedArray, T2 <: Real}(f::Func{2}, a::T, b::T2)
    expr = ntuple(i -> :( f( $(accessor(a, i, :a)), $(accessor(b, i, :b))) ), length(T))
    :($a($(expr...)))
end
stagedfunction map{T <: AbstractFixedArray, T2 <: Real}(f::Func{2}, a::T2, b::T)
    expr = ntuple(i -> :( f( $(accessor(a, i, :a)), $(accessor(b, i, :b))) ), length(T))
    :($a($(expr...)))
end


stagedfunction convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::T2)
    @assert sizeof(a) == sizeof(b) "Type $a ($(sizeof(a))) doesn't have the same bit size as type $b ($(sizeof(b)))"
    :(reinterpret(a, [b])[1])
end
stagedfunction convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(a) == sizeof(b) "Type $a ($(sizeof(a))) doesn't have the same bit size as type $b ($(sizeof(b)))"
    :(reinterpret(a, b))
end