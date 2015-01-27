immutable IMax <: Base.Func{2} end
call(::IMax, x, y) = max(x,y)
immutable IMin <: Base.Func{2} end
call(::IMin, x, y) = min(x,y)
Base.maximum{T}(v::Vector{Vector3{T}}) = reduce(IMax(), Vector3(typemin(T)), v)
Base.minimum{T}(v::Vector{Vector3{T}}) = reduce(IMin(), Vector3(typemax(T)), v)