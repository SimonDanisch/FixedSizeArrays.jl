import Base: getindex, convert, eltype
abstract FixedSizeVector{T, SZ, NDIM}
abstract FixedSizeWrapper{T}

getindex(a::FixedSizeVector, i::Int)      = getfield(a, i)
getindex(a::FixedSizeWrapper, i::Int)     = getfield(a.data, i)

eltype{T <: FixedSizeVector}(A::Type{T})  = A.types[1]
eltype{T <: FixedSizeWrapper}(a::Type{T}) = eltype(T.types[1])


macro type_accessors(typ, fields)
    result = Any[]
    for elem in fields.args
        push!(result, quote 
            fieldindex{T <: $typ}(::Type{T}, ::Type{$(elem.args[1])}) = $(elem.args[2])
        end)
    end
    esc(Expr(:block, result...))
end




abstract Dimension{T} <: Number
Base.convert{T <: Real}(::Type{T}, a::Dimension) = convert(T, a.val)

immutable Red{T} <: Dimension{T}
    val::T
end
immutable Green{T} <: Dimension{T}
    val::T
end
immutable Blue{T} <: Dimension{T}
    val::T
end


immutable X{T} <: Dimension{T}
    val::T
end
immutable Y{T} <: Dimension{T}
    val::T
end
immutable Z{T} <: Dimension{T}
    val::T
end
immutable W{T} <: Dimension{T}
    val::T
end
stagedfunction getindex{T <: Dimension}(a::FixedSizeWrapper, key::Type{T})
    index = fieldindex(a, T)
    :(T(a.data[$index])) 
end
immutable D1{T} <: FixedSizeVector{T, (1,), 1}
    a::T
end
immutable D2{T} <: FixedSizeVector{T, (2,), 1}
    a::T
    b::T
end
immutable D3{T} <: FixedSizeVector{T, (3,), 1}
    a::T
    b::T
    c::T
end
immutable D4{T} <: FixedSizeVector{T, (4,), 1}
    a::T
    b::T
    c::T
    d::T
end

immutable RGB{T} <: FixedSizeWrapper{D3{T}}
    data::D3{T}
end
@type_accessors RGB (Red => 1, Green => 2, Blue => 3)

RGB{T}(a::T,b::T,c::T) = RGB(D3{T}(a,b,c))
RGB{T}(i::T)           = RGB(D3{T}(i,i,i))


abstract GeometricVector{T} <: FixedSizeWrapper{T}
@type_accessors GeometricVector (X => 1, Y => 2, Z => 3, W => 4)

immutable Point{T <: FixedSizeVector} <: GeometricVector{T}
    data::T
end

immutable Normal{T <: FixedSizeVector} <: GeometricVector{T}
    data::T
end

immutable Vertex{T <: FixedSizeVector} <: GeometricVector{T}
    data::T
end



stagedfunction getindex{V <: FixedSizeWrapper, D <: Dimension}(x::Array{V}, key::Type{D})
    
    typ   = symbol(string(D.name))
    eltyp = eltype(V)
    typ   = eval(:($typ{$eltyp}))

    :(reinterpret($typ, reshape([elem[D] for elem in x], size(x))))
end
#=
function Base.setindex!{T <: AbstractFixedVector, ElType}(a::Vector{T}, x::ElType, i::Integer, accessor::Integer)
  @assert eltype(T) == ElType # ugly workaround for not having triangular dispatch
  @assert length(a) >= i
  cardinality = length(T)
  @assert accessor <= cardinality
  ptr = convert(Ptr{ElType}, pointer(a))
  unsafe_store!(ptr, x, ((i-1)*cardinality)+accessor)
end
function Base.setindex!{T <: AbstractFixedVector, ElType}(a::Vector{T}, x::Vector{ElType}, i::Integer, accessor::UnitRange)
  @assert eltype(T) == ElType
  @assert length(a) >= i
  cardinality = length(T)
  @assert length(accessor) <= cardinality
  ptr = convert(Ptr{ElType}, pointer(a))
  unsafe_copy!(ptr + (sizeof(ElType)*((i-1)*cardinality)), pointer(x), length(accessor))
end
function setindex1D!{T <: AbstractFixedVector, ElType}(a::Union(Matrix{T}, Vector{T}), x::ElType, i::Integer, accessor::Integer)
    @assert length(a) >= i "Out of Bounds. 1D index: $i, Matrix: , $(typeof(a)), $length: $(length(a)) size: $(size(a))"

    cardinality = length(T)
    if length(accessor) > cardinality
        error("Out of Bounds. 1D index: ", i, " Matrix: ", typeof(a), " length: ", length(a), " size: ", size(a))
    end

  ptr = convert(Ptr{eltype(T)}, pointer(a))
  unsafe_store!(ptr, convert(eltype(T), x), ((i-1)*cardinality)+accessor)
end
function setindex1D!{T <: AbstractFixedVector, ElType}(a::Union(Matrix{T}, Vector{T}), x::Vector{ElType}, i::Integer, accessor::UnitRange)
   if length(a) < i
     error("Out of Bounds. 1D index: ", i, " Matrix: ", typeof(a), " length: ", length(a), " size: ", size(a))
   end
   cardinality = length(T)
   if length(accessor) > cardinality
     error("Out of Bounds. 1D index: ", i, " Matrix: ", typeof(a), " length: ", length(a), " size: ", size(a))
   end
   eltp  = eltype(T)
   x     = convert(Vector{eltp}, x)
   ptr   = convert(Ptr{eltp}, pointer(a))
   unsafe_copy!(ptr + (sizeof(eltp)*((i-1)*(cardinality)+first(accessor-1))), pointer(x), length(accessor))
end
=#
const v = D3(1f0, 2f0, 3f0)
@show const p = Point(v)
@show p[X]
const a = RGB(v)
@show a[1]

const ttt = [RGB(1f0) for i=1:100, j=1:100]
println(typeof(ttt[Red]))