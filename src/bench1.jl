abstract AbstractFixedArray{T,N,SZ}
typealias AbstractFixedVector{T, C} AbstractFixedArray{T, 1, (C,)}
typealias AbstractFixedMatrix{T, M, N} AbstractFixedArray{T, 2, (M, N)}
importall Base

eltype{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})           = T
length{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})           = prod(SZ)
ndims{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})            = N
size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})             = SZ
size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ}, d::Integer) = SZ[d]

# Ugly workaround for not having triangular dispatch:
eltype{T <: AbstractFixedArray}(A::Type{T})                 = A.types[1]
length{T <: AbstractFixedVector}(A::Type{T})                = length(A.types)
length{T <: AbstractFixedMatrix}(A::Type{T})                = prod(size(A))

ndims{T <: AbstractFixedVector}(A::Type{T})                 = 1
ndims{T <: AbstractFixedMatrix}(A::Type{T})                 = 2

size{T <: AbstractFixedVector}(A::Type{T})                  = (length(A),)
size{T <: AbstractFixedVector}(A::Type{T}, d::Integer)      = (length(A),) # should throw an error!?
size{T <: AbstractFixedMatrix}(A::Type{T})                  = (length(A.types), length(A.types[1]))
size{T <: AbstractFixedMatrix}(A::Type{T}, d::Integer)      = size(A)[d] 

getindex{T,C}(A::AbstractFixedVector{T, C}, i::Integer)                 = getfield(A, i)
getindex{T,M,N}(A::AbstractFixedMatrix{T, M,N}, i::Integer, j::Integer) = getfield(getfield(A, i),j)

columntype{T,N,SZ}(x::AbstractFixedMatrix{T,N,SZ}) = typeof(x[1])
columntype{T <: AbstractFixedMatrix}(x::Type{T})   = first(x.types)

immutable Vector4{T}<:AbstractFixedVector{T,4}
    e1::T
    e2::T
    e3::T
    e4::T
end
immutable Matrix4x4{T}<:AbstractFixedMatrix{T,4,4}
    c1::Vector4{T}
    c2::Vector4{T}
    c3::Vector4{T}
    c4::Vector4{T}
end
abstract Func{N}
immutable AddFun <: Func{2} end
call(::AddFun, x, y) = +(x,y)

const l = Vector4(1,2,3,4)
const m = Matrix4x4(l,l,l,l)

reduce{T}(f::Func{2}, a::AbstractFixedVector{T, 1}) = a
stagedfunction reduce(f::Func{2}, a::AbstractFixedVector)
    accessors = names(a)
    expr = Any[]
    gfa = parse("getfield(a, :$(accessors[1]))")
    gfb = parse("getfield(a, :$(accessors[2]))")
    push!(expr, :(s = call(f, $gfa, $gfb)))
    for i=3:length(accessors)
        gfa = parse("getfield(a, :$(accessors[i]))") # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(s = call(f, s, $gfa)))
    end
    Expr(:block, expr...)
end

println(reduce(AddFun(), l))

function test(l, n)
    for i=1:n
        reduce(AddFun(), l)
    end
end

test(l, 1)
@time test(l, 10^6)
@time test(l, 10^6)
@time test(l, 10^6)
@time test(l, 10^6)