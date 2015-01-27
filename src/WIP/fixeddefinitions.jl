module FixedSizeArrays
importall Base

export AbstractFixedArray
export AbstractFixedVector
export AbstractFixedMatrix


abstract AbstractFixedArray{T,N,SZ}
typealias AbstractFixedVector{T, C} AbstractFixedArray{T, 1, (C,)}
typealias AbstractFixedMatrix{T, Row, Column} AbstractFixedArray{T, 2, (Row, Column)}
importall Base

eltype{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})           = T
length{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})           = prod(SZ)
ndims{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})            = N
size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ})             = SZ
size{T,N,SZ}(A::AbstractFixedArray{T,N,SZ}, d::Integer) = SZ[d]

# Ugly workaround for not having triangular dispatch:
eltype{T <: AbstractFixedVector}(A::Type{T})                = A.types[1]
eltype{T <: AbstractFixedMatrix}(A::Type{T})                = eltype(A.types[1])

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
unit{T,C}(v::AbstractFixedVector{T,C}) = v/norm(v)
cross{T}(a::AbstractFixedVector{T, 2}, b::AbstractFixedVector{T, 2}) = a[1]*b[2]-a[2]*b[1]

cross{T}(a::AbstractFixedVector{T, 3}, b::AbstractFixedVector{T, 3}) = typeof(a)(a[2]*b[3]-a[3]*b[2], 
                                                     a[3]*b[1]-a[1]*b[3], 
                                                     a[1]*b[2]-a[2]*b[1])

stagedfunction zero{T <: AbstractFixedArray}(::Type{T}) 
    ttypes = T.types
    if !all(isleaftype, ttypes)
        error("please provide concrete type. Types given: ", accessors)
    end
    :($T($(map(x->:(zero($x)), ttypes)...)))
end

dot{T,C}(v1::AbstractFixedVector{T,C}, v2::AbstractFixedVector{T,C}) = sum(v1.*conj(v2))
norm{T,C}(v::AbstractFixedVector{T,C}) = sqrt(dot(v,v))
function norm{T,C}(v::AbstractFixedVector{T,C}, p::Number)
    if p == 1
        sum(abs(v))
    elseif p == 2
        norm(v)
    elseif p == Inf
        max(abs(v))
    else
        norm(copy(v), p)
    end
end

column{T, Cardinality}(v::AbstractFixedVector{T, Cardinality}) = v
column{T, Column, Row}(v::AbstractFixedMatrix{T, Row, Column}, i::Integer) = v[i]
row{T, Cardinality}(v::AbstractFixedVector{T, Cardinality}, i::Integer) = v[i]
stagedfunction row{T, Column, Row}(v::AbstractFixedMatrix{T, Row, Column}, i::Integer)
    ctype       = columntype(v)
    columnaccs  = names(v)
    expr = Any[]
    for colacc in columnaccs
        push!(expr, parse("getfield(v, :$(colacc))[i]"))
    end
    :($ctype($(expr...)))
end


abstract Func{N}

accessor_expr(variable, accessor_symbol) = parse("getfield($(variable), :$(accessor_symbol))")

reduce{T}(f::Func{2}, a::AbstractFixedVector{T, 1}) = a
stagedfunction reduce(f::Func{2}, a::AbstractFixedVector)
    accessors = names(a)
    expr = Any[]
    gfa = accessor_expr("a", accessors[1])
    gfb = accessor_expr("a", accessors[2])
    push!(expr, :(s = call(f, $gfa, $gfb)))
    for i=3:length(accessors)
        gfa = parse("getfield(a, :$(accessors[i]))") # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(s = call(f, s, $gfa)))
    end
    Expr(:block, expr...)
end
stagedfunction reduce(f::Func{2}, a::AbstractFixedMatrix)
    accessors = names(a)
    expr = Any[]
    gfa = accessor_expr("a", accessors[1])
    push!(expr, :(s = reduce(f, $gfa)))
    for i=3:length(accessors)
        gfa = parse("getfield(a, :$(accessors[i]))") # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(s = reduce(f, $gfa)))
    end
    Expr(:block, expr...)
end

stagedfunction map(f::Func{1}, a::AbstractFixedVector)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa)))
    end
    :($a($(expr...)))
end
stagedfunction map{T <: AbstractFixedVector}(f::Func{2}, a::T, b::T)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa, $gfb)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{2}, a::Real, b::AbstractFixedVector)
    accessors = names(b)
    expr = Any[]
    for elem in accessors
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, a, $gfb)))
    end
    :($b($(expr...)))
end
stagedfunction map(f::Func{2}, a::AbstractFixedVector, b::Real)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa, b)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{1}, a::AbstractFixedMatrix)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(map(f, $gfa)))
    end
    :($a($(expr...)))
end
stagedfunction map{T <: AbstractFixedMatrix}(f::Func{2}, a::T, b::T)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(map(f, $gfa, $gfb)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{2}, a::Real, b::AbstractFixedMatrix)
    accessors = names(b)
    expr = Any[]
    for elem in accessors
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(map(f, a, $gfb)))
    end
    :($a($(expr...)))
end

# operations
const unaryOps = (:-, :~, :conj, :abs, 
                  :sin, :cos, :tan, :sinh, :cosh, :tanh, 
                  :asin, :acos, :atan, :asinh, :acosh, :atanh,
                  :sec, :csc, :cot, :asec, :acsc, :acot,
                  :sech, :csch, :coth, :asech, :acsch, :acoth,
                  :sinc, :cosc, :cosd, :cotd, :cscd, :secd,
                  :sind, :tand, :acosd, :acotd, :acscd, :asecd,
                  :asind, :atand, :radians2degrees, :degrees2radians,
                  :log, :log2, :log10, :log1p, :exponent, :exp,
                  :exp2, :expm1, :cbrt, :sqrt, :square, :erf, 
                  :erfc, :erfcx, :erfi, :dawson, :ceil, :floor,
                  :trunc, :round, :significand, :lgamma, :hypot,
                  :gamma, :lfact, :frexp, :modf, :airy, :airyai,
                  :airyprime, :airyaiprime, :airybi, :airybiprime,
                  :besselj0, :besselj1, :bessely0, :bessely1,
                  :eta, :zeta, :digamma)

# vec-vec and vec-scalar
const binaryOps = (:.+, :.-,:.*, :./, :.\, :.^,:*,:/,
                   :.==, :.!=, :.<, :.<=, :.>, :.>=, :+, :-,
                   :min, :max,
                   :div, :fld, :rem, :mod, :mod1, :cmp,
                   :atan2, :besselj, :bessely, :hankelh1, :hankelh2, 
                   :besseli, :besselk, :beta, :lbeta)
# vec-vec only
const binaryOps2 = (:+,:-)
const reductions = ((:sum,:+),(:prod,:*),(:minimum,:min),(:maximum,:max))
for (callfun, reducefun) in reductions
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{2} end
        call(::$unicsymb, x, y) = $reducefun(x, y)
        $(callfun){T, N, SZ}(x::AbstractFixedArray{T, N, SZ}) = reduce($unicsymb(), x)
    end
end
for elem in unaryOps
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{1} end
        call(::$unicsymb, x) = $elem(x)
        $(elem){T, N, SZ}(x::AbstractFixedArray{T, N, SZ}) = map($unicsymb(), x)
    end
end
for elem in binaryOps2
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{2} end
        call(::$unicsymb, x, y) = $elem(x, y)
        $elem{T, N, SZ}(x::AbstractFixedArray{T, N, SZ}, y::AbstractFixedArray{T, N, SZ}) = map($unicsymb(), x, y)
    end
end
for elem in binaryOps
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{2} end
        call(::$unicsymb, x, y) = $elem(x, y)
        $elem{T, N, SZ}(x::AbstractFixedArray{T, N, SZ}, y::AbstractFixedArray{T, N, SZ}) = map($unicsymb(), x, y)
        $elem{T, N, SZ}(x::Real, y::AbstractFixedArray{T, N, SZ}) = map($unicsymb(), x, y)
        $elem{T, N, SZ}(x::AbstractFixedArray{T, N, SZ}, y::Real) = map($unicsymb(), x, y)
    end
end
end