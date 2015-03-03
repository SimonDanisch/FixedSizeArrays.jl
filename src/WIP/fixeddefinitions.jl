module FixedSizeArrays
importall Base

export FixedArray
export FixedVector
export FixedMatrix


abstract FixedArray{T,N,SZ}
typealias FixedVector{T, C} FixedArray{T, 1, (C,)}
typealias FixedMatrix{T, Row, Column} FixedArray{T, 2, (Row, Column)}
importall Base

eltype{T,N,SZ}(A::FixedArray{T,N,SZ})           = T
length{T,N,SZ}(A::FixedArray{T,N,SZ})           = prod(SZ)
ndims{T,N,SZ}(A::FixedArray{T,N,SZ})            = N
size{T,N,SZ}(A::FixedArray{T,N,SZ})             = SZ
size{T,N,SZ}(A::FixedArray{T,N,SZ}, d::Integer) = SZ[d]

# Ugly workaround for not having triangular dispatch:
eltype{T <: FixedVector}(A::Type{T})                = A.types[1]
eltype{T <: FixedMatrix}(A::Type{T})                = eltype(A.types[1])

length{T <: FixedVector}(A::Type{T})                = length(A.types)
length{T <: FixedMatrix}(A::Type{T})                = prod(size(A))

ndims{T <: FixedVector}(A::Type{T})                 = 1
ndims{T <: FixedMatrix}(A::Type{T})                 = 2

size{T <: FixedVector}(A::Type{T})                  = (length(A),)
size{T <: FixedVector}(A::Type{T}, d::Integer)      = (length(A),) # should throw an error!?
size{T <: FixedMatrix}(A::Type{T})                  = (length(A.types), length(A.types[1]))
size{T <: FixedMatrix}(A::Type{T}, d::Integer)      = size(A)[d] 

getindex{T,C}(A::FixedVector{T, C}, i::Integer)                 = getfield(A, i)
getindex{T,M,N}(A::FixedMatrix{T, M,N}, i::Integer, j::Integer) = getfield(getfield(A, i),j)

columntype{T,N,SZ}(x::FixedMatrix{T,N,SZ}) = typeof(x[1])
columntype{T <: FixedMatrix}(x::Type{T})   = first(x.types)
unit{T,C}(v::FixedVector{T,C}) = v/norm(v)
cross{T}(a::FixedVector{T, 2}, b::FixedVector{T, 2}) = a[1]*b[2]-a[2]*b[1]

cross{T}(a::FixedVector{T, 3}, b::FixedVector{T, 3}) = typeof(a)(a[2]*b[3]-a[3]*b[2], 
                                                     a[3]*b[1]-a[1]*b[3], 
                                                     a[1]*b[2]-a[2]*b[1])

stagedfunction zero{T <: FixedArray}(::Type{T}) 
    ttypes = T.types
    if !all(isleaftype, ttypes)
        error("please provide concrete type. Types given: ", accessors)
    end
    :($T($(map(x->:(zero($x)), ttypes)...)))
end

dot{T,C}(v1::FixedVector{T,C}, v2::FixedVector{T,C}) = sum(v1.*conj(v2))
norm{T,C}(v::FixedVector{T,C}) = sqrt(dot(v,v))
function norm{T,C}(v::FixedVector{T,C}, p::Number)
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

column{T, Cardinality}(v::FixedVector{T, Cardinality}) = v
column{T, Column, Row}(v::FixedMatrix{T, Row, Column}, i::Integer) = v[i]
row{T, Cardinality}(v::FixedVector{T, Cardinality}, i::Integer) = v[i]
stagedfunction row{T, Column, Row}(v::FixedMatrix{T, Row, Column}, i::Integer)
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

reduce{T}(f::Func{2}, a::FixedVector{T, 1}) = a
stagedfunction reduce(f::Func{2}, a::FixedVector)
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
stagedfunction reduce(f::Func{2}, a::FixedMatrix)
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

stagedfunction map(f::Func{1}, a::FixedVector)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa)))
    end
    :($a($(expr...)))
end
stagedfunction map{T <: FixedVector}(f::Func{2}, a::T, b::T)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa, $gfb)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{2}, a::Real, b::FixedVector)
    accessors = names(b)
    expr = Any[]
    for elem in accessors
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, a, $gfb)))
    end
    :($b($(expr...)))
end
stagedfunction map(f::Func{2}, a::FixedVector, b::Real)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(call(f, $gfa, b)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{1}, a::FixedMatrix)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(map(f, $gfa)))
    end
    :($a($(expr...)))
end
stagedfunction map{T <: FixedMatrix}(f::Func{2}, a::T, b::T)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem) # Am I silly, or is there no other way to integrate a symbol?
        gfb = accessor_expr("b", elem) # Am I silly, or is there no other way to integrate a symbol?
        push!(expr, :(map(f, $gfa, $gfb)))
    end
    :($a($(expr...)))
end
stagedfunction map(f::Func{2}, a::Real, b::FixedMatrix)
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
        $(callfun){T, N, SZ}(x::FixedArray{T, N, SZ}) = reduce($unicsymb(), x)
    end
end
for elem in unaryOps
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{1} end
        call(::$unicsymb, x) = $elem(x)
        $(elem){T, N, SZ}(x::FixedArray{T, N, SZ}) = map($unicsymb(), x)
    end
end
for elem in binaryOps2
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{2} end
        call(::$unicsymb, x, y) = $elem(x, y)
        $elem{T, N, SZ}(x::FixedArray{T, N, SZ}, y::FixedArray{T, N, SZ}) = map($unicsymb(), x, y)
    end
end
for elem in binaryOps
    unicsymb = gensym()
    @eval begin 
        immutable $unicsymb <: Func{2} end
        call(::$unicsymb, x, y) = $elem(x, y)
        $elem{T, N, SZ}(x::FixedArray{T, N, SZ}, y::FixedArray{T, N, SZ}) = map($unicsymb(), x, y)
        $elem{T, N, SZ}(x::Real, y::FixedArray{T, N, SZ}) = map($unicsymb(), x, y)
        $elem{T, N, SZ}(x::FixedArray{T, N, SZ}, y::Real) = map($unicsymb(), x, y)
    end
end
end