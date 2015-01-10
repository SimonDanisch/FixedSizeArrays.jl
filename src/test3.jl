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



function generate_arrays(VectorName, MatrixName, maxSz::Integer)

    # expression functions
    vecTyp(n) = symbol(string(VectorName,n))
    vecTypT(n) = Expr(:curly, vecTyp(n), :T)
    matTyp(r,c) = symbol(string(MatrixName, r,"x",c))
    matTypT(r,c) = Expr(:curly, matTyp(r,c,), :T)
    elt(i) = symbol(string("e",i))
    col(i) = symbol(string("c",i))
    mem(s,e) = Expr(:.,s,Expr(:quote,e))
    velt(v,i) = mem(v,elt(i))
    melt(m,i,j) = mem(mem(m,col(j)),elt(i))

    # vector types
    for sz = 1:maxSz
        local Typ = vecTyp(sz)
        local TypT = vecTypT(sz)

        # the body of the type definition
        local defn = :(immutable $TypT <: AbstractFixedVector{T, $sz} end)

        # the members of the type
        for i = 1:sz
            local e = elt(i)
            push!(defn.args[3].args, :($e::T))
        end

        # instantiate the type definition
        eval(defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::T) = $TypT())
        for i = 1:sz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::T))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        eval(ctorn)
        eval(ctor1)
    end

    # matrix types
    for rSz = 1:maxSz, cSz = 1:maxSz
        local Typ = matTyp(rSz,cSz)
        local TypT = matTypT(rSz,cSz)
        local ColTyp = vecTyp(rSz)
        local ColTypT = vecTypT(rSz)

        # the body of the type definition
        local defn = :(immutable $TypT <: AbstractFixedMatrix{T, $rSz, $cSz} end)

        # the members of the type
        for i = 1:cSz
            local c = col(i)
            push!(defn.args[3].args, :($c::$ColTypT))
        end

        # instantiate the type definition
        eval(defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::$ColTypT) = $TypT())
        for i = 1:cSz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::$ColTypT))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        eval(ctorn)
        eval(ctor1)

        # construction from a scalar
        @eval $TypT(a::T) = $Typ($ColTyp(a))
    end
end

generate_arrays("Vector", "Matrix", 4)
typealias Vec2d Vector2{Float64}
typealias Vec3d Vector3{Float64}
typealias Vec4d Vector4{Float64}
typealias Vec3f Vector3{Float32}

v1 = Vec3d(1.0,2.0,3.0)
v2 = Vec3d(6.0,5.0,4.0)
# indexing
@assert v1[1] == 1.0
@assert v1[2] == 2.0
@assert v1[3] == 3.0
@assert try v1[-1]; false; catch e; isa(e,BoundsError); end
@assert try v1[0];  false; catch e; isa(e,BoundsError); end
@assert try v1[4];  false; catch e; isa(e,BoundsError); end

# negation
@assert -v1 == Vec3d(-1.0,-2.0,-3.0)
@assert isa(-v1,Vec3d)

# addition
@assert v1+v2 == Vec3d(7.0,7.0,7.0)

# subtraction
@assert v2-v1 == Vec3d(5.0,3.0,1.0)

# multiplication
@assert v1.*v2 == Vec3d(6.0,10.0,12.0)

# division
@assert v1 ./ v1 == Vec3d(1.0,1.0,1.0)

# scalar operations
@assert 1.0 + v1 == Vec3d(2.0,3.0,4.0)
@assert 1.0 .+ v1 == Vec3d(2.0,3.0,4.0)
@assert v1 + 1.0 == Vec3d(2.0,3.0,4.0)
@assert v1 .+ 1.0 == Vec3d(2.0,3.0,4.0)
@assert 1 + v1 == Vec3d(2.0,3.0,4.0)
@assert 1 .+ v1 == Vec3d(2.0,3.0,4.0)
@assert v1 + 1 == Vec3d(2.0,3.0,4.0)
@assert v1 .+ 1 == Vec3d(2.0,3.0,4.0)

@assert v1 - 1.0 == Vec3d(0.0,1.0,2.0)
@assert v1 .- 1.0 == Vec3d(0.0,1.0,2.0)
@assert 1.0 - v1 == Vec3d(0.0,-1.0,-2.0)
@assert 1.0 .- v1 == Vec3d(0.0,-1.0,-2.0)
@assert v1 - 1 == Vec3d(0.0,1.0,2.0)
@assert v1 .- 1 == Vec3d(0.0,1.0,2.0)
@assert 1 - v1 == Vec3d(0.0,-1.0,-2.0)
@assert 1 .- v1 == Vec3d(0.0,-1.0,-2.0)

@assert 2.0 * v1 == Vec3d(2.0,4.0,6.0)
@assert 2.0 .* v1 == Vec3d(2.0,4.0,6.0)
@assert v1 * 2.0 == Vec3d(2.0,4.0,6.0)
@assert v1 .* 2.0 == Vec3d(2.0,4.0,6.0)
@assert 2 * v1 == Vec3d(2.0,4.0,6.0)
@assert 2 .* v1 == Vec3d(2.0,4.0,6.0)
@assert v1 * 2 == Vec3d(2.0,4.0,6.0)
@assert v1 .* 2 == Vec3d(2.0,4.0,6.0)

@assert v1 / 2.0 == Vec3d(0.5,1.0,1.5)
@assert v1 ./ 2.0 == Vec3d(0.5,1.0,1.5)
@assert v1 / 2 == Vec3d(0.5,1.0,1.5)
@assert v1 ./ 2 == Vec3d(0.5,1.0,1.5)

@assert 12.0 ./ v1 == Vec3d(12.0,6.0,4.0)
@assert 12 ./ v1 == Vec3d(12.0,6.0,4.0)

@assert v1.^2 == Vec3d(1.0,4.0,9.0)
@assert v1.^2.0 == Vec3d(1.0,4.0,9.0)
@assert 2.0.^v1 == Vec3d(2.0,4.0,8.0)
@assert 2.^v1 == Vec3d(2.0,4.0,8.0)

# vector norm
@assert norm(Vec3d(1.0,2.0,2.0)) == 3.0

# cross product
@assert cross(v1,v2) == Vec3d(-7.0,14.0,-7.0)
@assert isa(cross(v1,v2),Vec3d)
x = Vec2d(1,0)
y = Vec2d(0,1)
@assert cross(x,y) == 1.0
@assert cross(y,x) == -1.0


#=
# basis vectors
e1 = unit(Vec4d,1)
e2 = unit(Vec4d,2)
e3 = unit(Vec4d,3)
e4 = unit(Vec4d,4)
@assert e1 == Vec4d(1.0,0.0,0.0,0.0)
@assert e2 == Vec4d(0.0,1.0,0.0,0.0)
@assert e3 == Vec4d(0.0,0.0,1.0,0.0)
@assert e4 == Vec4d(0.0,0.0,0.0,1.0)
=#

# type conversion
#=
@assert isa(convert(Vec3f,v1),Vec3f)
@assert Vector3([1.0,2.0,3.0]) == v1
@assert convert(Vec3d,[1.0,2.0,3.0]) == v1
@assert isa(convert(Vector{Float64},v1),Vector{Float64})
@assert convert(Vector{Float64},v1) == [1.0,2.0,3.0]
=#

# matrix operations

typealias Mat1d Matrix1x1{Float64}
typealias Mat2d Matrix2x2{Float64}
typealias Mat3d Matrix3x3{Float64}
typealias Mat4d Matrix4x4{Float64}

@assert zero(Mat2d) == Mat2d(Vec2d(0.0,0.0),Vec2d(0.0,0.0))
@assert zero(Vec4d) == row(zero(Mat4d), 1)

#=
v = Vec4d(1.0,2.0,3.0,4.0)
r = row(v)
c = column(v)

@assert c*r == Mat4d(Vec4d(1.0,2.0,3.0,4.0),
                     Vec4d(2.0,4.0,6.0,8.0),
                     Vec4d(3.0,6.0,9.0,12.0),
                     Vec4d(4.0,8.0,12.0,16.0))
@assert r*c == Matrix1x1(30.0)
@assert r' == c
@assert c' == r
@assert row(r,1) == v
@assert column(c,1) == v
@assert row(r+c',1) == 2*v
@assert sum(r) == sum(v)
@assert prod(c) == prod(v)
@assert eye(Mat3d) == Mat3d(Vec3d(1.0,0.0,0.0),Vec3d(0.0,1.0,0.0),Vec3d(0.0,0.0,1.0))
@assert v*eye(Mat4d)*v == 30.0
@assert -r == -1.0*r
@assert diag(diagm(v)) == v

# type conversion
@assert isa(convert(Matrix1x4{Float32},r),Matrix1x4{Float32})
jm = rand(4,4)
im = Matrix4x4(jm)
@assert isa(im,Mat4d)
@assert jm == im
im = convert(Mat4d,jm)
@assert isa(im,Mat4d)
@assert jm == im
jm2 = convert(Array{Float64,2},im)
@assert isa(jm2,Array{Float64,2})
@assert jm == jm2
=#
