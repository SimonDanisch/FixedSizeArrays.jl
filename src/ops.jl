include("core.jl")
include("staged.jl")
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

const reductions = ((:sum,:+),(:prod,:*),(:minimum,:min),(:maximum,:max))

function gen_functor(func::Symbol, unary::Int)
    functor_name  = gensym()
    arguments     = ntuple(i->symbol("arg$i"), unary)
    functor_expr  = quote 
        immutable $functor_name <: Func{$unary} end
        call(::$functor_name, $(arguments...)) = $func($(arguments...))
    end
    return (functor_name, functor_expr)
end

for (callfun, reducefun) in reductions
    functor_name, functor_expr = gen_functor(reducefun, 2)
    eval(quote 
        $functor_expr
        $(callfun){T, N, SZ}(x::AbstractFixedArray{T, N, SZ}) = reduce($functor_name(), x)
    end)
end
for op in unaryOps
    functor_name, functor_expr = gen_functor(op, 1)
    eval(quote 
        $functor_expr
        $(op){T, N, SZ}(x::AbstractFixedArray{T, N, SZ}) = map($functor_name(), x)
    end)
end

for op in binaryOps
    functor_name, functor_expr = gen_functor(op, 2)
    eval(quote 
        $functor_expr
        $op{T, N, SZ}(x::AbstractFixedArray{T, N, SZ}, y::AbstractFixedArray{T, N, SZ}) = map($functor_name(), x, y)
        $op{T, N, SZ}(x::Real,                         y::AbstractFixedArray{T, N, SZ}) = map($functor_name(), x, y)
        $op{T, N, SZ}(x::AbstractFixedArray{T, N, SZ}, y::Real)                         = map($functor_name(), x, y)
    end)
end

immutable RandFunc <: Func{1} end
call(::Type{RandFunc}, x) = rand(x)
rand{T <: AbstractFixedArray}(x::Type{T}) =  T([rand(eltype(x)) for i=1:length(x)]...)

function convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::T2)
    @assert sizeof(a) == sizeof(b) "Type $a ($(sizeof(a))) doesn't have the same bit size as type $b ($(sizeof(b)))"
    reinterpret(a, [b])[1]
end


dot(a::AbstractFixedVector, b::AbstractFixedVector) = sum(a.*b)
function convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(b) % sizeof(a) == 0 "Type $a ($(sizeof(a))) doesn't have the same bit size as type $b ($(sizeof(b)))"
    println(a)
    reinterpret(a, b, (div(sizeof(b), sizeof(a)),))
end


# Matrix
stagedfunction (*){T, M, N, K}(a::AbstractFixedMatrix{T, M, N}, b::AbstractFixedMatrix{T, N, K})
    :(AbstractFixedMatrix{$T, $M, $K}( 
         $([:(dot(a[$i, :], b[:, $j])) for i=1:M, j=1:K]...)
    ))
end


immutable RGB{T} <: AbstractFixedVector{T, 3}
    r::T
    g::T
    b::T
end
immutable Vec3{T} <: AbstractFixedVector{T, 3}
    x::T
    y::T
    z::T
end
immutable Vec4{T} <: AbstractFixedVector{T, 4}
    x::T
    y::T
    z::T
    w::T
end
immutable Vec2{T} <: AbstractFixedVector{T, 2}
    x::T
    y::T
end

immutable Mat4x4{T} <: AbstractFixedMatrix{T, 4,4}
    c1::T
    c2::T
    c3::T
    c4::T
    c5::T
    c6::T
    c7::T
    c8::T
    c9::T
    c10::T
    c11::T
    c12::T
    c13::T
    c14::T
    c15::T
    c16::T
end

dot(a::RGB, b::RGB) = a.r*b.r + a.g*b.g + a.b*b.b


const a = RGB(1f0,2f0,3f0)
b = RGB(7f0,3f0,0f0)
const c = Vec3(7f0,3f0,0f0)
println(a+a)
println(max(a,b))
println(maximum(a))
println(digamma(a))
println(rand(RGB{Float32}))

@show convert(Vec3{Float32}, a)
@show convert(RGB{Float32},  c)
@show convert(RGB{Float32},  Vec3{Float32}[c for i=1:10])

@show const A = Mat4x4(ntuple(x->float32(x), 16)...)
@show const p = Vec4(1f0,2f0,4f0,1f0)
@show A[1,2]
@show A*A

@show b[RGB(1,1,1)]