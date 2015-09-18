# operations
const unaryOps = (:-, :~, :conj, :abs,
                  :sin, :cos, :tan, :sinh, :cosh, :tanh,
                  :asin, :acos, :atan, :asinh, :acosh, :atanh,
                  :sec, :csc, :cot, :asec, :acsc, :acot,
                  :sech, :csch, :coth, :asech, :acsch, :acoth,
                  :sinc, :cosc, :cosd, :cotd, :cscd, :secd,
                  :sind, :tand, :acosd, :acotd, :acscd, :asecd,
                  :asind, :atand, :rad2deg, :deg2rad,
                  :log, :log2, :log10, :log1p, :exponent, :exp,
                  :exp2, :expm1, :cbrt, :sqrt, :erf,
                  :erfc, :erfcx, :erfi, :dawson, :ceil, :floor,
                  :trunc, :round, :significand, :lgamma,
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

const reductions = ((:sum,:+), (:prod,:*), (:minimum,:min), (:maximum,:max))

function gen_functor(func::Symbol, unary::Int)
    functor_name  = gensym()
    arguments     = ntuple(i->symbol("arg$i"), unary)
    functor_expr  = quote
        immutable $functor_name <: Func{$unary} end
        @inline call(::$functor_name, $(arguments...)) = $func($(arguments...))
    end
    return (functor_name, functor_expr)
end

for (callfun, reducefun) in reductions
    functor_name, functor_expr = gen_functor(reducefun, 2)
    eval(quote
        $functor_expr
        @inline $(callfun){T <: FixedArray}(x::T) = reduce($functor_name(), x)
    end)
end
for op in unaryOps
    functor_name, functor_expr = gen_functor(op, 1)
    eval(quote
        $functor_expr
        @inline $(op){T <: FixedArray}(x::T) = map($functor_name(), x)
    end)
end

for op in binaryOps
    functor_name, functor_expr = gen_functor(op, 2)
    eval(quote
        $functor_expr
        @inline $op{T1 <: FixedArray, T2 <: FixedArray}(x::T1, y::T2) = $op(promote(x, y)...)
        @inline $op{T1 <: FixedArray}(x::Number, y::T1) = eltype(y) == typeof(x) ? map($functor_name(), x, y) : $op(promote(x, y)...)
        @inline $op{T1 <: FixedArray}(x::T1, y::Number) = eltype(x) == typeof(y) ? map($functor_name(), x, y) : $op(promote(x, y)...)

        @inline $op{T <: FixedArray}(x::T,    y::T) = map($functor_name(), x, y)
    end)
end

@inline function promote{T1 <: FixedArray, T2 <: FixedArray}(a::T1, b::T2)
    T = promote_type(eltype(T1), eltype(T2))
    map(T, a), map(T, b)
end
@inline function promote{T1 <: FixedArray}(a::T1, b::Number)
    T = promote_type(eltype(T1), typeof(b))
    map(T, a), T(b)
end
@inline function promote{T1 <: FixedArray}(a::Number, b::T1)
    T = promote_type(typeof(a), eltype(T1))
    T(a), map(T, b)
end
@inline ctranspose{R, C, T}(a::Mat{R, C, T}) = Mat(ntuple(CRowFunctor(a), Val{R}))

@inline Base.hypot{T}(v::FixedVector{2,T}) = hypot(v[1],v[2])

immutable DotFunctor <: Func{2} end
call(::DotFunctor, a, b) = a'*b
@inline dot{T <: Union(FixedArray, Tuple)}(a::T, b::T) = sum(map(DotFunctor(), a, b))

@inline dot{T}(a::NTuple{1,T}, b::NTuple{1,T}) = @inbounds return a[1]'*b[1]
@inline dot{T}(a::NTuple{2,T}, b::NTuple{2,T}) = @inbounds return (a[1]'*b[1] + a[2]'*b[2])
@inline dot{T}(a::NTuple{3,T}, b::NTuple{3,T}) = @inbounds return (a[1]'*b[1] + a[2]'*b[2] + a[3]'*b[3])
@inline dot{T}(a::NTuple{4,T}, b::NTuple{4,T}) = @inbounds return (a[1]'*b[1] + a[2]'*b[2] + a[3]'*b[3]+a[4]'*b[4])

#cross{T}(a::FixedVector{2, T}, b::FixedVector{2, T}) = a[1]*b[2]-a[2]*b[1] # not really used!?
@inline cross{T<:Real}(a::FixedVector{3, T}, b::FixedVector{3, T}) = @inbounds return typeof(a)(
    a[2]*b[3]-a[3]*b[2],
    a[3]*b[1]-a[1]*b[3],
    a[1]*b[2]-a[2]*b[1]
)

@inline norm{T, N}(a::FixedVector{T, N})     = sqrt(dot(a,a))
@inline normalize{FSA <: FixedArray}(a::FSA) = a / norm(a)

#Matrix
@inline det{T}(A::FixedMatrix{1, 1, T}) = @inbounds return ( A[1] )
@inline det{T}(A::FixedMatrix{2, 2, T}) = @inbounds return ( A[1,1]*A[2,2] - A[1,2]*A[2,1])
@inline det{T}(A::FixedMatrix{3, 3, T}) = @inbounds return (
    A[1,1]*(A[2,2]*A[3,3]-A[2,3]*A[3,2]) - 
    A[1,2]*(A[2,1]*A[3,3]-A[2,3]*A[3,1]) + 
    A[1,3]*(A[2,1]*A[3,2]-A[2,2]*A[3,1])
)
@inline det{T}(A::FixedMatrix{4, 4, T}) = @inbounds return (
    A[13] * A[10]  * A[7]  * A[4]  - A[9] * A[14] * A[7]  * A[4]   -
    A[13] * A[6]   * A[11] * A[4]  + A[5] * A[14] * A[11] * A[4]   +
    A[9]  * A[6]   * A[15] * A[4]  - A[5] * A[10] * A[15] * A[4]   -
    A[13] * A[10]  * A[3]  * A[8]  + A[9] * A[14] * A[3]  * A[8]   +
    A[13] * A[2]   * A[11] * A[8]  - A[1] * A[14] * A[11] * A[8]   -
    A[9]  * A[2]   * A[15] * A[8]  + A[1] * A[10] * A[15] * A[8]   +
    A[13] * A[6]   * A[3]  * A[12] - A[5] * A[14] * A[3]  * A[12]  -
    A[13] * A[2]   * A[7]  * A[12] + A[1] * A[14] * A[7]  * A[12]  +
    A[5]  * A[2]   * A[15] * A[12] - A[1] * A[6]  * A[15] * A[12]  -
    A[9]  * A[6]   * A[3]  * A[16] + A[5] * A[10] * A[3]  * A[16]  +
    A[9]  * A[2]   * A[7]  * A[16] - A[1] * A[10] * A[7]  * A[16]  -
    A[5]  * A[2]   * A[11] * A[16] + A[1] * A[6]  * A[11] * A[16]
)


@inline inv{T}(A::Mat{1, 1, T}) = @inbounds return Mat{1, 1, T}(inv(A[1]))
@inline function inv{T}(A::Mat{2, 2, T})
    determinant = det(A)
    @inbounds return Mat{2, 2, T}(
        (A[2,2] /determinant, -A[2,1]/determinant),
        (-A[1,2]/determinant, A[1,1] /determinant)
    )
end
@inline function inv{T}(A::Mat{3, 3, T})
    determinant = det(A)
    @inbounds return  Mat{3, 3, T}(
        ((A[2,2]*A[3,3]-A[2,3]*A[3,2]) /determinant,
        -(A[2,1]*A[3,3]-A[2,3]*A[3,1])/determinant,
        (A[2,1]*A[3,2]-A[2,2]*A[3,1]) /determinant),

        (-(A[1,2]*A[3,3]-A[1,3]*A[3,2])/determinant,
        (A[1,1]*A[3,3]-A[1,3]*A[3,1]) /determinant,
        -(A[1,1]*A[3,2]-A[1,2]*A[3,1])/determinant),

        ((A[1,2]*A[2,3]-A[1,3]*A[2,2]) /determinant,
        -(A[1,1]*A[2,3]-A[1,3]*A[2,1])/determinant,
        (A[1,1]*A[2,2]-A[1,2]*A[2,1]) /determinant)
    )
end


@inline function inv{T}(A::Mat{4, 4, T})
    determinant = det(A)
    @inbounds return Mat{4, 4, T}(
        ((A[2,3]*A[3,4]*A[4,2] - A[2,4]*A[3,3]*A[4,2] + A[2,4]*A[3,2]*A[4,3] - A[2,2]*A[3,4]*A[4,3] - A[2,3]*A[3,2]*A[4,4] + A[2,2]*A[3,3]*A[4,4]) / determinant,
        (A[2,4]*A[3,3]*A[4,1] - A[2,3]*A[3,4]*A[4,1] - A[2,4]*A[3,1]*A[4,3] + A[2,1]*A[3,4]*A[4,3] + A[2,3]*A[3,1]*A[4,4] - A[2,1]*A[3,3]*A[4,4]) / determinant,
        (A[2,2]*A[3,4]*A[4,1] - A[2,4]*A[3,2]*A[4,1] + A[2,4]*A[3,1]*A[4,2] - A[2,1]*A[3,4]*A[4,2] - A[2,2]*A[3,1]*A[4,4] + A[2,1]*A[3,2]*A[4,4]) / determinant,
        (A[2,3]*A[3,2]*A[4,1] - A[2,2]*A[3,3]*A[4,1] - A[2,3]*A[3,1]*A[4,2] + A[2,1]*A[3,3]*A[4,2] + A[2,2]*A[3,1]*A[4,3] - A[2,1]*A[3,2]*A[4,3]) / determinant),

        ((A[1,4]*A[3,3]*A[4,2] - A[1,3]*A[3,4]*A[4,2] - A[1,4]*A[3,2]*A[4,3] + A[1,2]*A[3,4]*A[4,3] + A[1,3]*A[3,2]*A[4,4] - A[1,2]*A[3,3]*A[4,4]) / determinant,
        (A[1,3]*A[3,4]*A[4,1] - A[1,4]*A[3,3]*A[4,1] + A[1,4]*A[3,1]*A[4,3] - A[1,1]*A[3,4]*A[4,3] - A[1,3]*A[3,1]*A[4,4] + A[1,1]*A[3,3]*A[4,4]) / determinant,
        (A[1,4]*A[3,2]*A[4,1] - A[1,2]*A[3,4]*A[4,1] - A[1,4]*A[3,1]*A[4,2] + A[1,1]*A[3,4]*A[4,2] + A[1,2]*A[3,1]*A[4,4] - A[1,1]*A[3,2]*A[4,4]) / determinant,
        (A[1,2]*A[3,3]*A[4,1] - A[1,3]*A[3,2]*A[4,1] + A[1,3]*A[3,1]*A[4,2] - A[1,1]*A[3,3]*A[4,2] - A[1,2]*A[3,1]*A[4,3] + A[1,1]*A[3,2]*A[4,3]) / determinant),


        ((A[1,3]*A[2,4]*A[4,2] - A[1,4]*A[2,3]*A[4,2] + A[1,4]*A[2,2]*A[4,3] - A[1,2]*A[2,4]*A[4,3] - A[1,3]*A[2,2]*A[4,4] + A[1,2]*A[2,3]*A[4,4]) / determinant,
        (A[1,4]*A[2,3]*A[4,1] - A[1,3]*A[2,4]*A[4,1] - A[1,4]*A[2,1]*A[4,3] + A[1,1]*A[2,4]*A[4,3] + A[1,3]*A[2,1]*A[4,4] - A[1,1]*A[2,3]*A[4,4]) / determinant,
        (A[1,2]*A[2,4]*A[4,1] - A[1,4]*A[2,2]*A[4,1] + A[1,4]*A[2,1]*A[4,2] - A[1,1]*A[2,4]*A[4,2] - A[1,2]*A[2,1]*A[4,4] + A[1,1]*A[2,2]*A[4,4]) / determinant,
        (A[1,3]*A[2,2]*A[4,1] - A[1,2]*A[2,3]*A[4,1] - A[1,3]*A[2,1]*A[4,2] + A[1,1]*A[2,3]*A[4,2] + A[1,2]*A[2,1]*A[4,3] - A[1,1]*A[2,2]*A[4,3]) / determinant),


        ((A[1,4]*A[2,3]*A[3,2] - A[1,3]*A[2,4]*A[3,2] - A[1,4]*A[2,2]*A[3,3] + A[1,2]*A[2,4]*A[3,3] + A[1,3]*A[2,2]*A[3,4] - A[1,2]*A[2,3]*A[3,4]) / determinant,
        (A[1,3]*A[2,4]*A[3,1] - A[1,4]*A[2,3]*A[3,1] + A[1,4]*A[2,1]*A[3,3] - A[1,1]*A[2,4]*A[3,3] - A[1,3]*A[2,1]*A[3,4] + A[1,1]*A[2,3]*A[3,4]) / determinant,
        (A[1,4]*A[2,2]*A[3,1] - A[1,2]*A[2,4]*A[3,1] - A[1,4]*A[2,1]*A[3,2] + A[1,1]*A[2,4]*A[3,2] + A[1,2]*A[2,1]*A[3,4] - A[1,1]*A[2,2]*A[3,4]) / determinant,
        (A[1,2]*A[2,3]*A[3,1] - A[1,3]*A[2,2]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,1]*A[2,2]*A[3,3]) / determinant)
    )
end


# Matrix
(*){T, M, N, O, K}(a::FixedMatrix{M, N, T}, b::FixedMatrix{O, K, T}) = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))


@generated function *{T, M, N}(a::Mat{M, N, T}, b::Vec{N,T})
   expr = [:(dot(row(a, $i), b.(1))) for i=1:M]
   return quote 
       $(Expr(:boundscheck, false))
       Vec($(expr...))
   end
end
@generated function *{T, M, N, R}(a::Mat{M, N, T}, b::Mat{N, R, T})
   expr = Expr(:tuple, [Expr(:tuple, [:(dot(row(a, $i), column(b,$j))) for i in 1:M]...) for j in 1:R]...)
   return quote 
       $(Expr(:boundscheck, false))
       Mat($(expr))
   end
end


@generated function (*){T, FSV <: FixedVector, C}(a::FSV, b::Mat{1, C, T})
    N = length(a)
    N != C && throw(DimensionMismatch("DimensionMissmatch: $N != $R for $(typeof(a)), $(typeof(b))"))
    expr = [:(tuple($([:(a[$i]*b[$j]) for j=1:C]...))) for i=1:C]
    :($(Expr(:boundscheck, false)); Mat(tuple($(expr...))))
end

(*){FSV <: FixedVector}(a::FSV, b::FSV) = Mat{1, 1, eltype(FSV)}(dot(a,b))

function (==)(a::FixedVectorNoTuple, b::FixedVectorNoTuple)
    s_a = size(a)
    s_b = size(b)
    s_a == s_b || return false
    @inbounds for i = 1:length(a)
        a[i] == b[i] || return false
    end
    true
end
(==)(a::FixedArray, b::FixedArray) = a.(1) == b.(1)

(==){R, T, FSA <: FixedVector}(a::FSA, b::Mat{R, 1, T}) = a.(1) == column(b,1)
(==){R, T, FSA <: FixedVector}(a::Mat{R, 1, T}, b::FSA) = column(a,1) == b.(1)
function (==)(a::FixedArray, b::AbstractArray)
    s_a = size(a)
    s_b = size(b)
    s_a == s_b || return false
    @inbounds for i = 1:length(a)
        a[i] == b[i] || return false
    end
    true
end

(==)(a::AbstractArray, b::FixedArray) = b == a


