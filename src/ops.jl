Base.promote_array_type{FSA <: FixedArray, T}(F, ::Type{T}, ::Type{FSA}) = FSA

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
const binaryOps = (:.+, :.-,:.*, :./, :.\, :.^, :+, :-,
                   :min, :max,
                   :div, :fld, :rem, :mod, :mod1, :cmp,
                   :atan2, :besselj, :bessely, :hankelh1, :hankelh2,
                   :besseli, :besselk, :beta, :lbeta)

# Pair by comparison of all elements and element wise comparisons
const comparisonOps = ((:(==), :(.==)), (:(!=), :(.!=)), (:(<), :(.<)), (:(<=),
                        :(.<=)), (:(>), :(.>)), (:(>=), :(.>=)))

const matrixOps = (:*, :/)

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
        @inline $op{T <: FixedArray}(x::T, y::T) = map($functor_name(), x, y)
        @inline $op{T1 <: FixedArray, T2 <: FixedArray}(x::T1, y::T2) = $op(promote(x, y)...)
        @inline $op{T <: Number}(x::T, y::FixedArray{T}) = map($functor_name(), x, y)
        @inline $op{T1 <: Number, T2}(x::T1, y::FixedArray{T2}) = $op(promote(x, y)...)
        @inline $op{T <: Number}(x::FixedArray{T}, y::T) = map($functor_name(), x, y)
        @inline $op{T1, T2 <: Number}(x::FixedArray{T1}, y::T2) = $op(promote(x, y)...)
    end)
end

for (op_all, op_el) in comparisonOps

    # Define element wise operations
    functor_name_el, functor_expr_el = gen_functor(op_el, 2)
    eval(quote
        $functor_expr_el
        @inline $op_el{T <: FixedArray}(x::T, y::T) = map($functor_name_el(), x, y)
        @inline $op_el{T1 <: FixedArray, T2 <: FixedArray}(x::T1, y::T2) = $op_el(promote(x, y)...)

        @inline $op_all{T1 <: FixedArray, T2 <: FixedArray}(x::T1, y::T2) = all($op_el(x, y))
    end)

end

for op in matrixOps
    functor_name, functor_expr = gen_functor(op, 2)
    eval(quote
        $functor_expr
        @inline $op{T <: Number}(x::T, y::FixedArray{T}) = map($functor_name(), x, y)
        @inline $op{T1 <: Number, T2}(x::T1, y::FixedArray{T2}) = $op(promote(x, y)...)
        @inline $op{T <: Number}(x::FixedArray{T}, y::T) = map($functor_name(), x, y)
        @inline $op{T1, T2 <: Number}(x::FixedArray{T1}, y::T2) = $op(promote(x, y)...)
    end)
end

@inline function promote{T1 <: FixedArray, T2 <: FixedArray}(a::T1, b::T2)
    T = promote_type(eltype(T1), eltype(T2))
    map(T, a), map(T, b)
end
function promote{T1, T2 <: Number}(a::FixedArray{T1}, b::T2)
    T = promote_type(T1, T2)
    map(T, a), T(b)
end
function promote{T1 <: Number, T2}(a::T1, b::FixedArray{T2})
    T = promote_type(T1, T2)
    T(a), map(T, b)
end

function promote_rule{N, T, X<:Number}(::Type{Vec{N,T}}, ::Type{X})
    Vec{N, promote_type(T, X)}
end

@inline ctranspose{R, C, T}(a::Mat{R, C, T}) = Mat(ntuple(CRowFunctor(a), Val{R}))
@generated function ctranspose{N,T}(b::Vec{N,T})
    expr = [:(b._[$i]',) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end
@inline transpose{R, C, T}(a::Mat{R, C, T}) = Mat(ntuple(RowFunctor(a), Val{R}))
@generated function transpose{N,T}(b::Vec{N,T})
    expr = [:(transpose(b._[$i]),) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end
@inline Base.hypot{T}(v::FixedVector{2,T}) = hypot(v[1],v[2])

immutable DotFunctor <: Func{2} end
call(::DotFunctor, a, b) = a'*b
@inline dot{T <:  FixedArray}(a::T, b::T) = sum(map(DotFunctor(), a, b))

immutable BilinearDotFunctor <: Func{2} end
call(::BilinearDotFunctor, a, b) = a*b
@inline bilindot{T <: Union{FixedArray, Tuple}}(a::T, b::T) = sum(map(DotFunctor(), a, b))

@inline bilindot{T}(a::NTuple{1,T}, b::NTuple{1,T}) = @inbounds return a[1]*b[1]
@inline bilindot{T}(a::NTuple{2,T}, b::NTuple{2,T}) = @inbounds return (a[1]*b[1] + a[2]*b[2])
@inline bilindot{T}(a::NTuple{3,T}, b::NTuple{3,T}) = @inbounds return (a[1]*b[1] + a[2]*b[2] + a[3]*b[3])
@inline bilindot{T}(a::NTuple{4,T}, b::NTuple{4,T}) = @inbounds return (a[1]*b[1] + a[2]*b[2] + a[3]*b[3]+a[4]*b[4])


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


trace(A::FixedMatrix{1,1}) = A[1,1]
trace(A::FixedMatrix{2,2}) = A[1,1] + A[2,2]
trace(A::FixedMatrix{3,3}) = A[1,1] + A[2,2] + A[3,3]
trace(A::FixedMatrix{4,4}) = A[1,1] + A[2,2] + A[3,3] + A[4,4]


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
(*){T, M, N, O}(a::FixedMatrix{M, N, T}, b::FixedVector{O, T}) = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))

@generated function *{T, N}(a::FixedVector{N, T}, b::FixedMatrix{1, N, T})
    expr = Expr(:tuple, [Expr(:tuple, [:(a[$i] * b[$j]) for i in 1:N]...) for j in 1:N]...)
    :( Mat($(expr)) )
end

@generated function *{T, M, N}(a::Mat{M, N, T}, b::Vec{N,T})
    expr = [:(bilindot(row(a, $i), b.(1))) for i=1:M]
    :( Vec($(expr...)) )
end
@generated function *{T, M, N, R}(a::Mat{M, N, T}, b::Mat{N, R, T})
    expr = Expr(:tuple, [Expr(:tuple, [:(bilindot(row(a, $i), column(b,$j))) for i in 1:M]...) for j in 1:R]...)
    :( Mat($(expr)) )
end


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

# To support @test_approx_eq
Base.Test.approx_full(a::FixedArray) = a
