function Base.promote_array_type{FSA <: FixedArray, N, T<:Number}(
        F, ::Type{T}, ::Type{Array{FSA, N}}
    )
    FSA
end
function Base.promote_array_type{FSA <: FixedArray, N, T<:Number}(
        F, ::Type{Array{FSA, N}}, ::Type{T}
    )
    FSA
end
function Base.promote_array_type{FSA <: FixedArray, T<:Number}(
        F, ::Type{T}, ::Type{FSA}
    )
    FSA
end
function Base.promote_array_type{FSA <: FixedArray, T<:Number}(
        F, ::Type{FSA}, ::Type{T}
    )
    FSA
end


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
                  :eta, :zeta, :digamma, :real, :imag)

# vec-vec and vec-scalar
const binaryOps = (:.+, :.-,:.*, :./, :.\, :.^,
                   :.==, :.!=, :.<, :.<=, :.>, :.>=, :+, :-,
                   :min, :max,
                   :div, :fld, :rem, :mod, :mod1, :cmp,
                   :atan2, :besselj, :bessely, :hankelh1, :hankelh2,
                   :besseli, :besselk, :beta, :lbeta)

const reductions = ((:sum,:+), (:prod,:*), (:minimum,:min), (:maximum,:max))

function gen_functor(func::Symbol, unary::Int)
    functor_name  = gensym()
    arguments     = ntuple(i->Symbol("arg$i"), unary)
    functor_expr  = quote
        immutable $functor_name <: Functor{$unary} end
        @compat @inline (::$functor_name)($(arguments...)) = $func($(arguments...))
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
        @inline $op{F1<:FixedArray, F2<:FixedArray}(x::F1, y::F2) = map($functor_name(), x, y)
        @inline $op{F<:FixedArray}(x::F, y::Number) = map($functor_name(), x, y)
        @inline $op{F<:FixedArray}(x::Number, y::F) = map($functor_name(), x, y)
        @inline $op{T,N}(x::FixedArray{T,N}, y::AbstractArray{T,N}) = map($functor_name(), x, y)
        @inline $op{T,N}(x::AbstractArray{T,N}, y::FixedArray{T,N}) = map($functor_name(), x, y)
    end)
end

# Non-dotted * and / of FSAs and scalars
mul_func_name, mul_func_expr = gen_functor(:*, 2)
div_func_name, div_func_expr = gen_functor(:/, 2)
@eval begin
    $div_func_expr
    $mul_func_expr
    @inline *(x::Number, y::FixedArray) = map($mul_func_name(), x, y)
    @inline *(x::FixedArray, y::Number) = map($mul_func_name(), x, y)
    @inline /(x::FixedArray, y::Number) = map($div_func_name(), x, y)
    function /(x::Number, y::FixedArray)
        Base.depwarn("/(x::Number, y::FixedArray) had unclear semantics and will be removed", :/)
        map($div_func_name(), x, y)
    end
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

@inline ctranspose{R, C, T}(a::Mat{R, C, T}) = Mat{C,R,T}(ntuple(CRowFunctor(a), Val{R}))
@generated function ctranspose{N,T}(b::Vec{N,T})
    expr = [:(Tuple(b)[$i]',) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end
@inline transpose{R, C, T}(a::Mat{R, C, T}) = Mat(ntuple(RowFunctor(a), Val{R}))
@generated function transpose{N,T}(b::Vec{N,T})
    expr = [:(transpose(Tuple(b)[$i]),) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end
@inline Base.hypot{T}(v::FixedVector{2,T}) = hypot(v[1],v[2])

immutable DotFunctor <: Functor{2} end
@compat (::DotFunctor)(a, b) = a'*b
@inline dot(a::FixedVector, b::FixedVector) = sum(map(DotFunctor(), a, b))

immutable BilinearDotFunctor <: Functor{2} end
@compat (::BilinearDotFunctor)(a, b) = a*b
@inline bilindot(a::FixedVector, b::FixedVector) = sum(map(BilinearDotFunctor(), a, b))


function cross{T}(a::FixedVector{2, T}, b::FixedVector{2, T})
    return a[1]*b[2]-a[2]*b[1]
end
@inline function cross{T1, T2}(a::FixedVector{3, T1}, b::FixedVector{3, T2})
    @inbounds elements = (a[2]*b[3]-a[3]*b[2],
                          a[3]*b[1]-a[1]*b[3],
                          a[1]*b[2]-a[2]*b[1])
    construct_similar(typeof(a), elements)
end

@inline norm{N, T}(a::FixedVector{N, T})     = sqrt(dot(a,a))
@inline normalize{FSA <: FixedArray}(a::FSA) = a / norm(a)

function norm{N, T}(a::FixedVector{N, T}, p)
    isinf(p) && return maxabs(a)
    ret = abs(a[1])^p
    for k = 2:N
        ret += abs(a[k])^p
    end
    ret^(1/p)
end

function Base.isnan(p::FixedArray)
    for elem in p
        isnan(elem) && return true
    end
    false
end

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
det(A::FixedMatrix) = det(Matrix(A))


trace(A::FixedMatrix{1,1}) = A[1,1]
trace(A::FixedMatrix{2,2}) = A[1,1] + A[2,2]
trace(A::FixedMatrix{3,3}) = A[1,1] + A[2,2] + A[3,3]
trace(A::FixedMatrix{4,4}) = A[1,1] + A[2,2] + A[3,3] + A[4,4]
trace(A::FixedMatrix) = trace(Matrix(A))

\{m,n,T1,T2}(mat::Mat{m,n,T1}, v::Vec{n,T2}) = inv(mat)*v

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

inv(A::FixedMatrix) = typeof(A)(inv(Matrix(A)))

lyap{T}(a::Mat{1, 1, T}, c::Mat{1, 1, T}) = Mat{1,1,T}(lyap(a[1,1],c[1,1]))
function lyap{T}(a::Mat{2, 2, T}, c::Mat{2, 2, T})
    d = det(a)
    t = trace(a)
     -(d*c  + (a - t*I)*c*(a-t*I)')/(2*d*t) # http://www.nber.org/papers/w8956.pdf
end
lyap{m,T}(a::Mat{m,m,T},c::Mat{m,m,T}) = Mat(lyap(Matrix(a),Matrix(c)))

chol{T<:Base.LinAlg.BlasFloat}(m::Mat{1, 1, T}) = Mat{1,1}(chol(m[1,1]))
function chol{T<:Base.LinAlg.BlasFloat}(m::Mat{2,2,T})
    m[1,2]==m[2,1]' || error("Matrix not symmetric")
    l11 = chol(m[1,1])
    @inbounds return Mat{2, 2, T}(
        (l11, zero(T)),
        (inv(l11)*m[1,2], chol(m[2,2] - m[2,1]*inv(m[1,1])*m[1,2]))
    )
end
chol{n,T<:Base.LinAlg.BlasFloat}(m::Mat{n,n,T}) = Mat{n,n,T}(full(Base.LinAlg.chol!(Matrix(m))))
chol!(m::Mat, ::Type{UpperTriangular}) = chol(m) # for 0.4
_chol!(m::Mat, ::Type{UpperTriangular}) = chol(m) # for 0.5
chol!(m::Mat, ::Type{Val{:U}}) = chol!(m, UpperTriangular) # for pre-0.5

# Matrix products
# General shape mismatched versions are errors
(*){T1, T2, M, N, O, K}(a::FixedMatrix{M, N, T1}, b::FixedMatrix{O, K, T2}) = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))
(*){T1, T2, M, N, O}(a::FixedMatrix{M, N, T1}, b::FixedVector{O, T2}) = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))
# vector * (row vector)
@generated function *{T1, T2, N, M}(a::FixedVector{N, T1}, b::FixedMatrix{1,M,T2})
    elements = Expr(:tuple, [Expr(:tuple, [:(a[$i] * b[$j]) for i in 1:N]...) for j in 1:M]...)
    :(construct_similar($b, $elements))
end
# matrix * vector
@generated function *{T1, T2, M, N}(a::FixedMatrix{M, N, T1}, b::FixedVector{N, T2})
    total_terms = M*N
    if total_terms <= 64
        # Full unrolling
        elements = Expr(:tuple, [Expr(:call, :+, [:(a[$i,$k]*b[$k]) for k = 1:N]...) for i in 1:M]...)
    else
        # Expand as a bunch of dot products
        elements = Expr(:tuple, [:(bilindot(Vec(row(a,$i)),b)) for i in 1:M]...)
    end
    :(construct_similar($b, $elements))
end
function *(a::AbstractMatrix, b::FixedVector)
    a*Vector(b)
end
@generated function *{T1, T2, M, N}(a::FixedMatrix{M, N, T1}, b::AbstractVector{T2})
    quote
        length(b) == $N || throw(DimensionMismatch("$b is wrong size - expecting vector of length $N"))
        return a*Vec{N,T2}(b)
    end
end
# matrix * matrix
@generated function *{T1, T2, M, N, P}(a::FixedMatrix{M, N, T1}, b::FixedMatrix{N, P, T2})
    total_terms = M*N*P
    if total_terms <= 64  # 4x4 * 4x4
        # Full unrolling
        elements = Expr(:tuple, [Expr(:tuple, [Expr(:call, :+, [:(a[$i,$k]*b[$k,$j]) for k = 1:N]...) for i in 1:M]...) for j in 1:P]...)
    else
        # Expand as a bunch of mat*vec expressions
        elements = Expr(:tuple, [:(Tuple(a*Vec(column(b,$j)))) for j = 1:P]...)
    end
    :(construct_similar($a, $elements))
end
function *(a::AbstractMatrix, b::FixedMatrix)
    a*Matrix(b)
end
function *(a::FixedMatrix, b::AbstractMatrix)
    Matrix(a)*b
end


# Predicates
function (==)(a::FixedVectorNoTuple, b::FixedVectorNoTuple)
    s_a = size(a)
    s_b = size(b)
    s_a == s_b || return false
    @inbounds for i = 1:length(a)
        a[i] == b[i] || return false
    end
    true
end
(==)(a::FixedArray, b::FixedArray) = Tuple(a) == Tuple(b)

(==){R, T, FSA <: FixedVector}(a::FSA, b::Mat{R, 1, T}) = Tuple(a) == column(b, 1)
(==){R, T, FSA <: FixedVector}(a::Mat{R, 1, T}, b::FSA) = column(a,1) == Tuple(b)
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

isinteger{FSA<:FixedArray}(A::FSA) = all(isinteger, A)
isinteger{T<:Integer}(A::FixedArray{T}) = true

isreal{FSA<:FixedArray}(A::FSA) = all(isreal, A)
isreal{T<:Real}(A::FixedArray{T}) = true

ishermitian{FSA<:FixedMatrix}(A::FSA) = (A == A')



# To support @test_approx_eq
Base.Test.approx_full(a::FixedArray) = a

# UniformScaling

*(J::Base.LinAlg.UniformScaling, A::FixedArray) = J.λ*A
*(A::FixedArray, J::Base.LinAlg.UniformScaling) = A*J.λ
/(A::FixedArray, J::Base.LinAlg.UniformScaling) = A/J.λ

+{m, n, T}(A::Mat{m,n, T}, J::Base.LinAlg.UniformScaling) = A + J.λ*eye(Mat{m,n,T})
+{m, n, T}(J::Base.LinAlg.UniformScaling, A::Mat{m,n, T}) = A + J
-{m, n, T}(A::Mat{m,n, T}, J::Base.LinAlg.UniformScaling) = A + (-J)
-{m, n, T}(J::Base.LinAlg.UniformScaling, A::Mat{m,n, T}) = J.λ*eye(Mat{m,n,T}) - A

# typed floor, ceil, round, trunc
Base.floor{T}(::Type{T}, A::FixedArray) = map(x->floor(T, x), A)
Base.ceil{T}( ::Type{T}, A::FixedArray) = map(x->ceil( T, x), A)
Base.trunc{T}(::Type{T}, A::FixedArray) = map(x->trunc(T, x), A)
Base.round{T}(::Type{T}, A::FixedArray) = map(x->round(T, x), A)

import Base.Cartesian.@nif
shift(t::Tuple) = Base.tail(t)
unshift(t::Tuple, val) = tuple(val, t...)
push(t::Tuple, val) = tuple(t..., val)
@generated function pop{N,T}(t::NTuple{N,T})
    Expr(:tuple, [:(t[$i]) for i in 1:N-1]...)
end

@generated function deleteat{N,T,i}(x::NTuple{N,T}, ::Type{Val{i}})
    (1 <= i <= N) || throw(BoundsError((x,Val{i})))
    Expr(:tuple, [:(x[$j]) for j in deleteat!([1:N...], i)]...)
end
@generated function deleteat{N,T}(x::NTuple{N,T}, i::Int)
    quote
        (1 <= i <= N) || throw(BoundsError(x,i))
        @nif $N d->(i==d) d-> deleteat(x, Val{d})
    end
end

@generated function insert{N,T,i}(t::NTuple{N,T}, ::Type{Val{i}}, item)
    (1 <= i <= N+1) || throw(BoundsError())
    args = Any[:(t[$k]) for k in 1:N]
    insert!(args, i, :item)
    Expr(:tuple, args...)
end
@generated function insert{N,T}(x::NTuple{N,T}, i::Int, item)
    quote
        (1 <= i <= N+1) || throw(BoundsError())
        @nif $(N+1) d->(i==d) d-> insert(x, Val{d}, item)
    end
end
insert(v::FixedVector, i, val) = construct_similar(typeof(v), insert(Tuple(v), i, val))

for f in (:shift, :pop)
    eval(quote
        ($f)(v::FixedVector) = construct_similar(typeof(v), ($f)(Tuple(v)))
    end)
end
for f in (:unshift, :deleteat, :push)
    eval(quote
        ($f)(v::FixedVector, val) = construct_similar(typeof(v), ($f)(Tuple(v), val))
    end)
end

function reverse{T, N}(x::FixedVector{N, T})
    map(typeof(x)) do i
        x[N-i+1]
    end
end
