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
dot(a::AbstractFixedArray, b::AbstractFixedArray) = sum(a.*b)

immutable RandFunc <: Func{1} end
call(::Type{RandFunc}, x) = rand(x)
rand{T <: AbstractFixedArray}(x::Type{T}) =  T([rand(eltype(x)) for i=1:length(x)]...)

function convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::T2)
    @assert sizeof(a) == sizeof(b) "Type $a ($(sizeof(a))) doesn't have the same bit size as type $b ($(sizeof(b)))"
    reinterpret(a, [b])[1] # why doesn't this work like reinterpret(Float32, Int32)
end
function convert{T1 <: AbstractFixedArray, T2 <: AbstractFixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(b) % sizeof(a) == 0 "Type $a, with size: ($(sizeof(a))) doesn't fit into the array b: $(length(b)) x $(sizeof(eltype(b)))"
    reinterpret(a, b, (div(sizeof(b), sizeof(a)),))
end







