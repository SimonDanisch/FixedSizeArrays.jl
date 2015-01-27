importall Base
abstract AbstractFixedArray1
abstract AbstractFixedArray2

getindex(A::AbstractFixedArray1,    i::Integer)         = A.(i) 
getindex(A::AbstractFixedArray2,    i::Integer)         = A.(i) 
length(A::AbstractFixedArray1) = 4
length(A::AbstractFixedArray2) = 4
length{T <: AbstractFixedArray1}(A::Type{T}) = 4
length{T <: AbstractFixedArray2}(A::Type{T}) = 4
accessor{T <: AbstractFixedArray1}(arg::Type{T}, i::Integer, name::Symbol) = :($name[$i])
accessor{T                      }(arg::Type{T}, i::Integer, name::Symbol) = :($name)
stagedfunction map{F <: Base.Func{2}, T <: AbstractFixedArray1}(f::Type{F}, t1::T, t2::T)
    expr = ntuple(i -> :( f( $(accessor(t1, i, :t1)), $(accessor(t2, i, :t2))) ), length(T))
    :($t1($(expr...)))
end
function map{F <: Base.Func{2}, T <: AbstractFixedArray2}(f::Type{F}, t1::T, t2::T)
    T(ntuple(i->f(t1[i], t2[i]), length(t1))...)
end

function gen_functor(func::Symbol, unary::Int)
    functor_name  = gensym()
    arguments     = ntuple(i->symbol("arg$i"), unary)
    functor_expr  = quote 
        immutable $functor_name <: Base.Func{$unary} end
        call(::Type{$functor_name}, $(arguments...)) = $func($(arguments...))
    end
    return (functor_name, functor_expr)
end

for op in [:+]
    functor_name, functor_expr = gen_functor(op, 2)
    eval(quote 
        $functor_expr
        $op(x::AbstractFixedArray1, y::AbstractFixedArray1) = map($functor_name, x, y)
        $op(x::AbstractFixedArray2, y::AbstractFixedArray2) = map($functor_name, x, y)
    end)
end
immutable F1{T} <: AbstractFixedArray1
  x::T
  y::T
  z::T
  w::T
end
immutable F2{T} <: AbstractFixedArray2
  x::T
  y::T
  z::T
  w::T
end
function test(a, N)
  result = a
  for i=1:N
    result = result + a
  end
end
const a1 = F1(1,2,3,4)
const a2 = F2(1,2,3,4)
const N = 10^7
@time test(a1, 1)
@time test(a1, N)

@time test(a2, 1)
@time test(a2, N)