abstract Func{T}

abstract AbstractFixedArray{T,N,SZ}

immutable Vec3{T} <: AbstractFixedArray{T, 1, (3,)}
    x::T 
    y::T
    z::T
end
accessor_expr(variable, accessor_symbol) = parse("getfield($(variable), :$(accessor_symbol))") # Am I silly, or is there no other way to integrate a symbol?
Base.getindex(x::Vec3, i) = getfield(x, i)
stagedfunction map(f::Func{2}, a::Vec3, b::Vec3)
    accessors = names(a)
    expr = Any[]
    for elem in accessors
        gfa = accessor_expr("a", elem)
        gfb = accessor_expr("b", elem) 
        push!(expr, :(call(f, $gfa, $gfb)))
    end
    :($a($(expr...)))
end
function map2(f::Func{2}, a::AbstractFixedArray, b::AbstractFixedArray)

    Vec3(ntuple(x-> call(f, a[x], b[x]), 3)...)
end
immutable AddFun <: Func{2} end
call(::AddFun, x, y) = x + y


function test(n, x)
    for i=1:n
        map(AddFun(), x,x)
    end
end
function test2(n, x)
    for i=1:n
        map2(AddFun(), x,x)
    end
end
const a = Vec3(1,2,3)
@time test(1, a)
@time test(10^7, a)
@time test(10^7, a)
println("liiiiil")
@time test2(1, a)
@time test2(10^7, a)
@time test2(10^7, a)