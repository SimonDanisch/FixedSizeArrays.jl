using GeometryTypes

const v = (0,0,0)
const v2 = [0, 0, 0]

 function f(v, n)
    x = v[1:3]
    for i = 1:n
        for row=1:3
            x = v[1:3]
        end
    end
    return x
end
const N = 10^5
@time f(v, N)
@time f(v, N)
@time f(v, N)


# 518.823 milliseconds (3000 k allocations: 78125 KB, 3.03% gc time)
@time f(v2, N)
@time f(v2, N)
@time f(v2, N)

@code_llvm getindex(v, 1)
println("###########################")
@code_llvm getindex(v2, 1)