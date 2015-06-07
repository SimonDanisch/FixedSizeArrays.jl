using GeometryTypes

const v = Vector3{Int}(0,0,0)
const v2 = [0, 0, 0]

 function f(v, n)
    x = 0
    for i = 1:n
        for row=1:3
            x += v[row] # even v.x is still slower than juliaarray[row]
        end
    end
    return x
end
const N = 1_000_0000
@time f(v, N)
@time f(v, N)
@time f(v, N)


# 518.823 milliseconds (3000 k allocations: 78125 KB, 3.03% gc time)
@time f(v2, N)
@time f(v2, N)
@time f(v2, N)
