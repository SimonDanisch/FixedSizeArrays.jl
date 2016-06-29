module Test
using BenchmarkTools

using FixedSizeArrays

function match1(v)
    v1,v2,v3,v4 = v
    v1,v2,v3,v4
end

function match2(v)
    v1,v2,v3,v4 = Tuple(v)
    v1,v2,v3,v4
end
v = Vec{4,Float64}((2., 5., 77., 3.))
@show @benchmark match1($v)
@show @benchmark match2($v)

end