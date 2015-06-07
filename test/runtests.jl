using GeometryTypes
using Base.Test
const a = ["1.909", "1.909", "1.909"]
@test Vector3{Float64}(1.909) == Vector3{Float64}(a)