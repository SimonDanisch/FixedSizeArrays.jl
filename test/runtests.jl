tic()
using FixedSizeArrays
toc()
using Base.Test

# write your own tests here
@test 1 == 1


immutable RGB{T} <: FixedVector{T, 3}
    r::T
    g::T
    b::T
end
#Base.call{T <: RGB}(::Type{T}, a,b,c) = RGB(a,b,c)

immutable Vec3{T} <: FixedVector{T, 3}
    x::T
    y::T
    z::T
end
immutable Vec4{T} <: FixedVector{T, 4}
    x::T
    y::T
    z::T
    w::T
end
immutable Vec2{T} <: FixedVector{T, 2}
    x::T
    y::T
end

immutable Mat4x4{T} <: FixedMatrix{T, 4,4}
    c1::T
    c2::T
    c3::T
    c4::T
    c5::T
    c6::T
    c7::T
    c8::T
    c9::T
    c10::T
    c11::T
    c12::T
    c13::T
    c14::T
    c15::T
    c16::T
end

const a = RGB(1f0,2f0,3f0)
b 		= RGB(7f0,3f0,0f0)
const c = Vec3(7f0,3f0,0f0)
const ca = [7f0,3f0,0f0]
println(Main.Vec3)
println(c+c)
println(max(a,b))
println(maximum(a))
println(digamma(a))
println(rand(RGB{Float32}))

@show convert(Vec3{Float32}, a)
@show convert(RGB{Float32},  c)
@show convert(RGB{Float32},  Vec3{Float32}[c for i=1:10])


const ARR = [float32(x) for x=1:16]
@show const A = Mat4x4(ARR...)
@show const p = Vec4(1f0,2f0,4f0,1f0)
@show A[1,2]
@show A[Vec3(1,2,3)]
@show nvec(1,2,3,4)
@show A[1,:]
@show A[:,1]
amul = A*A
const ARRM = reshape(ARR, (4,4))
mul = ARRM*ARRM

for i=1:4, j=1:4
	println(mul[i,j], " == ", amul[i,j])
end

@show b[RGB(1,1,1)]

println(dot(c,c))
println(dot(ca,ca))