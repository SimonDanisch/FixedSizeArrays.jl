module Test1
tic()
using FixedSizeArrays
elapsed_time = toq()
println("loading time FixedSizeArrays: ", elapsed_time)

immutable RGB{T} <: AbstractFixedVector{T, 3}
  r::T
  g::T
  b::T
end

immutable RGBA{T} <: AbstractFixedVector{T, 4}
  r::T
  g::T
  b::T
  a::T
end

immutable Normal3{T} <: AbstractFixedVector{T, 3}
  x::T
  y::T
  z::T
end

immutable Vector4{T} <: AbstractFixedVector{T, 4}
  x::T
  y::T
  z::T
  w::T
end

immutable Point2{T} <: AbstractFixedVector{T, 2}
  x::T
  y::T
end

immutable TransformationMatrix{T} <: AbstractFixedMatrix{T, 4, 4}
  c1::Vector4{T}
  c2::Vector4{T}
  c3::Vector4{T}
  c4::Vector4{T}
end

const a = RGB(0.0,0.0,1.0)
const b = RGB(1.0,2.0,0.0)
@show a + b
@show cross(a,b)

@show sin(b)

const c = [Vector4(i,i,i,i) for i=0:0.01:20]
function test_speed(a,b,n)
  for i=1:n
    a+b
  end
end
println("Vector + Vector speed:")
@time test_speed(c[rand(1:length(c))], c[rand(1:length(c))], 10^6)
@time test_speed(c[rand(1:length(c))], c[rand(1:length(c))], 10^6)
println()
t = TransformationMatrix(c[1:4]...)
println()
@show row(t, 1)
println("Matrix + Matrix speed:")
@time test_speed(t, t, 10^6)
@time test_speed(t, t, 10^6)
println()

@show t+t
@show t.*t
@show t.*0.5

end

module Test2
tic()
using ImmutableArrays
elapsed_time = toq()
println("loading time ImmutableArrays: ", elapsed_time)

a = Vector3(0.0,0.0,1.0)
b = Vector3(1.0,2.0,0.0)
@show a + b
@show cross(a,b)

@show sin(b)

const c = [Vector4(i,i,i,i) for i=0:0.01:20]
function test_speed(a,b,n)
  for i=1:n
    a+b
  end
end

println("Vector + Vector speed:")
@time test_speed(c[rand(1:length(c))], c[rand(1:length(c))], 10^6)
@time test_speed(c[rand(1:length(c))], c[rand(1:length(c))], 10^6)
println()
t = Matrix4x4(c[1:4]...)
println()

println("Matrix + Matrix speed:")
@time test_speed(t, t, 10^6)
@time test_speed(t, t, 10^6)
println()

@show t+t
@show t.*t
@show t.*0.5

end