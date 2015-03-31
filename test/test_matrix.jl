using FixedSizeArrays

type A{T} <: FixedVector{T, 3}
	x::T
	y::T 
	z::T
end
type B{T} <: FixedVector{T, 3}
	f::T
	d::T 
	g::T
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
const a = A(1,2,3)
const b = A(1,2,3)

@time a+a
@time a+a

@time b+b
@time b+b
const X = rand(4,4)
const Z = Mat4x4(X...)
@show sum(Z)
@time sum(Z)
@time sum(Z)

@time sum(X)
@time sum(X)
#elapsed time: 6.641e-6 seconds (92 bytes allocated)
#elapsed time: 7.245e-6 seconds (92 bytes allocated)
println(X)
println(Z)
println("################")
println(ctranspose(X))
println(ctranspose(Z))

@show Z*Z
@time Z*Z
@time Z*Z
@time Z*Z
@time Z*Z
@time Z*Z


const matmul = MatMulFunctor{Z, Z}
@time map(matmul, FixedMatrix{eltype(Z), 4, 4})
@time map(matmul, FixedMatrix{eltype(Z), 4, 4})
@time map(matmul, FixedMatrix{eltype(Z), 4, 4})
@time map(matmul, FixedMatrix{eltype(Z), 4, 4})
@code_native map(matmul, FixedMatrix{eltype(Z), 4, 4})
#elapsed time: 6.037e-6 seconds (224 bytes allocated)
#elapsed time: 7.848e-6 seconds (224 bytes allocated)
#elapsed time: 7.848e-6 seconds (224 bytes allocated)