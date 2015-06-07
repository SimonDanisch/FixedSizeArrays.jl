immutable FM4x4{T}
	x::NTuple{16, T}
end

immutable FM4x4S{T}
	x1::T
	x2::T
	x3::T
	x4::T
	x5::T
	x6::T
	x7::T
	x8::T
	x9::T
	x10::T
	x11::T
	x12::T
	x13::T
	x14::T
	x15::T
	x16::T
end
@inline Base.getindex{T}(a::FM4x4{T}, i::Int) = a.x[i]
@inline Base.getindex{T}(a::FM4x4S{T}, i::Int) = a.(i)

mdet(A) = @inbounds begin return (
        A[13] * A[10]  * A[7]  * A[4]  - A[9] * A[14] * A[7]  * A[4]   -
        A[13] * A[6]   * A[11] * A[4]  + A[5] * A[14] * A[11] * A[4]   +
        A[9]  * A[6]   * A[15] * A[4]  - A[5] * A[10] * A[15] * A[4]   -
        A[13] * A[10]  * A[3]  * A[8]  + A[9] * A[14] * A[3]  * A[8]   +
        A[13] * A[2]   * A[11] * A[8]  - A[1] * A[14] * A[11] * A[8]   -
        A[9]  * A[2]   * A[15] * A[8]  + A[1] * A[10] * A[15] * A[8]   +
        A[13] * A[6]   * A[3]  * A[12] - A[5] * A[14] * A[3]  * A[12]  -
        A[13] * A[2]   * A[7]  * A[12] + A[1] * A[14] * A[7]  * A[12]  +
        A[5]  * A[2]   * A[15] * A[12] - A[1] * A[6]  * A[15] * A[12]  -
        A[9]  * A[6]   * A[3]  * A[16] + A[5] * A[10] * A[3]  * A[16]  +
        A[9]  * A[2]   * A[7]  * A[16] - A[1] * A[10] * A[7]  * A[16]  -
        A[5]  * A[2]   * A[11] * A[16] + A[1] * A[6]  * A[11] * A[16])::Float64 end


const v2 = rand(4,4)
const v = FM4x4((v2...,))
const v3 = FM4x4S(v2...)

function f(v, n)
    x = 0.0
    for i = 1:n
        x += mdet(v)
    end
    return x
end
const N = 10^7
@time f(v, N)
@time f(v, N)
@time f(v, N)

# 518.823 milliseconds (3000 k allocations: 78125 KB, 3.03% gc time)
@time f(v2, N)
@time f(v2, N)
@time f(v2, N)

@time f(v3, N)
@time f(v3, N)
@time f(v3, N)

