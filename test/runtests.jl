module FSAtesting

using FixedSizeArrays
using FactCheck, Base.Test
using Compat

import FixedSizeArrays: similar_type

immutable Normal{N, T} <: FixedVector{N, T}
    values::NTuple{N, T}
end
immutable D3{N1, N2, N3, T} <: FixedArray{T, 3, Tuple{N1, N2, N3}}
    values::NTuple{N1, NTuple{N2, NTuple{N3, T}}}
end
immutable RGB{T} <: FixedVectorNoTuple{3, T}
    r::T
    g::T
    b::T
end

# subtyping:
immutable TestType{N,T} <: FixedVector{N,T}
    values::NTuple{N,T}
end

# Custom FSA with non-parameterized size and eltype
immutable Coord2D <: FixedVectorNoTuple{2,Float64}
    x::Float64
    y::Float64
end


typealias Vec1d Vec{1, Float64}
typealias Vec2d Vec{2, Float64}
typealias Vec3d Vec{3, Float64}
typealias Vec4d Vec{4, Float64}
typealias Vec3f Vec{3, Float32}
typealias Mat2d Mat{2,2, Float64}
typealias Mat3d Mat{3,3, Float64}
typealias Mat4d Mat{4,4, Float64}

# Compatibility hacks for 0.5 APL-style array slicing
if VERSION < v"0.5.0-dev+1195"
    # Remove all dimensions of size 1
    compatsqueeze(A) = squeeze(A,(find(collect(size(A)).==1)...))
else
    compatsqueeze(A) = A
end
function test()
facts("FixedSizeArrays") do

include("typeinf.jl")


context("fsa macro") do
    a = 1
    a1 = @fsa([a,2,3])
    a2 = @fsa([a 2 3])
    a3 = @fsa([a;2;3])
    a4 = @fsa([a 2;3 4])
    a5 = @fsa([a 2 3;4 5 6])
    a6 = @fsa([a 2;3 4;5 6])

    @fact a1 --> Vec(a,2,3)
    @fact a2 --> Mat((a,),(2,),(3,))
    @fact a3 --> Mat(((a,2,3),))
    @fact a4 --> Mat(((a,3),(2,4)))
    @fact a5 --> Mat(((a,4),(2,5),(3,6)))
    @fact a6 --> Mat(((a,3,5),(2,4,6)))
end

context("core") do
    context("ndims") do
        @fact ndims(D3) --> 3
        @fact ndims(Mat) --> 2
        @fact ndims(Vec) --> 1
        @fact ndims(Vec(1,2,3)) --> 1

        @fact ndims(D3{3,3,3}) --> 3
        @fact ndims(Mat{3,3}) --> 2
        @fact ndims(Vec{3}) --> 1

        @fact ndims(D3{3,3,3,Int}) --> 3
        @fact ndims(Mat{3,3,Int}) --> 2
        @fact ndims(Vec{3,Int}) --> 1
    end
    context("size_or") do
        @fact size_or(Mat, nothing) --> nothing
        @fact size_or(Mat{4}, nothing) --> nothing
        @fact size_or(Mat{4,4}, nothing) --> (4,4)
        @fact size_or(Mat{4,4, Float32}, nothing) --> (4,4)

        @fact size_or(Vec, nothing) --> nothing
        @fact size_or(Vec{4}, nothing) --> (4,)
        @fact size_or(Vec{4,Float32}, nothing) --> (4,)
        @fact size_or(FixedArray, nothing) --> nothing

    end
    context("eltype_or") do
        @fact eltype_or(Mat, nothing) --> nothing
        @fact eltype_or(Mat{4}, nothing) --> nothing
        @fact eltype_or(Mat{4,4}, nothing) --> nothing
        @fact eltype_or(Mat{4,4, Float32}, nothing) --> Float32

        @fact eltype_or(Vec, nothing) --> nothing
        @fact eltype_or(Vec{4}, nothing) --> nothing
        @fact eltype_or(Vec{4,Float32}, nothing) --> Float32

        @fact eltype_or(FixedArray, nothing) --> nothing

    end
    context("ndims_or") do
        @fact ndims_or(Mat, nothing) --> 2
        @fact ndims_or(Mat{4}, nothing) --> 2
        @fact ndims_or(Mat{4,4}, nothing) --> 2
        @fact ndims_or(Mat{4,4, Float32}, nothing) --> 2

        @fact ndims_or(Vec, nothing) --> 1
        @fact ndims_or(Vec{4}, nothing) --> 1
        @fact ndims_or(Vec{4, Float64}, nothing) --> 1

        @fact ndims_or(FixedArray, nothing) --> nothing
    end

    context("similar_type") do
        @fact similar_type(Vec{3,Int}, Float32) --> Vec{3, Float32}
        @fact similar_type(Vec{3}, Float32) --> Vec{3, Float32}
        @fact similar_type(Vec, Float32, (3,)) --> Vec{3, Float32}
        @fact similar_type(Vec, Float32, (1,2)) --> Mat{1,2, Float32}

        @fact similar_type(RGB, Float32) --> RGB{Float32}
        @fact similar_type(RGB{Float32}, Int) --> RGB{Int}
        @fact similar_type(RGB{Float32}, Int, (3,)) --> RGB{Int}
        @fact similar_type(RGB{Float32}, Int, (2,2)) --> Mat{2,2,Int}

        @fact similar_type(Mat{3,3,Int}, Float32) --> Mat{3,3,Float32}
        @fact similar_type(Mat, Float32, (3,3))   --> Mat{3,3,Float32}
        @fact similar_type(Mat{2,2,Int}, (3,3))   --> Mat{3,3,Int}

        @fact similar_type(Coord2D, Float64, (2,)) --> Coord2D
        @fact similar_type(Coord2D, Int, (2,))     --> Vec{2,Int}
        @fact similar_type(Coord2D, Float64, (3,)) --> Vec{3,Float64}
    end

    context("construct_similar") do
        @fact construct_similar(Vec{3,Int}, (1.0f0,2)) --> exactly(Vec{2,Float32}(1,2))
        @fact construct_similar(Vec{2}, (1,2,3))       --> exactly(Vec{3,Int}(1,2,3))
        @fact construct_similar(Vec, (1.0,2))          --> exactly(Vec{2,Float64}(1,2))

        @fact construct_similar(RGB, (1,2,3))                --> exactly(RGB{Int}(1,2,3))
        @fact construct_similar(RGB{Float32}, (1.0,2.0,3.0)) --> exactly(RGB{Float64}(1.0,2.0,3.0))

        @fact construct_similar(Mat{3,3,Int}, ((1.0f0,2),(1.0,2))) --> exactly(Mat{2,2,Float64}((1,2),(1,2)))
        @fact construct_similar(Mat, ((1,2),))                     --> exactly(Mat{2,1,Int}(((1,2),)))
    end

    context("nan") do
        for (p, r) in (
                (Point{2, Float32}(NaN, 1), true),
                (Point{2, Float64}(1, NaN), true),
                (Vec{11, Float64}(NaN), true),
                (Point{2, Float32}(1, 1), false),
                (RGB{Float32}(NaN), true),
            )
            @fact isnan(p) --> r
        end
    end

end


context("Array of FixedArrays") do

    N = 100
    a = Point{3, Float32}[Point{3, Float32}(0.7132) for i=1:N]
    b = RGB{Float32}[RGB{Float32}(52.293) for i=1:N]

    c = Point{3, Float64}[Point{3, Float64}(typemin(Float64)), a..., Point{3, Float64}(typemax(Float64))]
    d = RGB{Float64}[RGB(typemin(Float64)), b..., RGB(typemax(Float64))]

    context("reduce") do
        sa = sum(a)
        ma = mean(a)
        sb = sum(b)
        mb = mean(b)
        for i=1:3
            @fact sa[i]  --> roughly(Float32(0.7132*N))
            @fact ma[i]  --> roughly(Float32(0.7132*N)/ N)

            @fact sb[i]  --> roughly(Float32(52.293*N))
            @fact mb[i]  --> roughly(Float32(52.293*N)/ N)
        end

        @fact maximum(c) --> Point{3, Float32}(typemax(Float64))
        @fact minimum(c) --> Point{3, Float32}(typemin(Float64))

        @fact maximum(d) --> RGB(typemax(Float64))
        @fact minimum(d) --> RGB(typemin(Float64))

        @fact extrema(c) --> (minimum(c), maximum(c))
    end

    context("array ops") do
        for op in (.+, .-,.*, ./, .\, +, -)
            @fact typeof(op(a, 1f0)) --> typeof(a)
            @fact typeof(op(1f0, a)) --> typeof(a)
        end

        af = a + 1f0
        bf = b + 1f0
        aff = a + Point{3, Float32}(1)
        bff = b + RGB{Float32}(1)
        afd = a .+ 1f0
        bfd = b .+ 1f0
        @inferred(b .* 1f0)
        for i=1:N
            @fact a[1] + 1f0 --> af[i]
            @fact b[1] + 1f0 --> bf[i]
            @fact a[1] + 1f0 --> aff[i]
            @fact b[1] + 1f0 --> bff[i]
            @fact a[1] + 1f0 --> afd[i]
            @fact b[1] + 1f0 --> bfd[i]
        end
    end
    context("Show") do
        m1 = rand(Mat4d, 2)
        m2 = rand(RGB{Float32}, 2)
        m3 = rand(Vec3f, 2)
        println(m1)
        println(m2)
        println(m3)
        showcompact(Point(1,2,3))
    end
end


context("Constructor FixedVectorNoTuple") do
    for T=[Float32, Float64, Int, UInt, UInt32, UInt8]
        context("$T") do
            r = rand(T)
            x = RGB{Int}[RGB(1) for i=1:10]
            @fact RGB{Float32}(["0.222", "9.8822", "29.999"]) --> RGB{Float32}(0.222, 9.8822, 29.999)
            @fact RGB(["0.222", "9.8822", "29.999"]) --> RGB{Float64}(0.222, 9.8822, 29.999)
            @fact typeof(map(RGB{Float32}, x))  --> Vector{RGB{Float32}}
            @fact RGB{T}(r)                     --> RGB(r,r,r)
            @fact RGB{T}(Vec(r,r,r))            --> RGB(r,r,r)
            @fact RGB{T}([r,r,r])               --> RGB(r,r,r)
            @fact length(RGB{T}([r,r,r]))       --> 3
            @fact length(RGB{T})                --> 3
            @fact eltype(RGB{T}([r,r,r]))       --> T
            @fact eltype(RGB{T})                --> T
            @fact typeof(RGB(r,r,r))            --> RGB{T}
            @fact typeof(RGB{T}(1))             --> RGB{T}
            @fact typeof(RGB{T}(1,2,3))         --> RGB{T}
            @fact ndims(RGB{T}(1,2,3))          --> 1

            @fact typeof(RGB{T}(1f0))           --> RGB{T}
            @fact typeof(RGB{T}(1f0,2f0,3f0))   --> RGB{T}
            @fact typeof(RGB{T}(1f0, 2, 3.0))   --> RGB{T}
            @fact typeof(RGB(1f0, 2, 3.0))      --> RGB{Float64}
            @fact typeof(RGB{Int}(1f0, 2, 3.0)) --> RGB{Int}
            @fact_throws DimensionMismatch RGB((1,2,3), (2,3,4))
        end
    end
end

# A little brutal, but hey.... Better redudantant tests, than not enough tests
context("Constructor ") do
    context("Rand") do
        #Win32 seems to fail for rand(Vec4d)
        @fact typeof(rand(Vec4d)) --> Vec4d
        @fact typeof(rand(Mat4d)) --> Mat4d

        @fact typeof(rand(Mat{4,2, Int})) --> Mat{4,2, Int}
        @fact typeof(rand(Vec{7, Int})) --> Vec{7, Int}
        @fact typeof(rand(Vec{7, Int}, 1:7)) --> Vec{7, Int}
        @fact typeof(rand(Mat4d, -20f0:0.192f0:230f0)) --> Mat4d
        @fact typeof(rand(Mat{4,21,Float32}, -20f0:0.192f0:230f0)) --> Mat{4,21,Float32}

        x = rand(D3{4,4,4, Float32})
        @fact typeof(x) --> D3{4,4,4, Float32}
        @fact eltype(x) --> Float32
        @fact size(x) --> (4,4,4)
        @fact typeof(rand(Vec4d, 5,5)) --> Matrix{Vec4d}
    end
    context("Randn") do
        @fact typeof(randn(Base.Random.GLOBAL_RNG, Vec4d)) --> Vec4d
        @fact typeof(randn(Vec4d)) --> Vec4d
        @fact typeof(randn(Mat4d)) --> Mat4d
        @fact typeof(randn(Mat{4,2, Complex{Float64}})) --> Mat{4,2, Complex{Float64}}
        @fact typeof(randn(Vec{7, Complex{Float64}})) --> Vec{7, Complex{Float64}}
    end
    context("Zero") do
        @fact typeof(zero(Vec4d)) --> Vec4d
        @fact typeof(zero(Mat4d)) --> Mat4d

        @fact typeof(zero(Mat{4,2, Int})) --> Mat{4,2, Int}
        @fact typeof(zero(Vec{7, Int})) --> Vec{7, Int}
        @fact zero(Vec((1,2))) --> Vec((0,0))
        @fact zero(Vec((1.0,2.0))) --> Vec((0.0,0.0))
    end

    context("eye") do
        @fact typeof(eye(Mat4d)) --> Mat4d
        @fact typeof(eye(Mat{4,2, Int})) --> Mat{4,2, Int}
    end
    context("one") do
        x = one(Mat{4,2, Int})
        @fact typeof(one(Mat4d)) --> Mat4d
        @fact typeof(x) --> Mat{4,2, Int}
        @fact all(x-> x==1, x) --> true
    end

    context("unit") do
        u4 = unit(Vec4d, 1)
        u7 = unit(Vec{7, Int}, 7)
        @fact typeof(u4) --> Vec4d
        @fact typeof(u7) --> Vec{7, Int}
        @fact u4[1] --> 1.0
        @fact u4[2:end] --> (0.,0.,0.)

        @fact u7[end] --> 1
        @fact u7[1:end-1] --> (0,0,0,0,0,0)
    end
    for N=(1,10)
        context("construction, conversion, $N") do
            for VT=[Point, Vec], VT2=[Normal, Vec], ET=[Float32, Int, UInt], ET2=[Float64, UInt, Float32]
                rand_range  = ET(1):ET(10)
                rand_range2 = ET2(1):ET2(10)
                rn = rand(rand_range, N)
                v0 = VT(rn)
                # parse constructor:
                @fact VT{N, ET}(map(string, rn)) --> v0
                # multi constructor
                v1 = VT{N, ET}(rn...)
                @fact v1 --> v0
                @fact typeof(v1) --> VT{N, ET}
                @fact length(v1) --> N
                @fact eltype(v1) --> ET
                @fact ndims(v1) --> 1

                @fact length(typeof(v1)) --> N
                @fact eltype(typeof(v1)) --> ET

                for i=1:N
                    @fact v1[i] --> rn[i]
                end
                # from other FSA without parameters
                v2 = VT2(v1)

                @fact typeof(v2) --> VT2{N, ET}
                @fact length(v2) --> N
                @fact eltype(v2) --> ET
                for i=1:N
                    @fact v2[i] --> v1[i]
                end
                # from other FSA with parameters
                for i=1:N
                    v3 = VT2{i, ET2}(v1)
                    @fact typeof(v3) --> VT2{i, ET2}
                    @fact length(v3) --> i
                    @fact eltype(v3) --> ET2
                    for i=1:i
                        @fact v3[i] --> ET2(v2[i])
                    end
                end
                # from single
                r  = rand(rand_range)
                r2 = rand(rand_range2)
                v1 = VT{N, ET}(r)
                v2 = VT{N, ET2}(r)
                v3 = VT{N, ET}(r2)
                v4 = VT{N, ET2}(r2)

                for i=1:N
                    @fact v1[i] --> r
                    @fact v2[i] --> ET2(r)
                    @fact v3[i] --> r2
                    @fact v4[i] --> ET2(r2)
                end
                x = VT{N, ET}[VT{N, ET}(1) for i=1:10]
                x1 = VT2{N, ET}[VT{N, ET}(1) for i=1:10]
                x2 = map(VT2, x)
                x3 = map(VT, x2)
                @fact typeof(x)  --> Vector{VT{N, ET}}
                @fact typeof(x1) --> Vector{VT2{N, ET}}
                @fact typeof(x2) --> Vector{VT2{N, ET}}
                @fact typeof(x3) --> Vector{VT{N, ET}}
                @fact x3         --> x

                # Construction with only N, issue #56
                @fact VT{N}(ET(1)) --> Vec{N, ET}(1)
                @fact VT{N}(ntuple(x->ET(1), N)...) --> Vec{N, ET}(1)
            end
        end
    end
end


context("Constructors") do
    context("FixedVector: unary, from FixedVector") do
        @fact typeof(Vec3f(1,1,1))     --> Vec{3, Float32}
        @fact typeof(Vec3f(1,1f0,1))   --> Vec{3, Float32}
        @fact typeof(Vec3f(1f0,1,1.0)) --> Vec{3, Float32}
        @fact eltype(Vec3f(1f0,1,1.0)) --> Float32

        @fact typeof(Vec3f(1))      --> Vec{3, Float32}
        @fact typeof(Vec3f(0))      --> Vec{3, Float32}
        @fact Vec3f(1.0)             --> Vec(1f0,1f0,1f0)
        @fact Vec3f(1.0f0)             --> Vec(1f0,1f0,1f0)
        @fact Vec3f(1.0f0)             --> Vec3f(1)
        @fact Vec(1.0, 1.0, 1.0)     --> Vec3d(1)
        @fact Vec2d(Vec3d(1))         --> Vec(1.0, 1.0)
        @fact Vec(Vec3d(1), 1.0)     --> Vec4d(1)
        @fact Vec(Vec3d(1), 1)         --> Vec4d(1)
        @fact Vec3d(Vec3f(1.0))     --> Vec3d(1.0)
    end
end


context("map") do
    context("Vec and AbstractVector") do
        # Unary, binary & ternary map with specified output type
        @fact map(-, Vec{3,Float64}, Vec(1,2,3)) --> exactly(Vec{3,Float64}(-1,-2,-3))
        @fact map(+, Vec{3,Float64}, [1,2,3], Vec(1,2,3)) --> exactly(Vec{3,Float64}(2,4,6))
        @fact map(+, Vec{3,Float64}, [1,2,3], Vec(1,2,3), 1:3) --> exactly(Vec{3,Float64}(3,6,9))

        # Unary and binary map with deduced output types
        @fact map(-, Vec(1,2,3)) --> exactly(Vec{3,Int}(-1,-2,-3))
        @fact map(+, Vec(1,2,3), [1,2,3]) --> exactly(Vec{3,Int}(2,4,6))
        @fact map(+, [1,2,3], Vec(1,2,3)) --> exactly(Vec{3,Int}(2,4,6))
        @fact map(+, Vec(1,2,3), Vec(1,2,3)) --> exactly(Vec{3,Int}(2,4,6))
        # Some other `AbstractArray`s
        @fact map(+, Vec(1,2,3), 1:3) --> exactly(Vec{3,Int}(2,4,6))
        @fact map(+, 1:3, Vec(1,2,3)) --> exactly(Vec{3,Int}(2,4,6))

        # Binary map with mixed types
        @fact map(>, Vec(0.0,2.0), Vec(1,1)) --> exactly(Vec{2,Bool}(false,true))
        @fact map(+, Vec(0.0,0.0), Vec(1,1)) --> exactly(Vec{2,Float64}(1.0,1.0))
    end

    context("FixedVectorNoTuple") do
        # RGB with specified output
        @fact map(-, RGB{Float64}, RGB(1.0, 2.0, 3.0)) --> exactly(RGB{Float64}(-1.0, -2.0, -3.0))
        @fact map(-, RGB{Float64}, [1.0, 2.0, 3.0]) --> exactly(RGB{Float64}(-1.0, -2.0, -3.0))

        # RGB and AbstractVector interop
        @fact map(+, RGB(1.0, 2.0, 3.0), RGB(1.0, 2.0, 3.0)) --> exactly(RGB{Float64}(2.0, 4.0, 6.0))
        @fact map(+, RGB(1.0, 2.0, 3.0), [1.0, 2.0, 3.0]) --> exactly(RGB{Float64}(2.0, 4.0, 6.0))
        @fact map(+, [1.0, 2.0, 3.0], RGB(1.0, 2.0, 3.0)) --> exactly(RGB{Float64}(2.0, 4.0, 6.0))
        @fact map(+, RGB{Int}(1, 2, 3), RGB(1.0, 2.0, 3.0)) --> exactly(RGB{Float64}(2.0, 4.0, 6.0))
    end

    context("Mat and AbstractMatrix") do
        @fact map(+, Mat{2,2,Int}(((1,2),(3,4))), Mat{2,2,Int}(((1,2),(3,4)))) --> exactly(Mat{2,2,Int}(((2,4),(6,8))))
        @fact map(+, Mat{2,2,Int}(((1,2),(3,4))), [1 3; 2 4]) --> exactly(Mat{2,2,Int}(((2,4),(6,8))))
        @fact map(+, [1 3; 2 4], Mat{2,2,Int}(((1,2),(3,4)))) --> exactly(Mat{2,2,Int}(((2,4),(6,8))))
    end

    context("Size checking") do
        @fact_throws DimensionMismatch map(+, Vec(1,2,3), Vec(1,1))
        @fact_throws DimensionMismatch map(+, Vec(1,1), Vec(1,2,3))
        @fact_throws DimensionMismatch map(+, Vec(1,2,3), [1,1])
        @fact_throws DimensionMismatch map(+, [1,1], Vec(1,2,3))
        @fact_throws DimensionMismatch map(+, Vec(1,2,3), 1:2)
        @fact_throws DimensionMismatch map(+, 1:2, Vec(1,2,3))
        @fact_throws DimensionMismatch map(+, Vec(1,2,3), [1 2 3])
    end

    context("Broadcast of scalars") do
        # Arguably not the right thing for map(), but neither do we have a full
        # broadcast implementation yet...
        @fact map(+, Vec{3,Float64}, Vec(1,2,3), 1.0) --> exactly(Vec{3,Float64}(2,3,4))
        @fact map(+, Vec{3,Float64}, 1.0, Vec(1,2,3)) --> exactly(Vec{3,Float64}(2,3,4))
        @fact map(+, 1.0, Vec(1,2,3)) --> exactly(Vec{3,Float64}(2,3,4))
        @fact map(+, Vec(1,2,3), 1.0) --> exactly(Vec{3,Float64}(2,3,4))
    end
end


v2 = Vec(6.0,5.0,4.0)
v1 = Vec(1.0,2.0,3.0)
vi = Vec(1,2,3)
v2 = Vec(6.0,5.0,4.0)
v1c = Vec(6.0+3.0im,5.0-2im,4.0+0.0im)
v2c = v1 + v2*im
v2c = Vec(1.0 + 6.0im, 2.0 + 5.0im, 3.0 + 4.0im)

context("Complex Ops") do
    context("dot product") do
        @fact dot(v1c,v2c) --> dot([6.0+3.0im,5.0-2im,4.0+0.0im], [1.0,2.0,3.0] + [6.0,5.0,4.0]*im)
        @fact Vector(transpose(v1c)*v2c) --> [6.0+3.0im 5.0-2im 4.0+0.0im]*([1.0,2.0,3.0] + [6.0,5.0,4.0]*im)
        @fact Matrix(v2c*transpose(v1c)) --> ([1.0,2.0,3.0] + [6.0,5.0,4.0]*im)*[6.0+3.0im 5.0-2im 4.0+0.0im]
    end
end

context("Destructure") do
    rgb_ref = Int[1 2 3 4;
               2 4 6 8;
               3 6 9 12]
    rgb_ref_set = Int[1 10 10 4;
                   2 10 10 8;
                   3 10 10 12]
    # Test destructure
    rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:4]
    @fact destructure(rgb) --> rgb_ref
    destructure(rgb)[:,2:end-1] = 10
    @fact destructure(rgb) --> rgb_ref_set

    # Explicitly test DestructuredArray.  This wrapper type isn't used by
    # destructure() for plain old dense arrays, since a reinterpret is faster.
    rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:4]
    @fact FixedSizeArrays.DestructuredArray(rgb) --> rgb_ref
    destructure(rgb)[:,2:end-1] = 10
    @fact FixedSizeArrays.DestructuredArray(rgb) --> rgb_ref_set

    # destructure() with 2D FSA
    A = [@fsa([i 2*i; 3*i 4*i]) for i=1:2]
    @fact destructure(A) --> cat(3, [1 2; 3 4], [2 4; 6 8])
end

context("Indexing") do
    context("FixedVector") do
        @fact setindex(v1, 88.9, 1) --> Vec(88.9,2.0,3.0)
        @fact v1[1] --> 1.0
        @fact v1[2] --> 2.0
        @fact v1[3] --> 3.0
        @fact v1[1:3] --> (1.0, 2.0, 3.0)
        @fact v1[1:2] --> (1.0, 2.0)
        @fact v1[1:1] --> (1.0,)
        @fact v1[(1,2)] --> (1.0,2.0)
        @fact v1[(2,1)] --> (2.0,1.0)
        @fact_throws BoundsError v1[-1]
        @fact_throws BoundsError v1[0]
        @fact_throws BoundsError v1[4]
        @fact row(v1, 1) --> (1.0,)
    end
    m = Mat{4,4,Int}(
        (1,2,3,4),
        (5,6,7,8),
        (9,10,11,12),
        (13,14,15,16)
    )
    context("FixedMatrix") do
        @fact setindex(m, 42.0, 2,2) --> Mat{4,4,Int}(
            (1,2,3,4),
            (5,42.0,7,8),
            (9,10,11,12),
            (13,14,15,16)
        )
        @fact m[1] --> 1
        @fact m[2] --> 2
        @fact m[10] --> 10
        @fact m[2,2] --> 6
        @fact m[3,4] --> 15
        @fact m[1:4, 1] --> (1,5,9,13)
        @fact m[1, 1:4] --> (1,2,3,4)
        @fact_throws BoundsError m[-1]
        @fact_throws BoundsError m[0]
        @fact_throws BoundsError m[17]
        @fact_throws BoundsError m[5,1]
        @fact_throws BoundsError m[-1,1]
        @fact_throws BoundsError m[0,0]

        @fact row(m, 1) --> (1,5,9,13)



    end

    context("fslice") do
        context("getindex") do
            rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:10]

            # Plain indexing
            @fact @fslice(rgb[1,2]) --> rgb[2].r
            @fact @fslice(rgb[2,5]) --> rgb[5].g
            @fact @fslice(rgb[3,8]) --> rgb[8].b

            # Slicing along fixed dims
            @fact @fslice(rgb[:,1]) --> rgb[1]
            @fact @fslice(rgb[:,end]) --> rgb[end]

            # Slicing across fixed dims
            @fact compatsqueeze(@fslice(rgb[1,:]))  --> [c.r for c in rgb]
            @fact compatsqueeze(@fslice(rgb[2,:]))  --> [c.g for c in rgb]
            @fact compatsqueeze(@fslice(rgb[3,:]))  --> [c.b for c in rgb]
            # Slicing across fixed dims with field names
            @fact compatsqueeze(@fslice(rgb[:r,:])) --> [c.r for c in rgb]
            @fact compatsqueeze(@fslice(rgb[:g,:])) --> [c.g for c in rgb]
            @fact compatsqueeze(@fslice(rgb[:b,:])) --> [c.b for c in rgb]

            # Slicing FSAs with two fixed dimensions
            N = 3
            A = Mat{2,2,Int}[@fsa([i 2*i; 3*j 4*j]) for i=1:N, j=1:N]
            for i=1:N,j=1:N
                @fact compatsqueeze(@fslice(A[:,:,i,j])) --> A[i,j]
            end
            @fact compatsqueeze(@fslice(A[1,1,:,1])) --> [A[i,1][1,1] for i=1:N]
            @fact compatsqueeze(@fslice(A[end,end,end,:])) --> [A[end,j][end,end] for j=1:N]
            @fact compatsqueeze(@fslice(A[1,[1,end],1,1])) --> [A[1,1][1,1], A[1,1][1,end]]
        end

        context("setindex") do
            rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:10]

            @fslice rgb[:r,:] = -1
            @fslice rgb[:g,:] .+= 1
            @fslice rgb[3,:] = -3
            @fact rgb --> RGB{Int}[RGB(-1,2*i+1,-3) for i=1:10]
        end
    end
end



context("Ops") do
    context("revers") do
        @fact reverse(Vec(1,2,3,4,5)) --> Vec(5,4,3,2,1)
    end
    context("Negation") do
        @fact @inferred(-v1) --> Vec(-1.0,-2.0,-3.0)
        @fact isa(-v1, Vec3d) --> true
    end

    context("Addition") do
        @fact @inferred(v1+v2) --> Vec3d(7.0,7.0,7.0)
        @fact @inferred(RGB(1,2,3) + RGB(2,2,2)) --> exactly(RGB{Int}(3,4,5))
        @fact @inferred(Coord2D(1,2) + Coord2D(3,4)) --> exactly(Coord2D(4,6))
    end
    context("Subtraction") do
        @fact @inferred(v2-v1) --> Vec3d(5.0,3.0,1.0)
        @fact @inferred(RGB(1,2,3) - RGB(2,2,2)) --> exactly(RGB{Int}(-1,0,1))
        @fact @inferred(Coord2D(1,2) - Coord2D(3,4)) --> exactly(Coord2D(-2,-2))
    end
    context("Multiplication") do
        @fact @inferred(v1.*v2) --> Vec3d(6.0,10.0,12.0)
    end
    context("Mixed Type Multiplication") do
        @fact @inferred(vi.*v2) --> Vec3d(6.0,10.0,12.0)
    end
    context("Division") do
        @fact @inferred(v1 ./ v1) --> Vec3d(1.0,1.0,1.0)
    end

    context("Relational") do
        @fact Vec(1,3) .< Vec(2,2) --> exactly(Vec{2,Bool}(true,false))
        @fact RGB(1,2,3) .< RGB(2,2,2) --> exactly(RGB{Bool}(true,false,false))
        @fact Coord2D(1,3) .< Coord2D(2,2) --> exactly(Vec{2,Bool}(true,false))
        end

    context("Scalar") do
        @fact @inferred(1.0 + v1) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(1.0 .+ v1) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(v1 + 1.0) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(v1 .+ 1.0) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(1 + v1) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(1 .+ v1) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(v1 + 1) --> Vec3d(2.0,3.0,4.0)
        @fact @inferred(v1 .+ 1) --> Vec3d(2.0,3.0,4.0)

        @fact @inferred(v1 - 1.0) --> Vec3d(0.0,1.0,2.0)
        @fact @inferred(v1 .- 1.0) --> Vec3d(0.0,1.0,2.0)
        @fact @inferred(1.0 - v1) --> Vec3d(0.0,-1.0,-2.0)
        @fact @inferred(1.0 .- v1) --> Vec3d(0.0,-1.0,-2.0)
        @fact @inferred(v1 - 1) --> Vec3d(0.0,1.0,2.0)
        @fact @inferred(v1 .- 1) --> Vec3d(0.0,1.0,2.0)
        @fact @inferred(1 - v1) --> Vec3d(0.0,-1.0,-2.0)
        @fact @inferred(1 .- v1) --> Vec3d(0.0,-1.0,-2.0)

        @fact @inferred(2.0 * v1) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(2.0 .* v1) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(v1 * 2.0) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(v1 .* 2.0) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(2 * v1) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(2 .* v1) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(v1 * 2) --> Vec3d(2.0,4.0,6.0)
        @fact @inferred(v1 .* 2) --> Vec3d(2.0,4.0,6.0)

        @fact @inferred(v1 / 2.0) --> Vec3d(0.5,1.0,1.5)
        @fact @inferred(v1 ./ 2.0) --> Vec3d(0.5,1.0,1.5)
        @fact @inferred(v1 / 2) --> Vec3d(0.5,1.0,1.5)
        @fact @inferred(v1 ./ 2) --> Vec3d(0.5,1.0,1.5)

        @fact @inferred(12.0 ./ v1) --> Vec3d(12.0,6.0,4.0)
        @fact @inferred(12 ./ v1) --> Vec3d(12.0,6.0,4.0)

        @fact @inferred((v1 .^ 2)) --> Vec3d(1.0,4.0,9.0)
        @fact @inferred((v1 .^ 2.0)) --> Vec3d(1.0,4.0,9.0)
        @fact @inferred((2.0 .^ v1)) --> Vec3d(2.0,4.0,8.0)
        @fact @inferred((2 .^ v1)) --> Vec3d(2.0,4.0,8.0)

                a = Vec(3.2f0)
                @fact @inferred(a+0.2) --> Vec1d(3.2f0+0.2)
                @fact @inferred(0.2+a) --> Vec1d(3.2f0+0.2)
                @fact @inferred(a*0.2) --> Vec1d(3.2f0*0.2)
                @fact @inferred(0.2*a) --> Vec1d(3.2f0*0.2)
                @fact @inferred(a+0.2f0) --> Vec{1,Float32}(3.4f0)
                @fact @inferred(0.2f0+a) --> Vec{1,Float32}(3.4f0)
                @fact @inferred(a*0.2f0) --> Vec{1,Float32}(3.2f0*0.2f0)
                @fact @inferred(0.2f0*a) --> Vec{1,Float32}(3.2f0*0.2f0)
    end
    context("vector norm+cross product") do

        @fact norm(Vec3d(1.0,2.0,2.0)) --> 3.0
        @fact norm(Vec3d(1.0,2.0,2.0),2) --> 3.0
        @fact norm(Vec3d(1.0,2.0,2.0),Inf) --> 2.0
        @fact norm(Vec3d(1.0,2.0,2.0),1) --> 5.0

        # cross product
        @fact cross(v1,v2) --> Vec3d(-7.0,14.0,-7.0)
        @fact isa(cross(v1,v2), Vec3d)  --> true

        @fact cross(vi,v2) --> Vec3d(-7.0,14.0,-7.0)
        @fact isa(cross(vi,v2),Vec3d)  --> true

        a,b = Vec2d(0,1), Vec2d(1,0)
        @fact cross(a,b) --> -1.0
        @fact isa(cross(a,b), Float64) --> true
    end

    context("hypot") do
        a = Vec{2,Int}(1,2)
        b = Vec{2,Float64}(1.,2.)
        @fact hypot(a) --> 2.23606797749979
        @fact hypot(b) --> 2.23606797749979
        @fact hypot(a) == hypot(b) --> true
    end
    context("normalize") do
        a = Vec(3,4)
        b = Vec(3.,4.)
        @fact normalize(a) --> Vec(0.6,0.8)
        @fact normalize(b) --> Vec(0.6,0.8)
    end

    context("reduce") do
        a = rand(Vec{7, Float32})
        x = reduce(+, a)
        y = 0f0
        for elem in a
            y += elem
        end
        @fact y --> x

        a = rand(Mat{7, 9, Cuint})
        x2 = reduce(+, a)
        y2 = Cuint(0)
        for elem in a
            y2 += elem
        end
        @fact y2 --> x2
    end
end


context("Promotion") do
    @fact promote_type(Vec{2,Float64}, Int) --> Vec{2,Float64}
end


# type conversion
context("Conversion 2") do
    @fact isa(convert(Vec3f,v1), Vec3f)  --> true

    @fact isa(convert(Vector{Float64}, v1), Vector{Float64})  --> true
    @fact convert(Vector{Float64}, v1) --> [1.0,2.0,3.0]
end

for T in [UInt, Int, Float32, Float64]
    context("Conversion to Vec{N,$T}") do
        X = map(T, (1,2,3,4,5))

        context("single value conversion") do
            x = X[1]
            for N in 1:4
                @fact convert(Vec{N,T}, x) --> Vec{N,T}(repeated(x,N)...)
            end
        end

        context("conversion from vararg, tuple & array") do
            for N in 1:4
                tup = X[1:N]
                arr = [tup...]
                @fact convert(Vec{N,T}, tup...) --> Vec{N,T}(tup...) "Vec{$N,$T} from vararg"
                @fact convert(Vec{N,T}, tup)    --> Vec{N,T}(tup...) "Vec{$N,$T} from tuple"
                @fact convert(Vec{N,T}, arr)    --> Vec{N,T}(tup...) "Vec{$N,$T} from array"
                @fact convert(Vec, tup...) --> Vec{N,T}(tup...) "Vec from vararg"
                @fact convert(Vec, tup)    --> Vec{N,T}(tup...) "Vec from tuple"
                @fact convert(Vec, arr)    --> Vec{N,T}(tup...) "Vec from array"
            end
        end

        context("conversion from too many args should fail") do
            for N in 1:4
                tup = X[1:N+1]
                arr = [tup...]
                @fact_throws convert(Vec{N,T}, tup...)
                @fact_throws convert(Vec{N,T}, tup)
                @fact_throws convert(Vec{N,T}, arr)
            end
        end

        context("conversion from too few args should fail") do
            for N in 3:5
                tup = X[1:N-1]
                arr = [tup...]
                @fact_throws convert(Vec{N,T}, tup...)
                @fact_throws convert(Vec{N,T}, tup)
                @fact_throws convert(Vec{N,T}, arr)
            end
        end
    end
end

# matrix operations

#typealias Mat1d Matrix1x1{Float64}


zeromat = Mat2d((0.0,0.0),(0.0,0.0))




context("Matrix") do
    @fact map(Float64, zeromat) --> zeromat
    @fact length(Mat2d) --> 4
    @fact length(zeromat) --> 4

    @fact size(Mat2d) --> (2,2)
    @fact size(zeromat) --> (2,2)

    @fact zero(Mat2d) --> zeromat

    for i=1:4, j=1:4
        x1 = rand(i,j)
        @fact @inferred(ctranspose(Mat(x1))) --> Mat(x1')
    end


    v = Vec(1.0,2.0,3.0,4.0)
    r = row(v)
    c = column(v)

    #@fact r' --> c
    #@fact c' --> r

    a = c*r

    b = Mat(
        (1.0,2.0,3.0,4.0),
        (2.0,4.0,6.0,8.0),
        (3.0,6.0,9.0,12.0),
        (4.0,8.0,12.0,16.0)
    )

    x = Mat(
        (1,1,1,),
        (2,2,2,),
        (3,3,3,),
    )
    @fact transpose(x) --> Mat(
        (1,2,3),
        (1,2,3),
        (1,2,3),
    )
    @fact transpose(b) --> b

    @fact length(b) --> 16

    @fact a-->b
    mat30 = Mat(((30.0,),))
    @fact r*c --> mat30


    #@fact row(r, 1) --> v
    #@fact column(c,1) --> v
    #@fact row(r+c',1) --> 2*v
    @fact sum(r) --> sum(v)
    @fact prod(c) --> prod(v)

    @fact eye(Mat3d) --> Mat((1.0,0.0,0.0),
                                (0.0,1.0,0.0),
                                (0.0,0.0,1.0))
    #@fact v*eye(Mat4d)*v --> 30.0
    @fact -r --> -1.0*r
    #@fact diag(diagm(v)) --> v

    # type conversion
    #@fact isa(convert(Matrix1x4{Float32},r),Matrix1x4{Float32})
    jm = rand(4,4)
    im = Mat(jm)
    for i=1:4*2
        @fact jm[i] --> im[i]
    end
    #im = Matrix4x4(jm)
    @fact isa(im, Mat4d)  --> true

    jm2 = convert(Array{Float64,2}, im)
    @fact isa(jm2, Array{Float64,2})  --> true
    @fact jm --> jm2

    #Single valued constructor
    @fact Mat4d(0.0) --> zero(Mat4d)

    a = Vec4d(0)
    b = Vec4d(0,0,0,0)
    @fact a --> b

    v = rand(4)
    m = rand(4,4)
    vfs = Vec(v)
    mfs = Mat(m)
    @fact typeof(vfs) --> Vec4d
    @fact typeof(mfs) --> Mat4d

    # issue #65, wrong
    a = Mat((1,2), (3,4))
    @fact Mat(a) --> a
    b = Mat([1,2,3,4])
    @fact b --> Mat((1,2,3,4))
    @fact b --> Mat([1,2,3,4]'')
end
context("Matrix Math") do
    for i=1:5, j=1:5
        v = rand(j)
        m = rand(i,j)
        m2 = rand(i,j)
        mc = rand(i,j) + im*rand(i,j)
        vfs = Vec(v)
        mfs = Mat(m)
        m2fs = Mat(m2)
        mfsc = Mat(mc)

        vi = randperm(j)
        mi = reshape(randperm(i*j), i, j)
        mi2 = reshape(randperm(i*j), i, j)
        vifs = Vec(vi)
        mifs = Mat(mi)
        mi2fs = Mat(mi2)

        context("Matrix{$i, $j} * Vector{$j}") do
            vm = m * v
            @fact isapprox(@inferred(mfs * vfs), vm)  --> true
            @fact isapprox(@inferred(Matrix(mfs) * vfs), vm)  --> true
            @fact isapprox(@inferred(mfs * Vector(vfs)), vm)  --> true
        end
        context("Matrix{$i, $j} * Matrix{$j, $i}") do
            mm = m * m2'
            @fact isapprox(@inferred(mfs * m2fs'), mm)  --> true
            @fact isapprox(@inferred(Matrix(mfs) * m2fs'), mm)  --> true
            @fact isapprox(@inferred(mfs * Matrix(m2fs')), mm)  --> true
        end
        context("Matrix{$i, $j}*(2I)") do
            mm = m*(2)
            @fact isapprox(@inferred(m*(2I)), mm)  --> true
        end

        # test different element types
        context("Matrix{$i, $j, T} * Vector{$j, U}") do
            vmi = mi * v
            @fact isapprox(@inferred(mifs * vfs), vmi)  --> true
            vmi = m * vi
            @fact isapprox(@inferred(mfs * vifs), vmi)  --> true
            # Custom vector types
            @fact @inferred(eye(Mat{3,3,Float64}) * RGB{Int}(1,2,3)) --> exactly(RGB{Float64}(1,2,3))
        end
        context("Matrix{$i, $j, T} * Matrix{$j, $i, U}") do
            mmi = mi * m2'
            @fact isapprox(@inferred(mifs * m2fs'), mmi)  --> true
            mmi = m * mi2'
            @fact isapprox(@inferred(mfs * mi2fs'), mmi)  --> true
        end

        if i == j
            context("(2*I + I*M)\\v") do
                mm = (2*I+I*m) \ v
                @fact isapprox(@inferred((2*I+I*mfs) \ vfs), mm)  --> true
            end
            context("det(M)") do
                mm = det(m)
                fmm = det(mfs)
                @fact isapprox(fmm, mm)  --> true
            end
            context("trace(M)") do
                mm = trace(m)
                fmm = trace(mfs)
                @fact isapprox(fmm, mm)  --> true
            end
            context("inv(M)") do
                mm = inv(m)
                fmm = inv(mfs)
                @fact isapprox(fmm, mm)  --> true
            end
            context("expm(M)") do
                mm = expm(m)
                fmm = expm(mfs)
                @fact isapprox(fmm, mm)  --> true

                mm = expm(mc)
                fmm = expm(mfsc)
                @fact isapprox(fmm, mm)  --> true
            end
            context("lyap(M,M2*M2')") do
                mm = lyap(m, m2*m2')
                fmm = lyap(mfs, m2fs*m2fs')
                @fact isapprox(fmm, mm) --> true
            end
            context("chol(M2*M2')") do
                mm = full(chol(m2*m2'))
                mm2 = full(chol(map(Mat, m2*m2'))) # Matrix of Mat's
                fmm = chol(m2fs*m2fs')
                @fact isapprox(fmm, mm) --> true
                @fact isapprox(mm, map(first, mm2)) --> true

            end

        else
            context("Matrix{$i, $j} * Matrix{$i, $j}") do
                @fact_throws DimensionMismatch mfs * mfs
            end
        end
        context("transpose M") do
            mm = m'
            fmm = mfs'
            @fact isapprox(fmm, mm)  --> true
        end

        context("ctranspose M") do
            mm = mc'
            fmm = mfsc'
            @fact isapprox(fmm, mm)  --> true
        end
    end
    context("expm(M::Mat{3,3,Float64})") do
        # in practice the precision is eps(), if m has not a triple eigenvalue
        for i in 1:30
            m = (rand(0:1,3,3).*randn(3,3) .+ rand(-3:3,3,3)) # some entries are natural numbers to have higher chance of multiple eigenvalues to trigger all branches
            @fact norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 --> true
            m = m + m' # symmetric
            @fact norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 --> true
            m = 1. *rand(-1:1,3,3) # eigenvalues equal with high probability to test worse case
            @fact norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 --> true
            m = m + m'
            @fact norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 --> true
        end
    end
    context("expm(M::Mat{3,3, BigFloat})") do
        @fact norm(Matrix(expm(Mat(big([0.0 0.0 1.0; -1.0 1.0 0.0; -1.0 0.0 2.0])))) - big([-0.0 0.0  1.0; -0.5  1.0  -0.5; -1.0  0.0  2.0])*e,1) <  10eps(big(1.)) --> true
    end

    context("Matrix * FixedVectorNoTuple") do
        rgb = rand(3)
        m = rand(3,3)
        rgbfs = RGB(rgb)
        mfs = Mat(m)
        @fact isapprox(mfs * rgbfs, m * rgb) --> true
    end

    context("Outer product  Vec{N} * Mat{1,M}") do
        v1 = Vec(1,2)
        v2 = Vec(1,2,3)
        @fact v1*v2' --> Vector(v1)*Vector(v2)'
    end

    context("Large matrix multiply") do
        M = rand(9,9)
        v = rand(9)
        @fact isapprox(Mat(M)*Vec(v), M*v) --> true
    end
end





ac = rand(3)
aci = randperm(3)
bc = rand(3)


a = rand(4)
b = rand(4)
c = rand(4,4)

d = cross(ac, bc)
di = cross(aci, bc)
d2 = a+b
f = c*a
g = c*b
h = c*f
i = dot(f, a)
j = dot(a, g)
k = abs(f)
l = abs(-f)

acfs = Vec(ac)
acifs = Vec(aci)
bcfs = Vec(bc)

afs = Vec(a)
bfs = Vec(b)
cfs = Mat(c)

dfs = cross(acfs, bcfs)
difs = cross(acifs, bcfs)
d2fs = afs+bfs
ffs = cfs*afs
gfs = cfs*bfs
hfs = cfs*ffs
ifs = dot(ffs, afs)
jfs = dot(afs, gfs)
kfs = abs(ffs)
lfs = abs(-ffs)

context("Meta") do
    sym, expr = FixedSizeArrays.gen_functor(:+, 2)
    @fact typeof(sym) --> Symbol
    @fact typeof(expr) --> Expr
end

context("Vector Math") do
    context("all") do
        @fact isapprox(acfs, ac)  --> true
        @fact isapprox(bcfs, bc)  --> true

        @fact isapprox(afs, a) --> true
        @fact isapprox(bfs, b) --> true
        @fact isapprox(cfs, c) --> true

        @fact isapprox(dfs, d) --> true
        @fact isapprox(difs, di) --> true
        @fact isapprox(d2fs, d2) --> true
        @fact isapprox(ffs, f) --> true
        @fact isapprox(gfs, g) --> true
        @fact isapprox(hfs, h) --> true
        @fact isapprox(ifs, i) --> true
        @fact isapprox(jfs, j) --> true
        @fact isapprox(kfs, k) --> true
        @fact isapprox(lfs, l) --> true
        @fact isapprox(lfs, lfs) --> true
    end
end

context("Equality") do
    @fact Vec{3, Int}(1) --> Vec{3, Float64}(1)
    @fact Vec{2, Int}(1) --> not(Vec{3, Float64}(1))
    @fact Vec(1,2,3) --> Vec(1.0,2.0,3.0)
    @fact Vec(1,2,3) --> not(Vec(1.0,4.0,3.0))
    @fact Vec(1,2,3) --> [1,2,3]
    @fact Mat((1,2),(3,4)) --> Mat((1,2),(3,4))
    @fact one(Mat{4,1, Float32}) --> one(Vec{4, Float32})
    @fact isapprox(Vec(1.0,0.0), Vec(1.0,1e-14)) --> true
end
#=
#don't have this yet
let
    a = rand(16)
    b = Mat4d(a)
    @fact b --> reshape(a, (4,4))
    @fact reshape(a, (4,4)) --> b
    @fact b --> not(reshape(a, (2,8)))
end
=#

const unaryOps = (
    -, ~, conj, abs,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    sec, csc, cot, asec, acsc, acot,
    sech, csch, coth, asech, acsch, acoth,
    sinc, cosc, cosd, cotd, cscd, secd,
    sind, tand, acosd, acotd, acscd, asecd,
    asind, atand, rad2deg, deg2rad,
    log, log2, log10, log1p, exponent, exp,
    exp2, expm1, cbrt, sqrt, erf,
    erfc, erfcx, erfi, dawson,

    trunc, round, ceil, floor,
    significand, lgamma,
    gamma, lfact, frexp, modf, airy, airyai,
    airyprime, airyaiprime, airybi, airybiprime,
    besselj0, besselj1, bessely0, bessely1,
    eta, zeta, digamma, real, imag
)

# vec-vec and vec-scalar
const binaryOps = (
    .+, .-, .*, ./, .\, /,
    .==, .!=, .<, .<=, .>, .>=, +, -,
    min, max,

    atan2, besselj, bessely, hankelh1, hankelh2,
    besseli, besselk, beta, lbeta
)




context("mapping operators") do
    context("binary: ") do
        test1 = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
        test2 = (Vec(1,0,typemax(Int)), Mat((typemin(Int),77,1), (2,typemax(Int),5), (-2,3,6)), Vec{4, Float32}(-23.2929))
        for op in binaryOps
            for i=1:length(test1)
                v1 = test1[i]
                v2 = test2[i]
                context("$op with $v1 and $v2") do
                    try # really bad tests, but better than nothing...
                        if applicable(op, v1[1], v2[1]) && typeof(op(v1[1], v2[1])) == eltype(v1)
                            r = op(v1, v2)
                            for j=1:length(v1)
                                @fact r[j] --> op(v1[j], v2[j])
                            end
                        end
                    end
                end
            end
        end
    end
    context("unary: ") do
        test = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
        for op in unaryOps
            for t in test
                context("$op with $t") do
                    try
                        if applicable(op, t[1]) && typeof(op(t[1])) == eltype(t)
                            v = op(t)
                            for i=1:length(v)
                                @fact v[i] --> op(t[i])
                            end
                        end
                    end
                end
            end
        end
    end
end

context("typed round/floor/ceil/trunc") do
    v = Vec((0.8,1.2,-0.3))
    @fact floor(Int, v) --> Vec((0,1,-1))
    @fact ceil( Int, v) --> Vec((1,2,0))
    @fact trunc(Int, v) --> Vec((0,1,0))
    @fact round(Int, v) --> Vec((1,1,0))
end


context("shift, push...") do
    v = Vec(1,2,3)
    p = Point(1,2.)
    @fact @inferred(shift(v)       ) --> Vec(2,3)
    @fact @inferred(shift(p)       ) --> Point(2.)

    @fact @inferred(unshift(v, 42) ) --> Vec(42, 1,2,3)
    @fact @inferred(unshift(v, 42.)) --> Vec(42., 1,2,3)
    @fact @inferred(unshift(p, 42.)) --> Point(42, 1, 2.)

    @fact @inferred(push(v, 42)    ) --> Vec(1,2,3, 42)
    @fact @inferred(push(v, 42.)   ) --> Vec(1,2,3, 42.)
    @fact @inferred(push(p, 42.)   ) --> Point(1,2,42.)

    @fact @inferred(pop(v)      ) --> Vec(1,2)
    @fact @inferred(pop(p)      ) --> Point(1.)

    @fact @inferred(deleteat(v,1)  ) --> Vec(2,3)
    @fact @inferred(deleteat(v,2)  ) --> Vec(1,3)
    @fact @inferred(deleteat(v,3)  ) --> Vec(1,2)
    @fact @inferred(deleteat(p,2)  ) --> Point(1.)

    @fact_throws BoundsError deleteat(v,5)
    @fact_throws BoundsError deleteat(v,-9)

    @fact @inferred(insert(v, 1, 42) ) --> Vec(42,1,2,3)
    @fact @inferred(insert(v, 2, 42) ) --> Vec(1,42,2,3)
    @fact @inferred(insert(v, 3, 42) ) --> Vec(1,2,42,3)
    @fact @inferred(insert(v, 4, 42) ) --> Vec(1,2,3,42)
    @fact @inferred(insert(p, 3, 42.)) --> Point(1,2,42.)

    @fact_throws BoundsError insert(v, 5, 42)
    @fact_throws BoundsError insert(v, 0, 42)
end


context("Base.Test") do
    a = rand(2)
    @fact Base.Test.@test_approx_eq(a, Vec(a)) --> nothing
end

facts("show for subtype") do

    Base.show(io::IO, x::TestType) = print(io, "$(Tuple(x))")  # show for new type

    x = TestType(1, 2)
    @fact string(x) --> "(1,2)"
end

end


end

test()

end





FactCheck.exitstatus()
