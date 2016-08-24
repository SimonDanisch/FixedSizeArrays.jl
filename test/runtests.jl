module FSAtesting

using FixedSizeArrays
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

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
function run_tests()
    @testset "FixedSizeArrays" begin

    include("typeinf.jl")


    @testset "fsa macro" begin
        a = 1
        a1 = @fsa([a,2,3])
        a2 = @fsa([a 2 3])
        a3 = @fsa([a;2;3])
        a4 = @fsa([a 2;3 4])
        a5 = @fsa([a 2 3;4 5 6])
        a6 = @fsa([a 2;3 4;5 6])

        @test a1 == Vec(a,2,3)
        @test a2 == Mat((a,),(2,),(3,))
        @test a3 == Mat(((a,2,3),))
        @test a4 == Mat(((a,3),(2,4)))
        @test a5 == Mat(((a,4),(2,5),(3,6)))
        @test a6 == Mat(((a,3,5),(2,4,6)))
    end

    @testset "core" begin
        @testset "ndims" begin
            @test ndims(D3) == 3
            @test ndims(Mat) == 2
            @test ndims(Vec) == 1
            @test ndims(Vec(1,2,3)) == 1

            @test ndims(D3{3,3,3}) == 3
            @test ndims(Mat{3,3}) == 2
            @test ndims(Vec{3}) == 1

            @test ndims(D3{3,3,3,Int}) == 3
            @test ndims(Mat{3,3,Int}) == 2
            @test ndims(Vec{3,Int}) == 1
        end
        @testset "size_or" begin
            @test size_or(Mat, nothing) == nothing
            @test size_or(Mat{4}, nothing) == nothing
            @test size_or(Mat{4,4}, nothing) == (4,4)
            @test size_or(Mat{4,4, Float32}, nothing) == (4,4)

            @test size_or(Vec, nothing) == nothing
            @test size_or(Vec{4}, nothing) == (4,)
            @test size_or(Vec{4,Float32}, nothing) == (4,)
            @test size_or(FixedArray, nothing) == nothing

        end
        @testset "eltype_or" begin
            @test eltype_or(Mat, nothing) == nothing
            @test eltype_or(Mat{4}, nothing) == nothing
            @test eltype_or(Mat{4,4}, nothing) == nothing
            @test eltype_or(Mat{4,4, Float32}, nothing) == Float32

            @test eltype_or(Vec, nothing) == nothing
            @test eltype_or(Vec{4}, nothing) == nothing
            @test eltype_or(Vec{4,Float32}, nothing) == Float32

            @test eltype_or(FixedArray, nothing) == nothing

        end
        @testset "ndims_or" begin
            @test ndims_or(Mat, nothing) == 2
            @test ndims_or(Mat{4}, nothing) == 2
            @test ndims_or(Mat{4,4}, nothing) == 2
            @test ndims_or(Mat{4,4, Float32}, nothing) == 2

            @test ndims_or(Vec, nothing) == 1
            @test ndims_or(Vec{4}, nothing) == 1
            @test ndims_or(Vec{4, Float64}, nothing) == 1

            @test ndims_or(FixedArray, nothing) == nothing
        end

        @testset "similar_type" begin
            @test similar_type(Vec{3,Int}, Float32) == Vec{3, Float32}
            @test similar_type(Vec{3}, Float32) == Vec{3, Float32}
            @test similar_type(Vec, Float32, (3,)) == Vec{3, Float32}
            @test similar_type(Vec, Float32, (1,2)) == Mat{1,2, Float32}

            @test similar_type(RGB, Float32) == RGB{Float32}
            @test similar_type(RGB{Float32}, Int) == RGB{Int}
            @test similar_type(RGB{Float32}, Int, (3,)) == RGB{Int}
            @test similar_type(RGB{Float32}, Int, (2,2)) == Mat{2,2,Int}

            @test similar_type(Mat{3,3,Int}, Float32) == Mat{3,3,Float32}
            @test similar_type(Mat, Float32, (3,3))   == Mat{3,3,Float32}
            @test similar_type(Mat{2,2,Int}, (3,3))   == Mat{3,3,Int}

            @test similar_type(Coord2D, Float64, (2,)) == Coord2D
            @test similar_type(Coord2D, Int, (2,))     == Vec{2,Int}
            @test similar_type(Coord2D, Float64, (3,)) == Vec{3,Float64}
        end

        @testset "construct_similar" begin
            @test construct_similar(Vec{3,Int}, (1.0f0,2)) === Vec{2,Float32}(1,2)
            @test construct_similar(Vec{2}, (1,2,3)) === Vec{3,Int}(1,2,3)
            @test construct_similar(Vec, (1.0,2)) === Vec{2,Float64}(1,2)

            @test construct_similar(RGB, (1,2,3)) === RGB{Int}(1,2,3)
            @test construct_similar(RGB{Float32}, (1.0,2.0,3.0)) === RGB{Float64}(1.0,2.0,3.0)

            @test construct_similar(Mat{3,3,Int}, ((1.0f0,2),(1.0,2))) === Mat{2,2,Float64}((1,2),(1,2))
            @test construct_similar(Mat, ((1,2),)) === Mat{2,1,Int}(((1,2),))
        end

    end


    @testset "Array of FixedArrays" begin

        N = 100
        a = Point{3, Float32}[Point{3, Float32}(0.7132) for i=1:N]
        b = RGB{Float32}[RGB{Float32}(52.293) for i=1:N]

        c = Point{3, Float64}[Point{3, Float64}(typemin(Float64)), a..., Point{3, Float64}(typemax(Float64))]
        d = RGB{Float64}[RGB(typemin(Float64)), b..., RGB(typemax(Float64))]

        @testset "reduce" begin
            sa = sum(a)
            ma = mean(a)
            sb = sum(b)
            mb = mean(b)
            for i=1:3
                @test sa[i] ≈ Float32(0.7132*N)
                @test ma[i] ≈ Float32(0.7132*N)/ N

                @test sb[i] ≈ Float32(52.293*N)
                @test mb[i] ≈ Float32(52.293*N)/ N
            end

            @test maximum(c) == Point{3, Float32}(typemax(Float64))
            @test minimum(c) == Point{3, Float32}(typemin(Float64))

            @test maximum(d) == RGB(typemax(Float64))
            @test minimum(d) == RGB(typemin(Float64))

            @test extrema(c) == (minimum(c), maximum(c))
        end

        @testset "array ops" begin
            for op in (.+, .-,.*, ./, .\, +, -)
                @test typeof(op(a, 1f0)) == typeof(a)
                @test typeof(op(1f0, a)) == typeof(a)
            end

            af = a + 1f0
            bf = b + 1f0
            aff = a + Point{3, Float32}(1)
            bff = b + RGB{Float32}(1)
            afd = a .+ 1f0
            bfd = b .+ 1f0
            @inferred(b .* 1f0)
            for i=1:N
                @test a[1] + 1f0 == af[i]
                @test b[1] + 1f0 == bf[i]
                @test a[1] + 1f0 == aff[i]
                @test b[1] + 1f0 == bff[i]
                @test a[1] + 1f0 == afd[i]
                @test b[1] + 1f0 == bfd[i]
            end
        end
        @testset "Show" begin
            m1 = rand(Mat4d, 2)
            m2 = rand(RGB{Float32}, 2)
            m3 = rand(Vec3f, 2)
            println(m1)
            println(m2)
            println(m3)
            showcompact(Point(1,2,3))
        end
    end


    @testset "Constructor FixedVectorNoTuple" begin
        for T=[Float32, Float64, Int, UInt, UInt32, UInt8]
            @testset "$T" begin
                r = rand(T)
                x = RGB{Int}[RGB(1) for i=1:10]
                @test RGB{Float32}(["0.222", "9.8822", "29.999"]) == RGB{Float32}(0.222, 9.8822, 29.999)
                @test RGB(["0.222", "9.8822", "29.999"]) == RGB{Float64}(0.222, 9.8822, 29.999)
                @test typeof(map(RGB{Float32}, x))  == Vector{RGB{Float32}}
                @test RGB{T}(r)                     == RGB(r,r,r)
                @test RGB{T}(Vec(r,r,r))            == RGB(r,r,r)
                @test RGB{T}([r,r,r])               == RGB(r,r,r)
                @test length(RGB{T}([r,r,r]))       == 3
                @test length(RGB{T})                == 3
                @test eltype(RGB{T}([r,r,r]))       == T
                @test eltype(RGB{T})                == T
                @test typeof(RGB(r,r,r))            == RGB{T}
                @test typeof(RGB{T}(1))             == RGB{T}
                @test typeof(RGB{T}(1,2,3))         == RGB{T}
                @test ndims(RGB{T}(1,2,3))          == 1

                @test typeof(RGB{T}(1f0))           == RGB{T}
                @test typeof(RGB{T}(1f0,2f0,3f0))   == RGB{T}
                @test typeof(RGB{T}(1f0, 2, 3.0))   == RGB{T}
                @test typeof(RGB(1f0, 2, 3.0))      == RGB{Float64}
                @test typeof(RGB{Int}(1f0, 2, 3.0)) == RGB{Int}
                @test_throws DimensionMismatch RGB((1,2,3), (2,3,4))
            end
        end
    end

    # A little brutal, but hey.... Better redudantant tests, than not enough tests
    @testset "Constructor " begin
        @testset "Rand" begin
            #Win32 seems to fail for rand(Vec4d)
            @test typeof(rand(Vec4d)) == Vec4d
            @test typeof(rand(Mat4d)) == Mat4d

            @test typeof(rand(Mat{4,2, Int})) == Mat{4,2, Int}
            @test typeof(rand(Vec{7, Int})) == Vec{7, Int}
            @test typeof(rand(Vec{7, Int}, 1:7)) == Vec{7, Int}
            @test typeof(rand(Mat4d, -20f0:0.192f0:230f0)) == Mat4d
            @test typeof(rand(Mat{4,21,Float32}, -20f0:0.192f0:230f0)) == Mat{4,21,Float32}

            x = rand(D3{4,4,4, Float32})
            @test typeof(x) == D3{4,4,4, Float32}
            @test eltype(x) == Float32
            @test size(x) == (4,4,4)
            @test typeof(rand(Vec4d, 5,5)) == Matrix{Vec4d}
        end
        @testset "Randn" begin
            @test typeof(randn(Base.Random.GLOBAL_RNG, Vec4d)) == Vec4d
            @test typeof(randn(Vec4d)) == Vec4d
            @test typeof(randn(Mat4d)) == Mat4d
            @test typeof(randn(Mat{4,2, Complex{Float64}})) == Mat{4,2, Complex{Float64}}
            @test typeof(randn(Vec{7, Complex{Float64}})) == Vec{7, Complex{Float64}}
        end
        @testset "Zero" begin
            @test typeof(zero(Vec4d)) == Vec4d
            @test typeof(zero(Mat4d)) == Mat4d

            @test typeof(zero(Mat{4,2, Int})) == Mat{4,2, Int}
            @test typeof(zero(Vec{7, Int})) == Vec{7, Int}
            @test zero(Vec((1,2))) == Vec((0,0))
            @test zero(Vec((1.0,2.0))) == Vec((0.0,0.0))
        end

        @testset "eye" begin
            @test typeof(eye(Mat4d)) == Mat4d
            @test typeof(eye(Mat{4,2, Int})) == Mat{4,2, Int}
        end
        @testset "one" begin
            x = one(Mat{4,2, Int})
            @test typeof(one(Mat4d)) == Mat4d
            @test typeof(x) == Mat{4,2, Int}
            @test all(x-> x==1, x) == true
        end

        @testset "unit" begin
            u4 = unit(Vec4d, 1)
            u7 = unit(Vec{7, Int}, 7)
            @test typeof(u4) == Vec4d
            @test typeof(u7) == Vec{7, Int}
            @test u4[1] == 1.0
            @test u4[2:end] == (0.,0.,0.)

            @test u7[end] == 1
            @test u7[1:end-1] == (0,0,0,0,0,0)
        end
        for N=(1,10)
            @testset "construction, conversion, $N" begin
                for VT=[Point, Vec], VT2=[Normal, Vec], ET=[Float32, Int, UInt], ET2=[Float64, UInt, Float32]
                    rand_range  = ET(1):ET(10)
                    rand_range2 = ET2(1):ET2(10)
                    rn = rand(rand_range, N)
                    v0 = VT(rn)
                    # parse constructor:
                    @test VT{N, ET}(map(string, rn)) == v0
                    # multi constructor
                    v1 = VT{N, ET}(rn...)
                    @test v1 == v0
                    @test typeof(v1) == VT{N, ET}
                    @test length(v1) == N
                    @test eltype(v1) == ET
                    @test ndims(v1) == 1

                    @test length(typeof(v1)) == N
                    @test eltype(typeof(v1)) == ET

                    for i=1:N
                        @test v1[i] == rn[i]
                    end
                    # from other FSA without parameters
                    v2 = VT2(v1)

                    @test typeof(v2) == VT2{N, ET}
                    @test length(v2) == N
                    @test eltype(v2) == ET
                    for i=1:N
                        @test v2[i] == v1[i]
                    end
                    # from other FSA with parameters
                    for i=1:N
                        v3 = VT2{i, ET2}(v1)
                        @test typeof(v3) == VT2{i, ET2}
                        @test length(v3) == i
                        @test eltype(v3) == ET2
                        for i=1:i
                            @test v3[i] == ET2(v2[i])
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
                        @test v1[i] == r
                        @test v2[i] == ET2(r)
                        @test v3[i] == r2
                        @test v4[i] == ET2(r2)
                    end
                    x = VT{N, ET}[VT{N, ET}(1) for i=1:10]
                    x1 = VT2{N, ET}[VT{N, ET}(1) for i=1:10]
                    x2 = map(VT2, x)
                    x3 = map(VT, x2)
                    @test typeof(x)  == Vector{VT{N, ET}}
                    @test typeof(x1) == Vector{VT2{N, ET}}
                    @test typeof(x2) == Vector{VT2{N, ET}}
                    @test typeof(x3) == Vector{VT{N, ET}}
                    @test x3         == x

                    # Construction with only N, issue #56
                    @test VT{N}(ET(1)) == Vec{N, ET}(1)
                    @test VT{N}(ntuple(x->ET(1), N)...) == Vec{N, ET}(1)
                end
            end
        end
    end


    @testset "Constructors" begin
        @testset "FixedVector: unary, from FixedVector" begin
            @test typeof(Vec3f(1,1,1))     == Vec{3, Float32}
            @test typeof(Vec3f(1,1f0,1))   == Vec{3, Float32}
            @test typeof(Vec3f(1f0,1,1.0)) == Vec{3, Float32}
            @test eltype(Vec3f(1f0,1,1.0)) == Float32

            @test typeof(Vec3f(1))      == Vec{3, Float32}
            @test typeof(Vec3f(0))      == Vec{3, Float32}
            @test Vec3f(1.0)             == Vec(1f0,1f0,1f0)
            @test Vec3f(1.0f0)             == Vec(1f0,1f0,1f0)
            @test Vec3f(1.0f0)             == Vec3f(1)
            @test Vec(1.0, 1.0, 1.0)     == Vec3d(1)
            @test Vec2d(Vec3d(1))         == Vec(1.0, 1.0)
            @test Vec(Vec3d(1), 1.0)     == Vec4d(1)
            @test Vec(Vec3d(1), 1)         == Vec4d(1)
            @test Vec3d(Vec3f(1.0))     == Vec3d(1.0)
        end
    end


    @testset "map" begin
        @testset "Vec and AbstractVector" begin
            # Unary, binary & ternary map with specified output type
            @test map(-, Vec{3,Float64}, Vec(1,2,3)) === Vec{3,Float64}(-1,-2,-3)
            @test map(+, Vec{3,Float64}, [1,2,3], Vec(1,2,3)) === Vec{3,Float64}(2,4,6)
            @test map(+, Vec{3,Float64}, [1,2,3], Vec(1,2,3), 1:3) === Vec{3,Float64}(3,6,9)

            # Unary and binary map with deduced output types
            @test map(-, Vec(1,2,3)) === Vec{3,Int}(-1,-2,-3)
            @test map(+, Vec(1,2,3), [1,2,3]) === Vec{3,Int}(2,4,6)
            @test map(+, [1,2,3], Vec(1,2,3)) === Vec{3,Int}(2,4,6)
            @test map(+, Vec(1,2,3), Vec(1,2,3)) === Vec{3,Int}(2,4,6)
            # Some other `AbstractArray`s
            @test map(+, Vec(1,2,3), 1:3) === Vec{3,Int}(2,4,6)
            @test map(+, 1:3, Vec(1,2,3)) === Vec{3,Int}(2,4,6)

            # Binary map with mixed types
            @test map(>, Vec(0.0,2.0), Vec(1,1)) === Vec{2,Bool}(false,true)
            @test map(+, Vec(0.0,0.0), Vec(1,1)) === Vec{2,Float64}(1.0,1.0)
        end

        @testset "FixedVectorNoTuple" begin
            # RGB with specified output
            @test map(-, RGB{Float64}, RGB(1.0, 2.0, 3.0)) === RGB{Float64}(-1.0, -2.0, -3.0)
            @test map(-, RGB{Float64}, [1.0, 2.0, 3.0]) === RGB{Float64}(-1.0, -2.0, -3.0)

            # RGB and AbstractVector interop
            @test map(+, RGB(1.0, 2.0, 3.0), RGB(1.0, 2.0, 3.0)) === RGB{Float64}(2.0, 4.0, 6.0)
            @test map(+, RGB(1.0, 2.0, 3.0), [1.0, 2.0, 3.0]) === RGB{Float64}(2.0, 4.0, 6.0)
            @test map(+, [1.0, 2.0, 3.0], RGB(1.0, 2.0, 3.0)) === RGB{Float64}(2.0, 4.0, 6.0)
            @test map(+, RGB{Int}(1, 2, 3), RGB(1.0, 2.0, 3.0)) === RGB{Float64}(2.0, 4.0, 6.0)
        end

        @testset "Mat and AbstractMatrix" begin
            @test map(+, Mat{2,2,Int}(((1,2),(3,4))), Mat{2,2,Int}(((1,2),(3,4)))) === Mat{2,2,Int}(((2,4),(6,8)))
            @test map(+, Mat{2,2,Int}(((1,2),(3,4))), [1 3; 2 4]) === Mat{2,2,Int}(((2,4),(6,8)))
            @test map(+, [1 3; 2 4], Mat{2,2,Int}(((1,2),(3,4)))) === Mat{2,2,Int}(((2,4),(6,8)))
        end

        @testset "Size checking" begin
            @test_throws DimensionMismatch map(+, Vec(1,2,3), Vec(1,1))
            @test_throws DimensionMismatch map(+, Vec(1,1), Vec(1,2,3))
            @test_throws DimensionMismatch map(+, Vec(1,2,3), [1,1])
            @test_throws DimensionMismatch map(+, [1,1], Vec(1,2,3))
            @test_throws DimensionMismatch map(+, Vec(1,2,3), 1:2)
            @test_throws DimensionMismatch map(+, 1:2, Vec(1,2,3))
            @test_throws DimensionMismatch map(+, Vec(1,2,3), [1 2 3])
        end

        @testset "Broadcast of scalars" begin
            # Arguably not the right thing for map(), but neither do we have a full
            # broadcast implementation yet...
            @test map(+, Vec{3,Float64}, Vec(1,2,3), 1.0) === Vec{3,Float64}(2,3,4)
            @test map(+, Vec{3,Float64}, 1.0, Vec(1,2,3)) === Vec{3,Float64}(2,3,4)
            @test map(+, 1.0, Vec(1,2,3)) === Vec{3,Float64}(2,3,4)
            @test map(+, Vec(1,2,3), 1.0) === Vec{3,Float64}(2,3,4)
        end
    end


    v2 = Vec(6.0,5.0,4.0)
    v1 = Vec(1.0,2.0,3.0)
    vi = Vec(1,2,3)
    v2 = Vec(6.0,5.0,4.0)
    v1c = Vec(6.0+3.0im,5.0-2im,4.0+0.0im)
    v2c = v1 + v2*im
    v2c = Vec(1.0 + 6.0im, 2.0 + 5.0im, 3.0 + 4.0im)

    @testset "Complex Ops" begin
        @testset "dot product" begin
            @test dot(v1c,v2c) == dot([6.0+3.0im,5.0-2im,4.0+0.0im], [1.0,2.0,3.0] + [6.0,5.0,4.0]*im)
            @test Vector(transpose(v1c)*v2c) == [6.0+3.0im 5.0-2im 4.0+0.0im]*([1.0,2.0,3.0] + [6.0,5.0,4.0]*im)
            @test Matrix(v2c*transpose(v1c)) == ([1.0,2.0,3.0] + [6.0,5.0,4.0]*im)*[6.0+3.0im 5.0-2im 4.0+0.0im]
        end
    end

    @testset "Destructure" begin
        rgb_ref = Int[1 2 3 4;
                   2 4 6 8;
                   3 6 9 12]
        rgb_ref_set = Int[1 10 10 4;
                       2 10 10 8;
                       3 10 10 12]
        # Test destructure
        rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:4]
        @test destructure(rgb) == rgb_ref
        destructure(rgb)[:,2:end-1] = 10
        @test destructure(rgb) == rgb_ref_set

        # Explicitly test DestructuredArray.  This wrapper type isn't used by
        # destructure() for plain old dense arrays, since a reinterpret is faster.
        rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:4]
        @test FixedSizeArrays.DestructuredArray(rgb) == rgb_ref
        destructure(rgb)[:,2:end-1] = 10
        @test FixedSizeArrays.DestructuredArray(rgb) == rgb_ref_set

        # destructure() with 2D FSA
        A = [@fsa([i 2*i; 3*i 4*i]) for i=1:2]
        @test destructure(A) == cat(3, [1 2; 3 4], [2 4; 6 8])
    end

    @testset "Indexing" begin
        @testset "FixedVector" begin
            @test setindex(v1, 88.9, 1) == Vec(88.9,2.0,3.0)
            @test v1[1] == 1.0
            @test v1[2] == 2.0
            @test v1[3] == 3.0
            @test v1[1:3] == (1.0, 2.0, 3.0)
            @test v1[1:2] == (1.0, 2.0)
            @test v1[1:1] == (1.0,)
            @test v1[(1,2)] == (1.0,2.0)
            @test v1[(2,1)] == (2.0,1.0)
            @test_throws BoundsError v1[-1]
            @test_throws BoundsError v1[0]
            @test_throws BoundsError v1[4]
            @test row(v1, 1) == (1.0,)
        end
        m = Mat{4,4,Int}(
            (1,2,3,4),
            (5,6,7,8),
            (9,10,11,12),
            (13,14,15,16)
        )
        @testset "FixedMatrix" begin
            @test setindex(m, 42.0, 2,2) == Mat{4,4,Int}(
                (1,2,3,4),
                (5,42.0,7,8),
                (9,10,11,12),
                (13,14,15,16)
            )
            @test m[1] == 1
            @test m[2] == 2
            @test m[10] == 10
            @test m[2,2] == 6
            @test m[3,4] == 15
            @test m[1:4, 1] == (1,5,9,13)
            @test m[1, 1:4] == (1,2,3,4)
            @test_throws BoundsError m[-1]
            @test_throws BoundsError m[0]
            @test_throws BoundsError m[17]
            @test_throws BoundsError m[5,1]
            @test_throws BoundsError m[-1,1]
            @test_throws BoundsError m[0,0]

            @test row(m, 1) == (1,5,9,13)



        end

        @testset "fslice" begin
            @testset "getindex" begin
                rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:10]

                # Plain indexing
                @test @fslice(rgb[1,2]) == rgb[2].r
                @test @fslice(rgb[2,5]) == rgb[5].g
                @test @fslice(rgb[3,8]) == rgb[8].b

                # Slicing along fixed dims
                @test @fslice(rgb[:,1]) == rgb[1]
                @test @fslice(rgb[:,end]) == rgb[end]

                # Slicing across fixed dims
                @test compatsqueeze(@fslice(rgb[1,:]))  == [c.r for c in rgb]
                @test compatsqueeze(@fslice(rgb[2,:]))  == [c.g for c in rgb]
                @test compatsqueeze(@fslice(rgb[3,:]))  == [c.b for c in rgb]
                # Slicing across fixed dims with field names
                @test compatsqueeze(@fslice(rgb[:r,:])) == [c.r for c in rgb]
                @test compatsqueeze(@fslice(rgb[:g,:])) == [c.g for c in rgb]
                @test compatsqueeze(@fslice(rgb[:b,:])) == [c.b for c in rgb]

                # Slicing FSAs with two fixed dimensions
                N = 3
                A = Mat{2,2,Int}[@fsa([i 2*i; 3*j 4*j]) for i=1:N, j=1:N]
                for i=1:N,j=1:N
                    @test compatsqueeze(@fslice(A[:,:,i,j])) == A[i,j]
                end
                @test compatsqueeze(@fslice(A[1,1,:,1])) == [A[i,1][1,1] for i=1:N]
                @test compatsqueeze(@fslice(A[end,end,end,:])) == [A[end,j][end,end] for j=1:N]
                @test compatsqueeze(@fslice(A[1,[1,end],1,1])) == [A[1,1][1,1], A[1,1][1,end]]
            end

            @testset "setindex" begin
                rgb = RGB{Int}[RGB(i,2*i,3*i) for i=1:10]

                @fslice rgb[:r,:] = -1
                @fslice rgb[:g,:] .+= 1
                @fslice rgb[3,:] = -3
                @test rgb == RGB{Int}[RGB(-1,2*i+1,-3) for i=1:10]
            end
        end
    end



    @testset "Ops" begin
        @testset "Negation" begin
            @test @inferred(-v1) == Vec(-1.0,-2.0,-3.0)
            @test isa(-v1, Vec3d) == true
        end

        @testset "Addition" begin
            @test @inferred(v1+v2) == Vec3d(7.0,7.0,7.0)
            @test @inferred(RGB(1,2,3) + RGB(2,2,2)) === RGB{Int}(3,4,5)
            @test @inferred(Coord2D(1,2) + Coord2D(3,4)) === Coord2D(4,6)
        end
        @testset "Subtraction" begin
            @test @inferred(v2-v1) == Vec3d(5.0,3.0,1.0)
            @test @inferred(RGB(1,2,3) - RGB(2,2,2)) === RGB{Int}(-1,0,1)
            @test @inferred(Coord2D(1,2) - Coord2D(3,4)) === Coord2D(-2,-2)
        end
        @testset "Multiplication" begin
            @test @inferred(v1.*v2) == Vec3d(6.0,10.0,12.0)
        end
        @testset "Mixed Type Multiplication" begin
            @test @inferred(vi.*v2) == Vec3d(6.0,10.0,12.0)
        end
        @testset "Division" begin
            @test @inferred(v1 ./ v1) == Vec3d(1.0,1.0,1.0)
        end

        @testset "Relational" begin
            @test Vec(1,3) .< Vec(2,2) === Vec{2,Bool}(true,false)
            @test RGB(1,2,3) .< RGB(2,2,2) === RGB{Bool}(true,false,false)
            @test Coord2D(1,3) .< Coord2D(2,2) === Vec{2,Bool}(true,false)
            end

        @testset "Scalar" begin
            @test @inferred(1.0 + v1) == Vec3d(2.0,3.0,4.0)
            @test @inferred(1.0 .+ v1) == Vec3d(2.0,3.0,4.0)
            @test @inferred(v1 + 1.0) == Vec3d(2.0,3.0,4.0)
            @test @inferred(v1 .+ 1.0) == Vec3d(2.0,3.0,4.0)
            @test @inferred(1 + v1) == Vec3d(2.0,3.0,4.0)
            @test @inferred(1 .+ v1) == Vec3d(2.0,3.0,4.0)
            @test @inferred(v1 + 1) == Vec3d(2.0,3.0,4.0)
            @test @inferred(v1 .+ 1) == Vec3d(2.0,3.0,4.0)

            @test @inferred(v1 - 1.0) == Vec3d(0.0,1.0,2.0)
            @test @inferred(v1 .- 1.0) == Vec3d(0.0,1.0,2.0)
            @test @inferred(1.0 - v1) == Vec3d(0.0,-1.0,-2.0)
            @test @inferred(1.0 .- v1) == Vec3d(0.0,-1.0,-2.0)
            @test @inferred(v1 - 1) == Vec3d(0.0,1.0,2.0)
            @test @inferred(v1 .- 1) == Vec3d(0.0,1.0,2.0)
            @test @inferred(1 - v1) == Vec3d(0.0,-1.0,-2.0)
            @test @inferred(1 .- v1) == Vec3d(0.0,-1.0,-2.0)

            @test @inferred(2.0 * v1) == Vec3d(2.0,4.0,6.0)
            @test @inferred(2.0 .* v1) == Vec3d(2.0,4.0,6.0)
            @test @inferred(v1 * 2.0) == Vec3d(2.0,4.0,6.0)
            @test @inferred(v1 .* 2.0) == Vec3d(2.0,4.0,6.0)
            @test @inferred(2 * v1) == Vec3d(2.0,4.0,6.0)
            @test @inferred(2 .* v1) == Vec3d(2.0,4.0,6.0)
            @test @inferred(v1 * 2) == Vec3d(2.0,4.0,6.0)
            @test @inferred(v1 .* 2) == Vec3d(2.0,4.0,6.0)

            @test @inferred(v1 / 2.0) == Vec3d(0.5,1.0,1.5)
            @test @inferred(v1 ./ 2.0) == Vec3d(0.5,1.0,1.5)
            @test @inferred(v1 / 2) == Vec3d(0.5,1.0,1.5)
            @test @inferred(v1 ./ 2) == Vec3d(0.5,1.0,1.5)

            @test @inferred(12.0 ./ v1) == Vec3d(12.0,6.0,4.0)
            @test @inferred(12 ./ v1) == Vec3d(12.0,6.0,4.0)

            @test @inferred((v1 .^ 2)) == Vec3d(1.0,4.0,9.0)
            @test @inferred((v1 .^ 2.0)) == Vec3d(1.0,4.0,9.0)
            @test @inferred((2.0 .^ v1)) == Vec3d(2.0,4.0,8.0)
            @test @inferred((2 .^ v1)) == Vec3d(2.0,4.0,8.0)

                    a = Vec(3.2f0)
                    @test @inferred(a+0.2) == Vec1d(3.2f0+0.2)
                    @test @inferred(0.2+a) == Vec1d(3.2f0+0.2)
                    @test @inferred(a*0.2) == Vec1d(3.2f0*0.2)
                    @test @inferred(0.2*a) == Vec1d(3.2f0*0.2)
                    @test @inferred(a+0.2f0) == Vec{1,Float32}(3.4f0)
                    @test @inferred(0.2f0+a) == Vec{1,Float32}(3.4f0)
                    @test @inferred(a*0.2f0) == Vec{1,Float32}(3.2f0*0.2f0)
                    @test @inferred(0.2f0*a) == Vec{1,Float32}(3.2f0*0.2f0)
        end
        @testset "vector norm+cross product" begin

            @test norm(Vec3d(1.0,2.0,2.0)) == 3.0

            # cross product
            @test cross(v1,v2) == Vec3d(-7.0,14.0,-7.0)
            @test isa(cross(v1,v2), Vec3d)  == true

            @test cross(vi,v2) == Vec3d(-7.0,14.0,-7.0)
            @test isa(cross(vi,v2),Vec3d)  == true

            a,b = Vec2d(0,1), Vec2d(1,0)
            @test cross(a,b) == -1.0
            @test isa(cross(a,b), Float64) == true
        end

        @testset "hypot" begin
            a = Vec{2,Int}(1,2)
            b = Vec{2,Float64}(1.,2.)
            @test hypot(a) == 2.23606797749979
            @test hypot(b) == 2.23606797749979
            @test hypot(a) == hypot(b) == true
        end
        @testset "normalize" begin
            a = Vec(3,4)
            b = Vec(3.,4.)
            @test normalize(a) == Vec(0.6,0.8)
            @test normalize(b) == Vec(0.6,0.8)
        end

        @testset "reduce" begin
            a = rand(Vec{7, Float32})
            x = reduce(+, a)
            y = 0f0
            for elem in a
                y += elem
            end
            @test y == x

            a = rand(Mat{7, 9, Cuint})
            x2 = reduce(+, a)
            y2 = Cuint(0)
            for elem in a
                y2 += elem
            end
            @test y2 == x2
        end
    end


    @testset "Promotion" begin
        @test promote_type(Vec{2,Float64}, Int) == Vec{2,Float64}
    end


    # type conversion
    @testset "Conversion 2" begin
        @test isa(convert(Vec3f,v1), Vec3f)  == true

        @test isa(convert(Vector{Float64}, v1), Vector{Float64})  == true
        @test convert(Vector{Float64}, v1) == [1.0,2.0,3.0]
    end

    for T in [UInt, Int, Float32, Float64]
        @testset "Conversion to Vec{N,$T}" begin
            X = map(T, (1,2,3,4,5))

            @testset "single value conversion" begin
                x = X[1]
                for N in 1:4
                    @test convert(Vec{N,T}, x) == Vec{N,T}(repeated(x,N)...)
                end
            end

            @testset "conversion from vararg, tuple & array" begin
                for N in 1:4
                    tup = X[1:N]
                    arr = [tup...]
                    @test convert(Vec{N,T}, tup...) == Vec{N,T}(tup...) "Vec{$N,$T} from vararg"
                    @test convert(Vec{N,T}, tup)    == Vec{N,T}(tup...) "Vec{$N,$T} from tuple"
                    @test convert(Vec{N,T}, arr)    == Vec{N,T}(tup...) "Vec{$N,$T} from array"
                    @test convert(Vec, tup...) == Vec{N,T}(tup...) "Vec from vararg"
                    @test convert(Vec, tup)    == Vec{N,T}(tup...) "Vec from tuple"
                    @test convert(Vec, arr)    == Vec{N,T}(tup...) "Vec from array"
                end
            end

            @testset "conversion from too many args should fail" begin
                for N in 1:4
                    tup = X[1:N+1]
                    arr = [tup...]
                    # @fact_throws convert(Vec{N,T}, tup...)
                    # @fact_throws convert(Vec{N,T}, tup)
                    # @fact_throws convert(Vec{N,T}, arr)
                end
            end

            @testset "conversion from too few args should fail" begin
                for N in 3:5
                    tup = X[1:N-1]
                    arr = [tup...]
                    # @fact_throws convert(Vec{N,T}, tup...)
                    # @fact_throws convert(Vec{N,T}, tup)
                    # @fact_throws convert(Vec{N,T}, arr)
                end
            end
        end
    end

    # matrix operations

    #typealias Mat1d Matrix1x1{Float64}


    zeromat = Mat2d((0.0,0.0),(0.0,0.0))




    @testset "Matrix" begin
        @test map(Float64, zeromat) == zeromat
        @test length(Mat2d) == 4
        @test length(zeromat) == 4

        @test size(Mat2d) == (2,2)
        @test size(zeromat) == (2,2)

        @test zero(Mat2d) == zeromat

        for i=1:4, j=1:4
            x1 = rand(i,j)
            @test @inferred(ctranspose(Mat(x1))) == Mat(x1')
        end


        v = Vec(1.0,2.0,3.0,4.0)
        r = row(v)
        c = column(v)

        #@test r' == c
        #@test c' == r

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
        @test transpose(x) == Mat(
            (1,2,3),
            (1,2,3),
            (1,2,3),
        )
        @test transpose(b) == b

        @test length(b) == 16

        @fact a-->b
        mat30 = Mat(((30.0,),))
        @test r*c == mat30


        #@test row(r, 1) == v
        #@test column(c,1) == v
        #@test row(r+c',1) == 2*v
        @test sum(r) == sum(v)
        @test prod(c) == prod(v)

        @test eye(Mat3d) == Mat((1.0,0.0,0.0),
                                    (0.0,1.0,0.0),
                                    (0.0,0.0,1.0))
        #@test v*eye(Mat4d)*v == 30.0
        @test -r == -1.0*r
        #@test diag(diagm(v)) == v

        # type conversion
        #@fact isa(convert(Matrix1x4{Float32},r),Matrix1x4{Float32})
        jm = rand(4,4)
        im = Mat(jm)
        for i=1:4*2
            @test jm[i] == im[i]
        end
        #im = Matrix4x4(jm)
        @test isa(im, Mat4d)  == true

        jm2 = convert(Array{Float64,2}, im)
        @test isa(jm2, Array{Float64,2})  == true
        @test jm == jm2

        #Single valued constructor
        @test Mat4d(0.0) == zero(Mat4d)

        a = Vec4d(0)
        b = Vec4d(0,0,0,0)
        @test a == b

        v = rand(4)
        m = rand(4,4)
        vfs = Vec(v)
        mfs = Mat(m)
        @test typeof(vfs) == Vec4d
        @test typeof(mfs) == Mat4d

        # issue #65, wrong
        a = Mat((1,2), (3,4))
        @test Mat(a) == a
        b = Mat([1,2,3,4])
        @test b == Mat((1,2,3,4))
        @test b == Mat([1,2,3,4]'')
    end
    @testset "Matrix Math" begin
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

            @testset "Matrix{$i, $j} * Vector{$j}" begin
                vm = m * v
                @test isapprox(@inferred(mfs * vfs), vm)  == true
                @test isapprox(@inferred(Matrix(mfs) * vfs), vm)  == true
                @test isapprox(@inferred(mfs * Vector(vfs)), vm)  == true
            end
            @testset "Matrix{$i, $j} * Matrix{$j, $i}" begin
                mm = m * m2'
                @test isapprox(@inferred(mfs * m2fs'), mm)  == true
                @test isapprox(@inferred(Matrix(mfs) * m2fs'), mm)  == true
                @test isapprox(@inferred(mfs * Matrix(m2fs')), mm)  == true
            end
            @testset "Matrix{$i, $j}*(2I)" begin
                mm = m*(2)
                @test isapprox(@inferred(m*(2I)), mm)  == true
            end

            # test different element types
            @testset "Matrix{$i, $j, T} * Vector{$j, U}" begin
                vmi = mi * v
                @test isapprox(@inferred(mifs * vfs), vmi)  == true
                vmi = m * vi
                @test isapprox(@inferred(mfs * vifs), vmi)  == true
                # Custom vector types
                @test @inferred(eye(Mat{3,3,Float64}) * RGB{Int}(1,2,3)) === RGB{Float64}(1,2,3)
            end
            @testset "Matrix{$i, $j, T} * Matrix{$j, $i, U}" begin
                mmi = mi * m2'
                @test isapprox(@inferred(mifs * m2fs'), mmi)  == true
                mmi = m * mi2'
                @test isapprox(@inferred(mfs * mi2fs'), mmi)  == true
            end

            if i == j
                @testset "(2*I + I*M)\\v" begin
                    mm = (2*I+I*m) \ v
                    @test isapprox(@inferred((2*I+I*mfs) \ vfs), mm)  == true
                end
                @testset "det(M)" begin
                    mm = det(m)
                    fmm = det(mfs)
                    @test isapprox(fmm, mm)  == true
                end
                @testset "trace(M)" begin
                    mm = trace(m)
                    fmm = trace(mfs)
                    @test isapprox(fmm, mm)  == true
                end
                @testset "inv(M)" begin
                    mm = inv(m)
                    fmm = inv(mfs)
                    @test isapprox(fmm, mm)  == true
                end
                @testset "expm(M)" begin
                    mm = expm(m)
                    fmm = expm(mfs)
                    @test isapprox(fmm, mm)  == true

                    mm = expm(mc)
                    fmm = expm(mfsc)
                    @test isapprox(fmm, mm)  == true
                end
                @testset "lyap(M,M2*M2')" begin
                    mm = lyap(m, m2*m2')
                    fmm = lyap(mfs, m2fs*m2fs')
                    @test isapprox(fmm, mm) == true
                end
                @testset "chol(M2*M2')" begin
                    mm = full(chol(m2*m2'))
                    mm2 = full(chol(map(Mat, m2*m2'))) # Matrix of Mat's
                    fmm = chol(m2fs*m2fs')
                    @test isapprox(fmm, mm) == true
                    @test isapprox(mm, map(first, mm2)) == true

                end

            else
                @testset "Matrix{$i, $j} * Matrix{$i, $j}" begin
                    @test_throws DimensionMismatch mfs * mfs
                end
            end
            @testset "transpose M" begin
                mm = m'
                fmm = mfs'
                @test isapprox(fmm, mm)  == true
            end

            @testset "ctranspose M" begin
                mm = mc'
                fmm = mfsc'
                @test isapprox(fmm, mm)  == true
            end
        end
        @testset "expm(M::Mat{3,3,Float64})" begin
            # in practice the precision is eps(), if m has not a triple eigenvalue
            for i in 1:30
                m = (rand(0:1,3,3).*randn(3,3) .+ rand(-3:3,3,3)) # some entries are natural numbers to have higher chance of multiple eigenvalues to trigger all branches
                @test norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 == true
                m = m + m' # symmetric
                @test norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 == true
                m = 1. *rand(-1:1,3,3) # eigenvalues equal with high probability to test worse case
                @test norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 == true
                m = m + m'
                @test norm(Matrix(expm(Mat(m))) -  expm(m))/norm(expm(m)) <= 1E-9 == true
            end
        end
        @testset "expm(M::Mat{3,3, BigFloat})" begin
            @test norm(Matrix(expm(Mat(big([0.0 0.0 1.0; -1.0 1.0 0.0; -1.0 0.0 2.0])))) - big([-0.0 0.0  1.0; -0.5  1.0  -0.5; -1.0  0.0  2.0])*e,1) <  10eps(big(1.)) == true
        end

        @testset "Matrix * FixedVectorNoTuple" begin
            rgb = rand(3)
            m = rand(3,3)
            rgbfs = RGB(rgb)
            mfs = Mat(m)
            @test isapprox(mfs * rgbfs, m * rgb) == true
        end

        @testset "Outer product  Vec{N} * Mat{1,M}" begin
            v1 = Vec(1,2)
            v2 = Vec(1,2,3)
            @test v1*v2' == Vector(v1)*Vector(v2)'
        end

        @testset "Large matrix multiply" begin
            M = rand(9,9)
            v = rand(9)
            @test isapprox(Mat(M)*Vec(v), M*v) == true
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

    @testset "Meta" begin
        sym, expr = FixedSizeArrays.gen_functor(:+, 2)
        @test typeof(sym) == Symbol
        @test typeof(expr) == Expr
    end

    @testset "Vector Math" begin
        @testset "all" begin
            @test isapprox(acfs, ac)  == true
            @test isapprox(bcfs, bc)  == true

            @test isapprox(afs, a) == true
            @test isapprox(bfs, b) == true
            @test isapprox(cfs, c) == true

            @test isapprox(dfs, d) == true
            @test isapprox(difs, di) == true
            @test isapprox(d2fs, d2) == true
            @test isapprox(ffs, f) == true
            @test isapprox(gfs, g) == true
            @test isapprox(hfs, h) == true
            @test isapprox(ifs, i) == true
            @test isapprox(jfs, j) == true
            @test isapprox(kfs, k) == true
            @test isapprox(lfs, l) == true
            @test isapprox(lfs, lfs) == true
        end
    end

    @testset "Equality" begin
        @test Vec{3, Int}(1) == Vec{3, Float64}(1)
        @test Vec{2, Int}(1) != Vec{3, Float64}(1)
        @test Vec(1,2,3) == Vec(1.0,2.0,3.0)
        @test Vec(1,2,3) != Vec(1.0,4.0,3.0)
        @test Vec(1,2,3) == [1,2,3]
        @test Mat((1,2),(3,4)) == Mat((1,2),(3,4))
        @test one(Mat{4,1, Float32}) == one(Vec{4, Float32})
        @test isapprox(Vec(1.0,0.0), Vec(1.0,1e-14)) == true
    end
    #=
    #don't have this yet
    let
        a = rand(16)
        b = Mat4d(a)
        @test b == reshape(a, (4,4))
        @test reshape(a, (4,4)) == b
        @test b != reshape(a, (2,8))
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




    @testset "mapping operators" begin
        @testset "binary: " begin
            test1 = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
            test2 = (Vec(1,0,typemax(Int)), Mat((typemin(Int),77,1), (2,typemax(Int),5), (-2,3,6)), Vec{4, Float32}(-23.2929))
            for op in binaryOps
                for i=1:length(test1)
                    v1 = test1[i]
                    v2 = test2[i]
                    @testset "$op with $v1 and $v2" begin
                        try # really bad tests, but better than nothing...
                            if applicable(op, v1[1], v2[1]) && typeof(op(v1[1], v2[1])) == eltype(v1)
                                r = op(v1, v2)
                                for j=1:length(v1)
                                    @test r[j] == op(v1[j], v2[j])
                                end
                            end
                        end
                    end
                end
            end
        end
        @testset "unary: " begin
            test = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
            for op in unaryOps
                for t in test
                    @testset "$op with $t" begin
                        try
                            if applicable(op, t[1]) && typeof(op(t[1])) == eltype(t)
                                v = op(t)
                                for i=1:length(v)
                                    @test v[i] == op(t[i])
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "typed round/floor/ceil/trunc" begin
        v = Vec((0.8,1.2,-0.3))
        @test floor(Int, v) == Vec((0,1,-1))
        @test ceil( Int, v) == Vec((1,2,0))
        @test trunc(Int, v) == Vec((0,1,0))
        @test round(Int, v) == Vec((1,1,0))
    end


    @testset "shift, push..." begin
        v = Vec(1,2,3)
        p = Point(1,2.)
        @test @inferred(shift(v)       ) == Vec(2,3)
        @test @inferred(shift(p)       ) == Point(2.)

        @test @inferred(unshift(v, 42) ) == Vec(42, 1,2,3)
        @test @inferred(unshift(v, 42.)) == Vec(42., 1,2,3)
        @test @inferred(unshift(p, 42.)) == Point(42, 1, 2.)

        @test @inferred(push(v, 42)    ) == Vec(1,2,3, 42)
        @test @inferred(push(v, 42.)   ) == Vec(1,2,3, 42.)
        @test @inferred(push(p, 42.)   ) == Point(1,2,42.)

        @test @inferred(pop(v)      ) == Vec(1,2)
        @test @inferred(pop(p)      ) == Point(1.)

        @test @inferred(deleteat(v,1)  ) == Vec(2,3)
        @test @inferred(deleteat(v,2)  ) == Vec(1,3)
        @test @inferred(deleteat(v,3)  ) == Vec(1,2)
        @test @inferred(deleteat(p,2)  ) == Point(1.)

        @test_throws BoundsError deleteat(v,5)
        @test_throws BoundsError deleteat(v,-9)

        @test @inferred(insert(v, 1, 42) ) == Vec(42,1,2,3)
        @test @inferred(insert(v, 2, 42) ) == Vec(1,42,2,3)
        @test @inferred(insert(v, 3, 42) ) == Vec(1,2,42,3)
        @test @inferred(insert(v, 4, 42) ) == Vec(1,2,3,42)
        @test @inferred(insert(p, 3, 42.)) == Point(1,2,42.)

        @test_throws BoundsError insert(v, 5, 42)
        @test_throws BoundsError insert(v, 0, 42)
    end


    @testset "Base.Test" begin
        a = rand(2)
        @test Base.Test.@test_approx_eq(a, Vec(a)) == nothing
    end

    @testset "show for subtype" begin

        Base.show(io::IO, x::TestType) = print(io, "$(Tuple(x))")  # show for new type

        x = TestType(1, 2)
        @test string(x) == "(1,2)"
    end

    end

end

run_tests()

end
