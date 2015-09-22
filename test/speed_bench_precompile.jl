using FixedSizeArrays
using FactCheck


immutable Normal{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
end
immutable RGB{T} <: FixedVectorNoTuple{3, T}
    r::T
    g::T
    b::T
end

typealias Vec1d Vec{1, Float64}
typealias Vec2d Vec{2, Float64}
typealias Vec3d Vec{3, Float64}
typealias Vec4d Vec{4, Float64}
typealias Vec3f Vec{3, Float32}
typealias Mat2d Mat{2,2, Float64}
typealias Mat3d Mat{3,3, Float64}
typealias Mat4d Mat{4,4, Float64}


function test()
a = 1
a1 = @fsa([a,2,3])
a2 = @fsa([a 2 3])
a3 = @fsa([a;2;3])
a4 = @fsa([a 2;3 4])
a5 = @fsa([a 2 3;4 5 6])
a6 = @fsa([a 2;3 4;5 6])

a1 == Vec(a,2,3)
a2 == Mat((a,),(2,),(3,))
a3 == Mat(((a,2,3),))
a4 == Mat(((a,3),(2,4)))
a5 == Mat(((a,4),(2,5),(3,6)))
a6 == Mat(((a,3,5),(2,4,6)))

N = 100
a = Point{3, Float32}[Point{3, Float32}(0.7132) for i=1:N]
b = RGB{Float32}[RGB{Float32}(52.293) for i=1:N]

c = Point{3, Float64}[Point{3, Float64}(typemin(Float64)), a..., Point{3, Float64}(typemax(Float64))]
d = RGB{Float64}[RGB(typemin(Float64)), b..., RGB(typemax(Float64))]

sa = sum(a)
ma = mean(a)
sb = sum(b)
mb = mean(b)

maximum(c) == Point{3, Float32}(typemax(Float64))
minimum(c) == Point{3, Float32}(typemin(Float64))

maximum(d) == RGB(typemax(Float64))
minimum(d) == RGB(typemin(Float64))
af = a + 1f0
bf = b + 1f0
for elem in af
    a[1] + 1f0 == elem
end
for elem in bf
    b[1] + 1f0 == elem
end


for T=[Float32, Float64, Int, UInt, UInt32, UInt8]
    r = rand(T)
    x = RGB{Int}[RGB(1) for i=1:10]
    RGB{Float32}(["0.222", "9.8822", "29.999"]) == RGB{Float32}(0.222, 9.8822, 29.999)
    typeof(map(RGB{Float32}, x))  == Vector{RGB{Float32}}
    RGB{T}(r)                     == RGB(r,r,r)
    RGB{T}([r,r,r])               == RGB(r,r,r)
    length(RGB{T}([r,r,r]))       == 3
    length(RGB{T})                == 3
    eltype(RGB{T}([r,r,r]))       == T
    eltype(RGB{T})                == T
    typeof(RGB(r,r,r))            == RGB{T}
    typeof(RGB{T}(1))             == RGB{T}
    typeof(RGB{T}(1,2,3))         == RGB{T}
    ndims(RGB{T}(1,2,3))          == 1

    typeof(RGB{T}(1f0))           == RGB{T}
    typeof(RGB{T}(1f0,2f0,3f0))   == RGB{T}
    typeof(RGB{T}(1f0, 2, 3.0))   == RGB{T}
    typeof(RGB(1f0, 2, 3.0))      == RGB{Float64}
    typeof(RGB{Int}(1f0, 2, 3.0)) == RGB{Int}
end

# A little brutal, but hey.... Better redudantant tests, than not enough tests
for N=1:3:10
    for VT=[Point, Vec, Normal], VT2=[Normal, Vec, Point], ET=[Float32, Int, UInt64, Float64], ET2=[Float64, UInt64, Int, Float32]
        rand_range  = ET(1):ET(10)
        rand_range2 = ET2(1):ET2(10)
        rn = rand(rand_range, N)
        v0 = VT(rn)
        # parse constructor:
        VT{N, ET}(map(string, rn)) == v0
        # multi constructor
        v1 = VT{N, ET}(rn...)
        v1 == v0
        typeof(v1) == VT{N, ET}
        length(v1) == N
        eltype(v1) == ET
        ndims(v1) == 1

        length(typeof(v1)) == N
        eltype(typeof(v1)) == ET

        for i=1:N
            v1[i] == rn[i]
        end
        # from other FSA without parameters
        v2 = VT2(v1)

        typeof(v2) == VT2{N, ET}
        length(v2) == N
        eltype(v2) == ET
        for i=1:N
            v2[i] == v1[i]
        end
        # from other FSA with parameters
        for i=1:N
            v3 = VT2{i, ET2}(v1)
            typeof(v3) == VT2{i, ET2}
            length(v3) == i
            eltype(v3) == ET2
            for i=1:i
                v3[i] == ET2(v2[i])
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
            v1[i] == r
            v2[i] == ET2(r)
            v3[i] == r2
            v4[i] == ET2(r2)
        end
        x = VT{N, ET}[VT{N, ET}(1) for i=1:10]
        x1 = VT2{N, ET}[VT{N, ET}(1) for i=1:10]
        x2 = map(VT2, x)
        x3 = map(VT, x2)
        typeof(x)  == Vector{VT{N, ET}}
        typeof(x1) == Vector{VT2{N, ET}}
        typeof(x2) == Vector{VT2{N, ET}}
        typeof(x3) == Vector{VT{N, ET}}
        x3         == x
    end
end


typeof(Vec3f(1,1,1))     == Vec{3, Float32}
typeof(Vec3f(1,1f0,1))   == Vec{3, Float32}
typeof(Vec3f(1f0,1,1.0)) == Vec{3, Float32}

typeof(Vec3f(1))  	== Vec{3, Float32}
typeof(Vec3f(0))  	== Vec{3, Float32}
Vec3f(1.0) 			== Vec(1f0,1f0,1f0)
Vec3f(1.0f0) 			== Vec(1f0,1f0,1f0)
Vec3f(1.0f0) 			== Vec3f(1)
Vec(1.0, 1.0, 1.0) 	== Vec3d(1)
Vec2d(Vec3d(1)) 		== Vec(1.0, 1.0)
Vec(Vec3d(1), 1.0) 	== Vec4d(1)
Vec(Vec3d(1), 1) 		== Vec4d(1)
Vec3d(Vec3f(1.0)) 	== Vec3d(1.0)

v2 = Vec(6.0,5.0,4.0)
v1 = Vec(1.0,2.0,3.0)
v2 = Vec(6.0,5.0,4.0)

setindex(v1, 88.9, 1) == Vec(88.9,2.0,3.0)
v1[1] == 1.0
v1[2] == 2.0
v1[3] == 3.0
v1[1:3] == (1.0, 2.0, 3.0)
v1[1:2] == (1.0, 2.0)
v1[1:1] == (1.0,)
v1[(1,2)] == (1.0,2.0)
v1[(2,1)] == (2.0,1.0)
m = Mat{4,4,Int}(
(1,2,3,4),
(5,6,7,8),
(9,10,11,12),
(13,14,15,16)
)
setindex(m, 42.0, 2,2) == Mat{4,4,Int}(
(1,2,3,4),
(5,42.0,7,8),
(9,10,11,12),
(13,14,15,16)
)
m[1] == 1
m[2] == 2
m[10] == 10
m[2,2] == 6
m[3,4] == 15
m[1:4, 1] == (1,5,9,13)
m[1, 1:4] == (1,2,3,4)


-v1 == Vec(-1.0,-2.0,-3.0)
isa(-v1, Vec3d) == true
v1+v2 == Vec3d(7.0,7.0,7.0)
v2-v1 == Vec3d(5.0,3.0,1.0)
v1.*v2 == Vec3d(6.0,10.0,12.0)
v1 ./ v1 == Vec3d(1.0,1.0,1.0)
1.0 + v1 == Vec3d(2.0,3.0,4.0)
1.0 .+ v1 == Vec3d(2.0,3.0,4.0)
v1 + 1.0 == Vec3d(2.0,3.0,4.0)
v1 .+ 1.0 == Vec3d(2.0,3.0,4.0)
1 + v1 == Vec3d(2.0,3.0,4.0)
1 .+ v1 == Vec3d(2.0,3.0,4.0)
v1 + 1 == Vec3d(2.0,3.0,4.0)
v1 .+ 1 == Vec3d(2.0,3.0,4.0)

v1 - 1.0 == Vec3d(0.0,1.0,2.0)
v1 .- 1.0 == Vec3d(0.0,1.0,2.0)
1.0 - v1 == Vec3d(0.0,-1.0,-2.0)
1.0 .- v1 == Vec3d(0.0,-1.0,-2.0)
v1 - 1 == Vec3d(0.0,1.0,2.0)
v1 .- 1 == Vec3d(0.0,1.0,2.0)
1 - v1 == Vec3d(0.0,-1.0,-2.0)
1 .- v1 == Vec3d(0.0,-1.0,-2.0)

2.0 * v1 == Vec3d(2.0,4.0,6.0)
2.0 .* v1 == Vec3d(2.0,4.0,6.0)
v1 * 2.0 == Vec3d(2.0,4.0,6.0)
v1 .* 2.0 == Vec3d(2.0,4.0,6.0)
2 * v1 == Vec3d(2.0,4.0,6.0)
2 .* v1 == Vec3d(2.0,4.0,6.0)
v1 * 2 == Vec3d(2.0,4.0,6.0)
v1 .* 2 == Vec3d(2.0,4.0,6.0)

v1 / 2.0 == Vec3d(0.5,1.0,1.5)
v1 ./ 2.0 == Vec3d(0.5,1.0,1.5)
v1 / 2 == Vec3d(0.5,1.0,1.5)
v1 ./ 2 == Vec3d(0.5,1.0,1.5)

12.0 ./ v1 == Vec3d(12.0,6.0,4.0)
12 ./ v1 == Vec3d(12.0,6.0,4.0)

(v1 .^ 2) == Vec3d(1.0,4.0,9.0)
(v1 .^ 2.0) == Vec3d(1.0,4.0,9.0)
(2.0 .^ v1) == Vec3d(2.0,4.0,8.0)
(2 .^ v1) == Vec3d(2.0,4.0,8.0)
norm(Vec3d(1.0,2.0,2.0)) == 3.0

# cross product
cross(v1,v2) == Vec3d(-7.0,14.0,-7.0)
isa(cross(v1,v2),Vec3d)  == true
a = Vec{2,Int}(1,2)
b = Vec{2,Float64}(1.,2.)
hypot(a) == 2.23606797749979
hypot(b) == 2.23606797749979
hypot(a) == hypot(b) == true




# type conversion
isa(convert(Vec3f,v1), Vec3f)  == true

isa(convert(Vector{Float64}, v1), Vector{Float64})  == true
convert(Vector{Float64}, v1) == [1.0,2.0,3.0]

# matrix operations

#typealias Mat1d Matrix1x1{Float64}

zeromat = Mat2d((0.0,0.0),(0.0,0.0))


length(Mat2d) == 4
length(zeromat) == 4

size(Mat2d) == (2,2)
size(zeromat) == (2,2)

zero(Mat2d) == zeromat

for i=1:4, j=1:4
    x1 = rand(i,j)
    Mat(x1') == Mat(x1)'
end


v = Vec(1.0,2.0,3.0,4.0)
r = row(v)
c = column(v)

#r' == c
#c' == r

a = c*r

b = Mat(
	(1.0,2.0,3.0,4.0),
	(2.0,4.0,6.0,8.0),
	(3.0,6.0,9.0,12.0),
	(4.0,8.0,12.0,16.0)
)

length(b) == 16

a==b
mat30 = Mat(((30.0,),))
r*c == mat30

#row(r, 1) == v
#column(c,1) == v
#row(r+c',1) == 2*v
sum(r) == sum(v)
prod(c) == prod(v)

eye(Mat3d) == Mat((1.0,0.0,0.0),
	(0.0,1.0,0.0),
	(0.0,0.0,1.0))
#v*eye(Mat4d)*v == 30.0
x = -r
y = -1.0*r
jm = rand(4,4)
lowl = Mat(jm)
jm2 = convert(Array{Float64,2}, lowl)

#Single valued constructor
Mat4d(0.0) 
a = Vec4d(0)
b = Vec4d(0,0,0,0)

v = rand(4)
m = rand(4,4)
vfs = Vec(v)
mfs = Mat(m)

for i=1:4, j=1:4
	vfs = rand(Vec{j, Float64})
	mfs = rand(Mat{i,j, Float64})
	vm = m * v
	fsvm = mfs * vfs
	if i == j
		fmm = mfs * mfs
		fmm = det(mfs)
		fmm = inv(mfs)
		fmm = expm(mfs)
    end
end
mm = m'
fmm = mfs'


ac = rand(3)
bc = rand(3)

a = rand(4)
b = rand(4)
c = rand(4,4)

d = cross(ac, bc)
d2 = a+b
f = c*a
g = c*b
h = c*f
i = dot(f, a)
j = dot(a, g)
k = abs(f)
l = abs(-f)

acfs = Vec(ac)
bcfs = Vec(bc)

afs = Vec(a)
bfs = Vec(b)
cfs = Mat(c)

dfs = cross(acfs, bcfs)
d2fs = afs+bfs
ffs = cfs*afs
gfs = cfs*bfs
hfs = cfs*ffs
ifs = dot(ffs, afs)
jfs = dot(afs, gfs)
kfs = abs(ffs)
lfs = abs(-ffs)

isapprox(acfs, ac)  == true
isapprox(bcfs, bc)  == true

isapprox(afs, a) == true
isapprox(bfs, b) == true
isapprox(cfs, c) == true

isapprox(dfs, d) == true
isapprox(d2fs, d2) == true
isapprox(ffs, f) == true
isapprox(gfs, g) == true
isapprox(hfs, h) == true
isapprox(ifs, i) == true
isapprox(jfs, j) == true
isapprox(kfs, k) == true
isapprox(lfs, l) == true

Vec{3, Int}(1) == Vec{3, Float64}(1)
Vec(1,2,3) == Vec(1.0,2.0,3.0)
Vec(1,2,3) == [1,2,3]
Mat((1,2),(3,4)) == Mat((1,2),(3,4))


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

    #trunc, round, ceil, floor, #see JuliaLang/julia#12163
    significand, lgamma,
    gamma, lfact, frexp, modf, airy, airyai,
    airyprime, airyaiprime, airybi, airybiprime,
    besselj0, besselj1, bessely0, bessely1,
    eta, zeta, digamma
)

# vec-vec and vec-scalar
const binaryOps = (
    .+, .-,.*, ./, .\, /,
    .==, .!=, .<, .<=, .>, .>=, +, -,
    min, max,

    atan2, besselj, bessely, hankelh1, hankelh2,
    besseli, besselk, beta, lbeta
)


test1 = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
test2 = (Vec(1,0,typemax(Int)), Mat((typemin(Int),77,1), (2,typemax(Int),5), (-2,3,6)), Vec{4, Float32}(-23.2929))
for op in binaryOps
    for i=1:length(test1)
        v1 = test1[i]
        v2 = test2[i]
        try # really bad tests, but better than nothing...
            if applicable(op, v1[1], v2[1]) && typeof(op(v1[1], v2[1])) == eltype(v1)
                r = op(v1, v2)
                for j=1:length(v1)
                    r[j] == op(v1[j], v2[j])
                end

            end
        end
    end
end

test = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
for op in unaryOps
    for t in test
        try
            if applicable(op, t[1]) && typeof(op(t[1])) == eltype(t)
                v = op(t)
                for i=1:length(v)
                    v[i] == op(t[i])
                end
            end
        end
    end
end


end

using SnoopCompile

snoop_on()
SnoopCompile.@snoop "/tmp/fsa_compiles.csv" begin
   test()
end
snoop_off()

data = SnoopCompile.read("/tmp/fsa_compiles.csv")
pc, discards = SnoopCompile.parcel(data[end:-1:1,2])
SnoopCompile.write(pc, "/tmp/precompile")

