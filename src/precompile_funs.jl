
macro fsa_precc(expr)
    caller = expr.args[1]
    args = expr.args[2:end]
    targs = map(x -> :(typeof($x)), args)
    :( Base.precompile($caller, tuple($(targs...))) )
end


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

a = 1
a1 = @fsa([a,2,3])
a2 = @fsa([a 2 3])
a3 = @fsa([a;2;3])
a4 = @fsa([a 2;3 4])
a5 = @fsa([a 2 3;4 5 6])
a6 = @fsa([a 2;3 4;5 6])

@fsa_precc Vec(a,2,3)
@fsa_precc Mat((a,),(2,),(3,))
@fsa_precc Mat(((a,2,3),))
@fsa_precc Mat(((a,3),(2,4)))
@fsa_precc Mat(((a,4),(2,5),(3,6)))
@fsa_precc Mat(((a,3,5),(2,4,6)))

N = 100
a = Point{3, Float32}[Point{3, Float32}(0.7132) for i=1:N]
b = RGB{Float32}[RGB{Float32}(52.293) for i=1:N]

c = Point{3, Float64}[Point{3, Float64}(typemin(Float64)), a..., Point{3, Float64}(typemax(Float64))]
d = RGB{Float64}[RGB(typemin(Float64)), b..., RGB(typemax(Float64))]

@fsa_precc sum(a)
@fsa_precc mean(a)
@fsa_precc sum(b)
@fsa_precc mean(b)

@fsa_precc maximum(c)
@fsa_precc minimum(c)

@fsa_precc maximum(d)
@fsa_precc minimum(d)
@fsa_precc a + 1f0
@fsa_precc b + 1f0



for T=[Float32, Float64, Int, Uint, Uint32, Uint8]
    r = rand(T)
    x = RGB{Int}[RGB(1) for i=1:10]
    @fsa_precc RGB{Float32}(["0.222", "9.8822", "29.999"])
    @fsa_precc RGB{Float32}(0.222, 9.8822, 29.999)
    @fsa_precc typeof(map(RGB{Float32}, x))
    @fsa_precc RGB{T}(r)
    @fsa_precc RGB{T}([r,r,r])              
    @fsa_precc RGB(r,r,r)
    @fsa_precc length(RGB{T}([r,r,r]))
    @fsa_precc length(RGB{T})             
    @fsa_precc eltype(RGB{T}([r,r,r]))    
    @fsa_precc eltype(RGB{T})              
    @fsa_precc  typeof(RGB(r,r,r))           
    @fsa_precc  typeof(RGB{T}(1))            
    @fsa_precc typeof(RGB{T}(1,2,3))        
    @fsa_precc ndims(RGB{T}(1,2,3))         

   @fsa_precc  typeof(RGB{T}(1f0))           
   @fsa_precc  typeof(RGB{T}(1f0,2f0,3f0))  
   @fsa_precc  typeof(RGB{T}(1f0, 2, 3.0))   
   @fsa_precc  typeof(RGB(1f0, 2, 3.0))    
   @fsa_precc RGB{Int}(1f0, 2, 3.0)
end

# A little brutal, but hey.... Better redudantant tests, than not enough tests
for N=1:3:10
    for VT=[Point, Vec, Normal], VT2=[Normal, Vec, Point], ET=[Float32, Int, Uint, Float64], ET2=[Float64, Uint, Int, Float32]
        rand_range  = ET(1):ET(10)
        rand_range2 = ET2(1):ET2(10)
        rn = rand(rand_range, N)
        @fsa_precc VT(rn)
        # parse constructor:
        @fsa_precc VT{N, ET}(map(string, rn))
        # multi constructor
        v1 = VT{N, ET}(rn...)

        @fsa_precc typeof(v1)
        @fsa_precc length(v1)
        @fsa_precc eltype(v1)
        @fsa_precc ndims(v1)

        @fsa_precc length(typeof(v1))
        @fsa_precc eltype(typeof(v1))

 
        # from other FSA without parameters
        @fsa_precc VT2(v1)
        v2 = VT2(v1)

        @fsa_precc typeof(v2)
        @fsa_precc length(v2)
        @fsa_precc eltype(v2)

        # from other FSA with parameters
        for i=1:N
            @fsa_precc VT2{i, ET2}(v1)
        end
        # from single
        r  = rand(rand_range)
        r2 = rand(rand_range2)
        @fsa_precc VT{N, ET}(r)
        @fsa_precc VT{N, ET2}(r)
        @fsa_precc VT{N, ET}(r2)
        @fsa_precc VT{N, ET2}(r2)

        for i=1:N
            @fsa_precc getindex(v1, i)
        end
        x = VT{N, ET}[VT{N, ET}(1) for i=1:10]
        x1 = VT2{N, ET}[VT{N, ET}(1) for i=1:10]
        @fsa_precc map(VT2, x)
        @fsa_precc map(VT, x1)
    end
end


@fsa_precc Vec3f(1,1,1)
@fsa_precc Vec3f(1,1f0,1)
@fsa_precc Vec3f(1f0,1,1.0)

@fsa_precc Vec3f(1)
@fsa_precc Vec3f(0)
@fsa_precc Vec3f(1.0) 			
@fsa_precc Vec3f(1.0f0) 			
@fsa_precc Vec3f(1.0f0)
@fsa_precc Vec(1.0, 1.0, 1.0) 
@fsa_precc Vec2d(Vec3d(1)) 	
@fsa_precc Vec(Vec3d(1), 1.0)
@fsa_precc Vec(Vec3d(1), 1) 
@fsa_precc Vec3d(Vec3f(1.0)) 

v2 = Vec(6.0,5.0,4.0)
v1 = Vec(1.0,2.0,3.0)
v2 = Vec(6.0,5.0,4.0)
@fsa_precc Vec(6.0+3.im,5.0-2im,4.0+0.im)
@fsa_precc v2*im
vim2 = v2*im
@fsa_precc v1 + vim2
v2c = Vec(1.0 + 6.0im, 2.0 + 5.0im, 3.0 + 4.0im)

@fsa_precc setindex(v1, 88.9, 1)
@fsa_precc getindex(v1, 1)

@fsa_precc getindex(v1,1:3)

@fsa_precc getindex(v1,(1,2))
m = Mat{4,4,Int}(
(1,2,3,4),
(5,6,7,8),
(9,10,11,12),
(13,14,15,16)
)
@fsa_precc setindex(m, 42.0, 2,2)
@fsa_precc Mat{4,4,Int}(
(1,2,3,4),
(5,42.0,7,8),
(9,10,11,12),
(13,14,15,16)
)
@fsa_precc getindex(m, 1)
m[2] == 2
m[10] == 10
@fsa_precc getindex(m, 2,2)
m[3,4] == 15
@fsa_precc getindex(m, 1:4, 1)
@fsa_precc getindex(m, 1, 1:4)


@fsa_precc -v1

@fsa_precc 2.0 * v1
@fsa_precc 2.0 .* v1
@fsa_precc v1 * 2.0 
@fsa_precc v1 .* 2.0 
@fsa_precc 2 * v1 == Vec3d(2.0,4.0,6.0)
@fsa_precc 2 .* v1 == Vec3d(2.0,4.0,6.0)
@fsa_precc v1 * 2 == Vec3d(2.0,4.0,6.0)
@fsa_precc v1 .* 2 == Vec3d(2.0,4.0,6.0)

@fsa_precc v1 / 2.0
@fsa_precc v1 ./ 2.0 
@fsa_precc v1 / 2
@fsa_precc v1 ./ 2 

@fsa_precc 12.0 ./ v1
@fsa_precc 12 ./ v1 
@fsa_precc (v1 .^ 2)
@fsa_precc (v1 .^ 2.0) 
@fsa_precc (2.0 .^ v1)
@fsa_precc (2 .^ v1)
@fsa_precc norm(Vec3d(1.0,2.0,2.0))

# cross product
@fsa_precc cross(v1,v2) 

a = Vec{2,Int}(1,2)
b = Vec{2,Float64}(1.,2.)
@fsa_precc hypot(a) 
@fsa_precc hypot(b)
@fsa_precc hypot(a) 



# type conversion
@fsa_precc convert(Vec3f,v1)

@fsa_precc convert(Vector{Float64}, v1)

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
    mfsc = mfs + im*mfs
	vm = m * v
	fsvm = mfs * vfs
	if i == j
		fmm = mfs * mfs
		fmm = det(mfs)
		fmm = inv(mfs)
		fmm = expm(mfs)
		fmm = expm(mfsc)
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


const unaryOps2 = (
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
const binaryOps2 = (
    .+, .-,.*, ./, .\, /,
    .==, .!=, .<, .<=, .>, .>=, +, -,
    min, max,

    atan2, besselj, bessely, hankelh1, hankelh2,
    besseli, besselk, beta, lbeta
)


test1 = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
test2 = (Vec(1,0,typemax(Int)), Mat((typemin(Int),77,1), (2,typemax(Int),5), (-2,3,6)), Vec{4, Float32}(-23.2929))
for op in binaryOps2
    for i=1:length(test1)
        v1 = test1[i]
        v2 = test2[i]
        try # really bad tests, but better than nothing...
            if applicable(op, v1[1], v2[1]) && typeof(op(v1[1], v2[1])) == eltype(v1)
                @fsa_precc op(v1, v2)

            end
        end
    end
end

test = (Vec(1,2,typemax(Int)), Mat((typemin(Int),2,5), (2,3,5), (-2,3,6)), Vec{4, Float32}(0.777))
for op in unaryOps2
    for t in test
        try
            if applicable(op, t[1]) && typeof(op(t[1])) == eltype(t)
                @fsa_precc op(t)
            end
        end
    end
end

