using FixedSizeArrays
using FactCheck

immutable Vec{N, T} <: FixedVector{N, T}
    _::NTuple{N, T}
end

immutable RGB{T} <: FixedVectorNoTuple{3, T}
    r::T
    g::T
    b::T
end
sleep(0.1)

facts("FixedVectorNoTuple") do
	context("Constructor") do 
		@fact typeof(RGB(1,2,3)) => RGB{Int}
		@fact typeof(RGB(1f0,2f0,3f0)) => RGB{Float32}
		@fact typeof(RGB{Float32}(1,2,3)) => RGB{Float32}
	end
end
typealias Vec2d Vec{2, Float64}
typealias Vec3d Vec{3, Float64}
typealias Vec4d Vec{4, Float64}
typealias Vec3f Vec{3, Float32}
Vec(Vec3d(1), 1.0)
#t1 = ["1.909", "1.909", "1.909"]
#@fact Vec{3, Float64}(1.909) => Vec{3, Float64}(t1)
#@fact length(t1) => 3

facts("Constructors") do
	context("FixedVector: unary, from FixedVector") do 
		@fact typeof(Vec3f(1,1,1))     => Vec{3, Float32}
		@fact typeof(Vec3f(1,1f0,1))   => Vec{3, Float32}
		@fact typeof(Vec3f(1f0,1,1.0)) => Vec{3, Float32}

		@fact typeof(Vec3f(1))  	=> Vec{3, Float32}
		@fact typeof(Vec3f(0))  	=> Vec{3, Float32}
		@fact Vec3f(1.0f0) 			=> Vec(1f0,1f0,1f0)
		@fact Vec3f(1.0f0) 			=> Vec(1f0,1f0,1f0)
		@fact Vec3f(1.0f0) 			=> Vec3f(1)
		@fact Vec(1.0, 1.0, 1.0) 	=> Vec3d(1)
		@fact Vec2d(Vec3d(1)) 		=> Vec(1.0, 1.0)
		@fact Vec(Vec3d(1), 1.0) 	=> Vec4d(1)
		@fact Vec(Vec3d(1), 1) 		=> Vec4d(1)
		@fact Vec3d(Vec3f(1.0)) 	=> Vec3d(1.0)
	end
end
v2 = Vec(6.0,5.0,4.0)
v1 = Vec(1.0,2.0,3.0)
v2 = Vec(6.0,5.0,4.0)
	
facts("Indexing") do
	
	context("FixedVector") do 
		@fact v1[1] => 1.0
		@fact v1[2] => 2.0
		@fact v1[3] => 3.0
		@fact_throws BoundsError v1[-1]
		@fact_throws BoundsError v1[0]
		@fact_throws BoundsError v1[4]
	end

end


facts("Ops") do
	context("Negation") do 
		@fact -v1 => Vec(-1.0,-2.0,-3.0)
		@fact isa(-v1, Vec3d) => true
	end

	context("Negation") do 
		@fact v1+v2 => Vec3d(7.0,7.0,7.0)
	end
	context("Negation") do 
		@fact v2-v1 => Vec3d(5.0,3.0,1.0)
	end
	context("Multiplication") do 
		@fact v1.*v2 => Vec3d(6.0,10.0,12.0)
	end
	context("Division") do 
		@fact v1 ./ v1 => Vec3d(1.0,1.0,1.0)
	end

	context("Scalar") do 
		@fact 1.0 + v1 => Vec3d(2.0,3.0,4.0)
		@fact 1.0 .+ v1 => Vec3d(2.0,3.0,4.0)
		@fact v1 + 1.0 => Vec3d(2.0,3.0,4.0)
		@fact v1 .+ 1.0 => Vec3d(2.0,3.0,4.0)
		@fact 1 + v1 => Vec3d(2.0,3.0,4.0)
		@fact 1 .+ v1 => Vec3d(2.0,3.0,4.0)
		@fact v1 + 1 => Vec3d(2.0,3.0,4.0)
		@fact v1 .+ 1 => Vec3d(2.0,3.0,4.0)

		@fact v1 - 1.0 => Vec3d(0.0,1.0,2.0)
		@fact v1 .- 1.0 => Vec3d(0.0,1.0,2.0)
		@fact 1.0 - v1 => Vec3d(0.0,-1.0,-2.0)
		@fact 1.0 .- v1 => Vec3d(0.0,-1.0,-2.0)
		@fact v1 - 1 => Vec3d(0.0,1.0,2.0)
		@fact v1 .- 1 => Vec3d(0.0,1.0,2.0)
		@fact 1 - v1 => Vec3d(0.0,-1.0,-2.0)
		@fact 1 .- v1 => Vec3d(0.0,-1.0,-2.0)

		@fact 2.0 * v1 => Vec3d(2.0,4.0,6.0)
		@fact 2.0 .* v1 => Vec3d(2.0,4.0,6.0)
		@fact v1 * 2.0 => Vec3d(2.0,4.0,6.0)
		@fact v1 .* 2.0 => Vec3d(2.0,4.0,6.0)
		@fact 2 * v1 => Vec3d(2.0,4.0,6.0)
		@fact 2 .* v1 => Vec3d(2.0,4.0,6.0)
		@fact v1 * 2 => Vec3d(2.0,4.0,6.0)
		@fact v1 .* 2 => Vec3d(2.0,4.0,6.0)

		@fact v1 / 2.0 => Vec3d(0.5,1.0,1.5)
		@fact v1 ./ 2.0 => Vec3d(0.5,1.0,1.5)
		@fact v1 / 2 => Vec3d(0.5,1.0,1.5)
		@fact v1 ./ 2 => Vec3d(0.5,1.0,1.5)

		@fact 12.0 ./ v1 => Vec3d(12.0,6.0,4.0)
		@fact 12 ./ v1 => Vec3d(12.0,6.0,4.0)

		@fact (v1 .^ 2) => Vec3d(1.0,4.0,9.0)
		@fact (v1 .^ 2.0) => Vec3d(1.0,4.0,9.0)
		@fact (2.0 .^ v1) => Vec3d(2.0,4.0,8.0)
		@fact (2 .^ v1) => Vec3d(2.0,4.0,8.0)
	end
end


# vector norm
@fact norm(Vec3d(1.0,2.0,2.0)) => 3.0

# cross product
@fact cross(v1,v2) => Vec3d(-7.0,14.0,-7.0)
@fact isa(cross(v1,v2),Vec3d)  => true


# type conversion
@fact isa(convert(Vec3f,v1), Vec3f)  => true

@fact isa(convert(Vector{Float64}, v1), Vector{Float64})  => true
@fact convert(Vector{Float64}, v1) => [1.0,2.0,3.0]


# matrix operations

#typealias Mat1d Matrix1x1{Float64}
typealias Mat2d Mat{2,2, Float64}
typealias Mat3d Mat{3,3, Float64}
typealias Mat4d Mat{4,4, Float64}
zeromat = Mat2d((0.0,0.0),(0.0,0.0))



@fact length(Mat2d) => 4
@fact length(zeromat) => 4

@fact size(Mat2d) => (2,2)
@fact size(zeromat) => (2,2)

@fact zero(Mat2d) => zeromat

for i=1:4, j=1:4
	x1 = rand(i,j)
	@fact Mat(x1') => Mat(x1)'
end



v = Vec(1.0,2.0,3.0,4.0)
r = row(v)
c = column(v)

#@fact r' => c
#@fact c' => r

a = c*r

b = Mat(
	(1.0,2.0,3.0,4.0),
	(2.0,4.0,6.0,8.0),
	(3.0,6.0,9.0,12.0),
	(4.0,8.0,12.0,16.0)
)

@fact length(b) => 16

@fact a=>b
mat30 = Mat(((30.0,),))
@fact r*c => mat30


#@fact row(r, 1) => v
#@fact column(c,1) => v
#@fact row(r+c',1) => 2*v
@fact sum(r) => sum(v)
@fact prod(c) => prod(v)

@fact eye(Mat3d) => Mat((1.0,0.0,0.0),
							(0.0,1.0,0.0),
							(0.0,0.0,1.0))
#@fact v*eye(Mat4d)*v => 30.0
@fact -r => -1.0*r
#@fact diag(diagm(v)) => v

# type conversion
#@fact isa(convert(Matrix1x4{Float32},r),Matrix1x4{Float32})
jm = rand(4,4)
im = Mat(jm)
for i=1:4*2
	@fact jm[i] => im[i]
end
#im = Matrix4x4(jm)
@fact isa(im, Mat4d)  => true

jm2 = convert(Array{Float64,2}, im)
@fact isa(jm2, Array{Float64,2})  => true
@fact jm => jm2

#Single valued constructor
Mat4d(0.0) => zeros(Mat4d)
a = Vec4d(0)
b = Vec4d(0,0,0,0)
@fact a => b

v = rand(4)
m = rand(4,4)
vfs = Vec(v)
mfs = Mat(m)
function Base.isapprox{FSA <: FixedArray}(a::FSA, b::Array)
	for i=1:length(a)
		!isapprox(a[i], b[i]) && return false
	end
	true
end
facts("Matrix Math") do
	for i=1:4, j=1:4
		v = rand(j)
		m = rand(i,j)
		vfs = Vec(v)
		mfs = Mat(m)
		
		context("Matrix{$i, $j} * Vector{$j}") do
			vm = m * v
			fsvm = mfs * vfs
			@fact isapprox(fsvm, vm)  => true
		end
		if i == j
			context("Matrix{$i, $j} * Matrix{$i, $j}") do
				mm = m * m
				#fmm = mfs * mfs
				#@fact isapprox(fmm, mm)  => true
			end
			context("det(M)") do
				mm = det(m)
				#fmm = det(mfs)
				#@fact isapprox(fmm, mm)  => true
			end
			context("inv(M)") do
				mm = inv(m)
				#fmm = inv(mfs)
				#@fact isapprox(fmm, mm)  => true
			end
		end
		
		context("transpose M") do
			mm = m'
			#fmm = mfs'
			#@fact isapprox(fmm, mm)  => true
		end
	end
end

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



@fact isapprox(acfs, ac)  => true
@fact isapprox(bcfs, bc)  => true

@fact isapprox(afs, a) => true
@fact isapprox(bfs, b) => true
@fact isapprox(cfs, c) => true

@fact isapprox(dfs, d) => true
@fact isapprox(d2fs, d2) => true
@fact isapprox(ffs, f) => true
@fact isapprox(gfs, g) => true
@fact isapprox(hfs, h) => true
@fact isapprox(ifs, i) => true
@fact isapprox(jfs, j) => true
@fact isapprox(kfs, k) => true
@fact isapprox(lfs, l) => true

# Equality
@fact Vec{3, Int}(1) => Vec{3, Float64}(1)
@fact Vec{2, Int}(1) => not(Vec{3, Float64}(1))
@fact Vec(1,2,3) => Vec(1.0,2.0,3.0)
@fact Vec(1,2,3) => not(Vec(1.0,4.0,3.0))
@fact Vec(1,2,3) => [1,2,3]
@fact Mat((1,2),(3,4)) => Mat((1,2),(3,4))
let
    a = rand(16)
    b = Mat4d(a)
    @fact b => reshape(a, (4,4))
    @fact reshape(a, (4,4)) => b
    @fact b => not(reshape(a, (2,8)))
end

println("SUCCESS")

