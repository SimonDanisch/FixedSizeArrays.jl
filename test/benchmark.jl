module FSABench

using BenchmarkTools
import FixedSizeArrays
import SIMD
include("missing_ops.jl")
# Define a parent BenchmarkGroup to contain our suite
suite = BenchmarkGroup()

function sum_index(a, b)
	r = a[1]
	@inbounds for i=1:length(a)
		r += a[i] + b[i]
	end
	r
end
@generated function splat{N}(a, ::Val{N})
	names = [Symbol("arg_$i") for i=1:N]
	tup = Expr(:tuple, names...)
	expr=Expr(:(=), tup, :a)
	expr=quote
		$expr
		$tup
	end
	expr
end

function rand_vector{X<:Vector}(A::Type{X}, N, T)
	rand(T, N)
end
function rand_vector{X<:FixedSizeArrays.FixedVector}(A::Type{X}, N, T)
	rand(A{N, T})
end
function rand_vector{X<:Tuple}(A::Type{X}, N, T)
	tuple(rand(T, N)...)
end
function rand_vector{X<:VecElement}(A::Type{X}, N, T)
	map(VecElement, tuple(rand(T, N)...))
end
function rand_vector{X<:SIMD.Vec}(A::Type{X}, N, T)
	SIMD.Vec(tuple(rand(T, N)...))
end

constructors = [
	("Vector", Vector), 
	("FSA", FixedSizeArrays.Vec), 
	("SIMD", SIMD.Vec),
	("VecTuple", VecElement),
	("Tuple", Tuple)
]
element_types = [Float32, Float64, Int64, Int32]
binary_funs = [sum_index, +, ./]
mixed_funs = [+, /, *]
unary_funs = [-, sum, prod]
vector_of_vector_funs = [-, sum, prod, mean]

for (name, c) in constructors
	
	for N in (3, 10), T in element_types
		a = rand_vector(c, N, T)
		b = rand_vector(c, N, T)
		splat(a, Val{N}())
		for fun in unary_funs
			(fun)(a)
		end
		for fun in binary_funs
			(fun)(a, b)
		end
		for fun in mixed_funs
			x = rand(T)
			fun(a, x)
		end
		for fun in vector_of_vector_funs
			ET = typeof(rand_vector(c, N, T))
			a = ET[rand_vector(c, N, T) for _=1:1000]
			@assert eltype(a) == ET
			(fun)(a)
		end
	end
end
println("Okay gang, lets start benchmarkin'")
for (name, c) in constructors
	suite[name] = BenchmarkGroup()
	for N in (3, 10), T in element_types
		a = rand_vector(c, N, T)
		b = rand_vector(c, N, T)

		suite[name][splat, N, T] = @benchmarkable splat($a, $(Val{N}()))

		for fun in unary_funs
			suite[name][fun, N, T] = @benchmarkable $(fun)($a)
		end
		for fun in binary_funs
			suite[name][fun, N, T] = @benchmarkable $(fun)($a, $b)
		end
		for fun in mixed_funs
			x = rand(T)
			suite[name]["mixed", fun, N, T] = @benchmarkable $(fun)($a, $x)
		end
		for fun in vector_of_vector_funs
			ET = typeof(rand_vector(c, N, T))
			a = ET[rand_vector(c, N, T) for _=1:50]
			@assert eltype(a) == ET
			suite[name]["vecvec", fun, N, T] = @benchmarkable $(fun)($a)
		end
	end
end
tune!(suite)
const result = run(suite)
end
