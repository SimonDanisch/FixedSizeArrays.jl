using FixedSizeArrays, ImmutableArrays

immutable BenchData <: FixedVector{Float64, 5}
	gc_bytes::Float64
	time_ns::Float64
	gc_time_ns::Float64
	gc_num_pause::Float64
	gc_num_full_sweep::Float64
end
function Base.show(io::IO, b::BenchData)
	print(io, """
		bytes allocated 	: $(b.gc_bytes)
		time_ns 			: $(b.time_ns)
		gc_time_ns 			: $(b.gc_time_ns)
		gc_num_pause 		: $(b.gc_num_pause)
		gc_num_full_sweep 	: $(b.gc_num_full_sweep)
		""")
end

macro timing(expr)
	esc(quote
		bts 	= Base.gc_bytes() # line 59:
		gctms 	= Base.gc_time_ns() # line 61:
		gcps 	= Base.gc_num_pause() # line 62:
		gcswps 	= Base.gc_num_full_sweep() # line 63:
		
		tms 	= Base.time_ns() # line 60:
		val 	= $(expr)
		tms1 	= Base.time_ns() # line 68:

		gcswps1 = Base.gc_num_full_sweep() # line 65:
		gcps1 	= Base.gc_num_pause() # line 66:
		gctms1 	= Base.gc_time_ns() # line 67:
		bts1 	= Base.gc_bytes() # line 69:
		BenchData(bts1-bts, tms1-tms, gctms1-gctms, gcps1-gcps, gcswps1-gcswps)
	end)
end
immutable LolVec{T} <: FixedVector{T, 4}
	x::T
	y::T
	z::T
	w::T
end
immutable Point{FSV <: FixedVector} <: FixedArrayWrapper{FSV}
	vec::FSV
end

function test(N)
	bench = Dict(
		(:julia, :vector_creation) 		=> BenchData[],
		(:julia, :matrix_creation) 		=> BenchData[],
		(:julia, :vector_creation_rand) => BenchData[],
		(:julia, :matrix_creation_rand) => BenchData[],
		(:fsa, 	 :vector_creation) 		=> BenchData[],
		(:fsa, 	 :matrix_creation) 		=> BenchData[],

		(:fsa, 	 		:vector_creation_fromj) => BenchData[],
		(:fsa, 	 		:matrix_creation_fromj) => BenchData[],
		(:fsawrapper, 	:vector_creation) 		=> BenchData[],
		(:immutablearrays, 	:vector_creation) 		=> BenchData[],
	)
	for i=1:N
		push!(bench[(:julia, :vector_creation)], (@timing vecdata = [1,2,3,4]))
		push!(bench[(:julia, :matrix_creation)], (@timing matdata = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]))

		push!(bench[(:julia, :vector_creation_rand)], (@timing vecrand = rand(4)))
		push!(bench[(:julia, :matrix_creation_rand)], (@timing matrand = rand(4,4)))
		nvec(1,2,3,4)
		nvec(1,2,3,4)

		push!(bench[(:fsa, :vector_creation)], (@timing vecrand = nvec(1,2,3,4)))
		push!(bench[(:fsa, :matrix_creation)], (@timing matrand = nvec((4,4), 1,2,3,4, 1,2,3,4, 1,2,3,4, 1,2,3,4)))
		
		push!(bench[(:fsa, :vector_creation_fromj)], (@timing vecrand = nvec(vecrand)))
		push!(bench[(:fsa, :matrix_creation_fromj)], (@timing matrand = nvec(matrand)))

		push!(bench[(:fsawrapper, :vector_creation)], (@timing vecrand = Point(1,2,3,4)))

		push!(bench[(:immutablearrays, :vector_creation)], (@timing vecrand = Vector4(1,2,3,4)))
	end
	bench
end
result = test(1000)

println(mean(result[(:fsa, :vector_creation)][100:end]))
println(mean(result[(:fsawrapper, :vector_creation)][100:end]))
println(mean(result[(:immutablearrays, :vector_creation)][100:end]))

println(@which NVec4(1,2,3,4))
function test2()
	println(@which nvec(1,2,3,4))
	@time a = nvec(1,2,3,4)
	@time a = nvec(1,2,3,4)
	@time a = nvec(1,2,3,4)

	@time c = NVec4(1,2,3,4)
	@time c = NVec4(1,2,3,4)
	@time b = Vector4(1,2,3,4)
	@time b = Vector4(1,2,3,4)
	a,b,c
end

a,b,c = test2()
println(a,b,c)