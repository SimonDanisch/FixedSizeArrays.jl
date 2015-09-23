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
	tms 	= gensym()
	bts 	= gensym()
	gctms 	= gensym()
	gcps 	= gensym()
	gcswps 	= gensym()

	tms1 	= gensym()
	gcswps1 = gensym()
	gcps1 	= gensym()
	gctms1 	= gensym()
	bts1 	= gensym()
	esc(quote
		$bts 	 = Base.gc_bytes()
		$gctms 	 = Base.gc_time_ns()
		$gcps 	 = Base.gc_num_pause()
		$gcswps  = Base.gc_num_full_sweep()
		$tms 	 = Base.time_ns()

		$expr

		$tms1 	 = Base.time_ns()
		$gcswps1 = Base.gc_num_full_sweep()
		$gcps1 	 = Base.gc_num_pause()
		$gctms1  = Base.gc_time_ns()
		$bts1 	 = Base.gc_bytes()

		BenchData($bts1-$bts, $tms1-$tms, $gctms1-$gctms, $gcps1-$gcps, $gcswps1-$gcswps)
	end)
end
type Vec4{T} <: FixedVector{T, 4}
	x::T
	y::T
	z::T
	w::T
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

#=
function test2(N)
	a,b,c,d = rand(Float64, 4)
	result =  LolVec(a,b,c,d)
	bench = 0.0
	for i=1:N
		a,b,c,d = rand(Float64, 4)
		tic()
		result = LolVec(a,b,c,d)
		b = toq()
		bench += b
	end

	result,bench
end

@time t= test2(10)
@time t,b= test2(10^4)
println(b)
#0.00622337 seconds (1 MB allocated)
#0.001912617 seconds (470 kB allocated)

#0.0036441109999999725 2mb
#0.0010228019999999729 1mb
#0.0012369930000000022 1mb
=#