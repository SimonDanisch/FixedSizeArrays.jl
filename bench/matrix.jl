using ImmutableArrays, FixedSizeArrays



function test(N, a, b)
	result = a
	for i=1:N
		result = a*b
	end
	return result
end



function bench(N)

	const AB 	= Float64[Float64(i*j) for i=1:4, j=1:4]
	const ABV 	= Float64[1.0, 1.2,2.2, 3.4]
	const AI 	= Matrix4x4{Float64}(AB)
	const AFS 	= nvec(AB)
	const FSV 	= nvec(ABV...)
	AFS * FSV
	@show AB * ABV

	@show AFS'
	@show AB'

	@show AB'
	bmul  = AB*AB
	imul  = AI*AI
	fsmul = matmulfs(AFS,AFS)
	println(@code_llvm matmulfs(AFS,AFS))

	for i=1:4, j=1:4
		@assert bmul[i,j] == imul[i,j] == fsmul[i,j]
	end

	println("julia:")
	@time result = test(1, AB, AB)
	@time result = test(N, AB, AB)

	println("Immutables:")
	@time result = test(1, AI, AI)
	@time result = test(N, AI, AI)
	
	println("Fsa:")
	@time result = test(1, AFS, AFS)
	@time result = test(N, AFS, AFS)
end
bench(10^7)