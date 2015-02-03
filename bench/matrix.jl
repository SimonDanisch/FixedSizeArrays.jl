using ImmutableArrays, FixedSizeArrays



const AB 	= Float32[float32(i*j) for i=1:4, j=1:4]
const AI 	= Matrix4x4{Float32}(AB)
const AFS 	= nvec(AB)

bmul  = AB*AB
imul  = AI*AI
fsmul = AFS*AFS

for i=1:4, j=1:4
	@assert bmul[i,j] == imul[i,j] == fsmul[i,j]
end

function test(N, a, b)
	result = a
	for i=1:N
		result = a*b
	end
	return result
end

@time test(1, AB, AB)
@time test(10^5, AB, AB)

@time test(1, AI, AI)
@time test(10^5, AI, AI)

@time test(1, AFS, AFS)
@time test(10^5, AFS, AFS)
@time test(10^5, AFS, AFS)
@time test(10^5, AFS, AFS)

AI*AI