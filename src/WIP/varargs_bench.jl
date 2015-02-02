immutable Vec4{T}
	a::T
	b::T
	c::T
	d::T
end
function sub2ind_test(a...)
	sub2ind((10,10))

function test(data...)
	Vec4(data...)
end

function test2(N)
	result = test(rand(Int, 4)...)
	for i=1:N
		result = test(i, 0,0,0)
	end
	result
end
function test3(N)
	result = Vec4(rand(Int, 4)...)
	for i=1:N
		result = Vec4(i, 0,0,0)
	end
	result
end

@time test2(1)
@time test2(10^8)

@time test3(1)
@time test3(10^8)