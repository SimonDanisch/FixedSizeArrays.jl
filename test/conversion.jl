using GeometryTypes
using FixedSizeArrays



function test1()
	a=Vector4(1,2,3,4)
	rsult = a[Tuple{1:2:4}]
	for i=1:10^6
		rsult = a[Tuple{1:2:4}]
	end
	rsult
end

function test2()
	a=Vector4(1,2,3,4)
	rsult = FixedSizeArrays.Vec2{Int64}(a[1], a[3])
	for i=1:10^6
		rsult = FixedSizeArrays.Vec2{Int64}(a[1], a[3])
	end
	rsult
end
function test6()
	a=(1,2,3,4)
	rsult = a[1:2:4]
	for i=1:10^6
		rsult = a[1:2:4]
	end
	rsult
end
@time test1()
@time test1()
@time test2()
@time test2()
@time test6()
@time test6()
@show a=Vector4(1,2,3,4)
const b = (1,2,3,4)
println("################")
@code_llvm b[1:2]

@show a[Tuple{1:2:4}]
@show Vector2(Vector3(1,2,3))


function test3()
	a=Vector4(1,2,3,4)
	rsult = Vector2(a)
	for i=1:10^6
		rsult = Vector2(a)
	end
	rsult
end
function test4()
	a=Vector4(1,2,3,4)
	rsult = Vector2(a[1], a[2])
	for i=1:10^6
		rsult = Vector2(a[1], a[2])
	end
	rsult
end
function test5()
	a=Vector4(1,2,3,4)
	rsult = Vector2{Int}(a)
	for i=1:10^6
		tic()
		rsult = Vector2{Int}(a)
