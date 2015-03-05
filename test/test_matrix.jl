using FixedSizeArrays

type A{T} <: FixedVector{T, 3}
	x::T
	y::T 
	z::T
end

function test()
	a = [A(1,2,3), A(1,2,3)]
	r = a[1]
	println(a.data)
	a = 0
	r
end

a = test()
gc()
sleep(1)
gc()
@show a