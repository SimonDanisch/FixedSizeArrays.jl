module Test

export test

function test()
	a = eval(Main, :RGB)
	a(1,2,3)
end

end


module Test2
export RGB

immutable RGB{T}
	x::T
	y::T
	z::T
end

end
using Test
using Test2


test()