immutable Matrix4x4{T}
	x1::T
	x2::T
	x3::T
	x4::T
	x5::T
	x6::T
	x7::T
	x8::T
	x9::T
	x10::T
	x11::T
	x12::T
	x13::T
	x14::T
	x15::T
	x16::T
end
@inline Base.getindex(x::Matrix4x4, i) = getfield(x, i)
@generated function (*){T}(a::Matrix4x4{T}, b::Matrix4x4{T})
	M,N,K = 4,4,4
    :(Matrix4x4(
        $([:(
            +($(ntuple(k->:(a[$(sub2ind((M,N), k, i))]*b[$(sub2ind((N,K), k, j))]), N)...))
        ) for i=1:M, j=1:K]...)
    ))
end



function test(a)
    r = a
    for i=1:10^6
        r = r*a
    end
    r
end

function test1()
	const a = Matrix4x4(
	    1,2,3,4,
	    1,2,3,4,
	    1,2,3,4,
	    1,2,3,4,
	)

	const b = [
	    1 2 3 4;
	    1 2 3 4;
	    1 2 3 4;
	    1 2 3 4;
	]
	@time test(a)
	@time test(a)

	@time test(b)
	@time test(b)
end

test1()