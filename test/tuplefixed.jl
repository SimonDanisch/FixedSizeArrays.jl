using GeometryTypes

function test(a)
    r = a
    for i=1:10^5
        r = r*a
    end
    r
end
function test1()
	const a = Matrix4x4{Int}(
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