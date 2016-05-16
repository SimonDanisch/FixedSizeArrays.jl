# workaround to return all variables, so that nothing get's optimized away,
# but don't use an array, since it will sprinkle the code with any's and grow_array

const syms = []

macro save(x)
    s = esc(gensym())
    push!(syms, s)
    :($s = $x)
end
macro ret()
    :(tuple($(syms...)))
end
immutable TestFunctor end
@compat (::TestFunctor)(a, b) = a+b

function use_operations()

    a = Vec{3, Float32}(0)
    b = Vec{3, Float32}(0.0)
    c = Vec{3, Float32}(0f0)

    a2 = Point{3}(0)
    b2 = Point{3}(0.0)
    c2 = Point{3}(0f0)

    a3 = Point(0f0, 1, 1.)
    b3 = Point{3, Float32}(0f0, 1, 1.)
    c3 = Point{3, Float32}((0f0, 1, 1.))

    m1 = @fsa([1 2;3 4;5 6])
    m2 = eye(Mat{3,3,Float32})

    @save dot(a, a)
    @save dot(a2, a2)

    @save a./b
    @save a.*b
    @save a == b
    @save m2*c
    @save m1.+m1
    @save m1./m1
    @save m2.*m2

    @save reduce(TestFunctor(), a3)
    @save reduce(TestFunctor(), m2)

    @save map(TestFunctor(), a3, a3)
    @save map(TestFunctor(), m2, m2)

    @save map(Float32, a3)
    @save map(Float64, m2)

    @save map(Float64, a3)
    @save map(Float32, m2)

    @ret
end

# This seems to be the easiest way to test,
# that all variables even if they come from inlining,
# are inferred as a concrete type.
# seems like it's hard to turn of inlining when doing coverage, though...
# so this must be executed locally to test for it.
context("type inference") do
    io = IOBuffer()
    use_operations()
    code_warntype(io, use_operations, ())
    str = takebuf_string(io)
    x = matchall(r"Any", str)
    @fact length(x) --> 0
end
