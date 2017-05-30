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


"""
    is_inferred(ex)

Test whether type inference has fully inferred all types in the AST `ex`, more
or less as done in `code_warntype()`.  That is, each AST node should have a
concrete type (`isleaftype()` should be `true`), except for some special
cases.
"""
function is_inferred(ex::Expr)
    inferred = true
    # We don't expect a type for certian `Expr`s - see the show_type deduction
    # in Base.show_unquoted() which is called from code_warntype()
    check_type = !(
        ex.head == :(=) ||
        ex.head == :boundscheck ||
        ex.head == :gotoifnot ||
        ex.head == :return ||
        (ex.head == :call && (in(ex.args[1], (GlobalRef(Base, :box), TopNode(:box), :throw)) ||
                             Base.ismodulecall(ex) ||
                             (ex.typ === Any && Base.is_intrinsic_expr(ex.args[1]))))
    )
    if check_type
        inferred &= isleaftype(ex.typ)
    end
    # Avoid traversing some expressions, since we don't care about their types.
    # It's not quite clear how to detect these reliably, since
    # Base.show_unquoted() has them in the else part of a large set of elseif
    # clauses.
    check_expr_args = !(
        ex.head == :boundscheck
    )
    if check_expr_args
        for arg in ex.args
            inferred &= is_inferred(arg)
        end
    end
    return inferred
end

# Assume AST nodes which aren't Exprs don't have meaningful type information
is_inferred(nonexpr) = true


# This seems to be the easiest way to test,
# that all variables even if they come from inlining,
# are inferred as a concrete type.
# seems like it's hard to turn of inlining when doing coverage, though...
# so this must be executed locally to test for it.
context("type inference") do
    ct = code_typed(use_operations, ())
    # Only inspect body of `use_operations()`
    body = ct[1].args[3]

    for arg in body.args
        @fact is_inferred(arg) --> true   "Failed type inference:  $arg"
    end
end
