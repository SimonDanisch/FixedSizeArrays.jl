
function reduce{FSA <: FixedArray}(f::Func{2}, a::FSA)
    red = f(a[1], a[2])
    @inbounds for i=3:length(a)
        red = f(red, a[i])
    end
    red
end
function Base.reduce(f::Base.Func{2}, a::FixedMatrix)
    red = reduce(f, a.(1)[1])
    @inbounds for i=2:size(a, 2)
        red = f(red, reduce(f, a.(1)[i]))
    end
    red
end

function reduce{FSA <: FixedArray}(f::Func{2}, a::FSA)
    red = f(a[1], a[2])
    for i=3:length(a)
        red = f(red, a[i])
    end
    red
end

index_expr{T <: Number}(::Type{T}, i::Int, inds::Int...) = :($(symbol("arg$i")))
index_expr{T <: FixedArray}(::Type{T}, i::Int, inds::Int...) = :($(symbol("arg$i"))[$(inds...)])

inner_expr{N}(args::NTuple{N, DataType}, inds::Int...) = :( F($(ntuple(i -> index_expr(args[i], i, inds...), N)...)) )

# This solves the combinational explosion from FixedVectorNoTuple while staying fast.
function constructor_expr{T <: FixedArray}(::Type{T}, tuple_expr::Expr)
    BaseName = T.name.name
    :(Main.$BaseName($tuple_expr))
end
function constructor_expr{T <: FixedVectorNoTuple}(::Type{T}, tuple_expr::Expr)
    BaseName = T.name.name
    :(Main.$BaseName($(tuple_expr)...))
end

@generated function map{FSA <: FixedArray}(F::Func{2}, arg1::FSA, arg2::FSA)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F::Func{2}, arg1::FSA, arg2::Number)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F::Func{2}, arg1::Number, arg2::FSA)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F::Func, ::Type{FSA})
    inner = fill_tuples_expr((inds...) -> :(F($(inds...))), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F::Func{1}, arg1::FSA)
    inner = fill_tuples_expr((inds...) -> :( F(arg1[$(inds...)]) ), size(FSA))
    constructor_expr(FSA, inner)
end
