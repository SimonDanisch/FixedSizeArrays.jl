@inline function reduce{FSA <: Union{FixedArray, Tuple}}(f, a::FSA)
    length(a) == 1 && return a[1]
    @inbounds begin
        red = f(a[1], a[2])
        for i=3:length(a)
            red = f(red, a[i])
        end
    end
    red
end

@inline function reduce(f::Base.Func{2}, a::Mat)
    length(a) == 1 && return a[1,1]
    @inbounds begin
        red = reduce(f, a.(1)[1])
        for i=2:size(a, 2)
            red = f(red, reduce(f, a.(1)[i]))
        end
    end
    red
end

index_expr{T <: Number}(::Type{T}, i::Int, inds::Int...) = :($(symbol("arg$i")))
index_expr{T <: FixedArray}(::Type{T}, i::Int, inds::Int...) = :($(symbol("arg$i"))[$(inds...)])
inner_expr{N}(args::NTuple{N, DataType}, inds::Int...) = :( F($(ntuple(i -> index_expr(args[i], i, inds...), N)...)) )


# This solves the combinational explosion from FixedVectorNoTuple while staying fast.
function constructor_expr{T <: FixedVector}(::Type{T}, tuple_expr::Expr)
    quote
        $(Expr(:boundscheck, false))
        $(Expr(:meta, :inline))
        rvalue = FSA($(tuple_expr))
        $(Expr(:boundscheck,:pop))
        rvalue
    end
end
constructor_expr{T <: Mat}(::Type{T}, tuple_expr::Expr) = quote
    $(Expr(:boundscheck, false))
    $(Expr(:meta, :inline))
    rvalue = Mat($(tuple_expr))
    $(Expr(:boundscheck,:pop))
    rvalue
end
constructor_expr{T <: FixedVectorNoTuple}(::Type{T}, tuple_expr::Expr) = quote
    $(Expr(:boundscheck, false))
    $(Expr(:meta, :inline))
    rvalue = FSA($(tuple_expr)...)
    $(Expr(:boundscheck,:pop))
    rvalue
end
@generated function map{FSA <: FixedArray}(F, arg1::FSA, arg2::FSA)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F, arg1::FSA, arg2::Number)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F, arg1::Number, arg2::FSA)
    inner = fill_tuples_expr((inds...) -> inner_expr((arg1, arg2), inds...), size(FSA))
    constructor_expr(FSA, inner)
end
@generated function map{FSA <: FixedArray}(F, ::Type{FSA})
    inner = fill_tuples_expr((inds...) -> :(F($(inds...))), size(FSA))
    :( FSA($inner) )
end

@generated function map{FSA <: FixedArray}(F, arg1::FSA)
    inner = fill_tuples_expr((inds...) -> :( F(arg1[$(inds...)]) ), size(FSA))
    constructor_expr(FSA, inner)
end

qualified_name{T}(::Type{T}) = parse(string("Main.", T.name))# is this actually a good idea?

@generated function similar{FSA <: FixedVector}(::Type{FSA}, ElType::DataType, s=size(FSA))
    name = qualified_name(FSA)
    :($name{s[1], ElType, $(FSA.parameters[3:end]...)})
end
similar{FSA <: Mat}(::Type{FSA}, ElType::DataType, s=Tuple{Int, Int}) = Mat{s..., ElType}
similar{FSA <: Mat}(::Type{FSA}, ElType::DataType, s=Tuple{Int}) = Mat{s[1],s[1], ElType}
similar{FSA <: Mat}(::Type{FSA}, ElType::DataType) = similar(FSA, ELType, size(FSA))

@generated function similar{FSA <: FixedVectorNoTuple}(::Type{FSA}, ElType::DataType)
    name = qualified_name(FSA)
    :($name{ElType, $(FSA.parameters[3:end]...)})
end

@generated function map{T}(::Type{T}, arg1::FixedArray)
    eltype(arg1) == T && return :(arg1)
    FSA = similar(arg1, T)
    inner = fill_tuples_expr((inds...) -> :( T(arg1[$(inds...)]) ), size(FSA))
    :( $FSA($(inner)) )
end

@inline map{R,C,T}(F::Type{T}, arg1::Mat{R,C,T}) = arg1
@generated function map{R,C,T}(F::DataType, arg1::Mat{R,C,T})
    inner = fill_tuples_expr((inds...) -> :( F(arg1[$(inds...)]) ), (R, C))
    :( Mat{R, C, F}($(inner)) )
end

@generated function map{FSA <: FixedVectorNoTuple}(F::DataType, arg1::FSA)
    eltype(FSA) == F && return :(arg1)
    inner = ntuple(i-> :(F(arg1[$i])), length(FSA))
    :( similar(FSA, F)($(inner...)) )
end
