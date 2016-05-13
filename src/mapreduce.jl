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

@inline function reduce(f::Functor{2}, a::Mat)
    length(a) == 1 && return a[1,1]
    @inbounds begin
        red = reduce(f, a._[1])
        for i=2:size(a, 2)
            red = f(red, reduce(f, a._[i]))
        end
    end
    red
end

# Get expression indexing collection `name` of type `T` for use in map()
index_expr{T <: Number}(::Type{T},     name, inds::Int...) = :($name)
index_expr{T <: FixedArray}(::Type{T}, name, inds::Int...) = :($name[$(inds...)])
index_expr{T <: Array}(::Type{T}, name, inds::Int...) = :($name[$(inds...)])

# Get expression checking size of collection `name` against `SIZE`
sizecheck_expr{T <: Number}(::Type{T}, name, SIZE) = :nothing
function sizecheck_expr{FSA<:FixedArray}(::Type{FSA}, name, SIZE)
    size(FSA) == SIZE || :(throw(DimensionMismatch(string($FSA)*" is wrong size")))
    :nothing
end
function sizecheck_expr{A<:Array}(::Type{A}, name, SIZE)
    quote
        # Note - should be marked with @boundscheck in 0.5
        size($name) == $SIZE || throw(DimensionMismatch(string($A)*" is wrong size"))
    end
end

# Get expression to construct FSA from a nested tuple
constructor_expr{FSA <: FixedArray}(::Type{FSA}, tuple_expr::Expr) = :($FSA($tuple_expr))
constructor_expr{FSA <: FixedVectorNoTuple}(::Type{FSA}, tuple_expr::Expr) = :($FSA($tuple_expr...))

# Generate an unrolled nested tuple of calls mapping funcname across the input,
# constructing an OutFSA to store the result.
#
# julia> FixedSizeArrays.unrolled_map_expr(:f, Vec{2,Bool}, (Vec{2,Int},Int), (:A,:b))
#
# generates, after cleaning up:
#
#    (FixedSizeArrays.Vec{2,Bool})(
#        tuple(tuple(f(A[1,1],b), f(A[2,1],b)),
#              tuple(f(A[1,2],b), f(A[2,2],b)))
#    )
function unrolled_map_expr(funcname, OutFSA, argtypes, argnames)
    SIZE = size(OutFSA)
    sizecheck = [sizecheck_expr(T,n,SIZE) for (T,n) in zip(argtypes,argnames)]
    tuple_expr = fill_tuples_expr(SIZE) do inds...
        Expr(:call, funcname,
            [index_expr(argtypes[i], argnames[i], inds...) for i=1:length(argtypes)]...
        )
    end
    quote
        $(Expr(:meta, :inline))
        $(sizecheck...)
        $(Expr(:boundscheck, false))
        rvalue = $(constructor_expr(OutFSA, tuple_expr))
        $(Expr(:boundscheck,:pop))
        rvalue
    end
end

# Versions of map() with a given `OutFSA` output container as the second
# argument.  Unary and binary versions are written explicitly since this
# generates better code in 0.4.
#
# TODO: Nullary version is inconsistent, since it maps indices through the
# functor.
@generated function map{OutFSA<:FixedArray}(func, ::Type{OutFSA}, arg1)
    unrolled_map_expr(:func, OutFSA, (arg1,), (:arg1,))
end
@generated function map{OutFSA<:FixedArray}(func, ::Type{OutFSA}, arg1, arg2)
    unrolled_map_expr(:func, OutFSA, (arg1,arg2), (:arg1,:arg2))
end
@generated function map{OutFSA<:FixedArray}(func, ::Type{OutFSA}, args...)
    argexprs = ntuple(i->:(args[$i]), length(args))
    unrolled_map_expr(:func, OutFSA, args, argexprs)
end


# Unary version
@inline map{FSA <: FixedArray}(F, arg1::FSA) = map(F, FSA, arg1)

immutable ConstructTypeFun{T}; end
call{T}(::ConstructTypeFun{T}, x) = T(x)

# Unary version for type conversion.  Need to override this explicitly for
# output type inference and to prevent conflict with Base.
@inline map{T,N,S}(::Type{T}, arg1::FixedArray{T,N,S}) = arg1 # nop version
@inline map{T,FSA<:FixedArray}(::Type{T}, arg1::FSA) = map(ConstructTypeFun{T}(), similar(FSA, T), arg1)

# Nullary special case version for constructing FSAs
# TODO: This is inconsistent with the above
@generated function map{FSA <: FixedArray}(F, ::Type{FSA})
    tuple_expr = fill_tuples_expr((inds...) -> :(F($(inds...))), size(FSA))
    constructor_expr(FSA, tuple_expr)
end

#immutable InferFSAConstructor{FSA}; end

# Type-inferred versions of map()
function unrolled_map_expr2(funcname, templateFSA, SIZE, argtypes, argnames)
    sizecheck = [sizecheck_expr(T,n,SIZE) for (T,n) in zip(argtypes,argnames)]
    tuple_expr = fill_tuples_expr(SIZE) do inds...
        Expr(:call, funcname,
            [index_expr(argtypes[i], argnames[i], inds...) for i=1:length(argtypes)]...
        )
    end
    quote
        $(Expr(:meta, :inline))
        $(sizecheck...)
        $(Expr(:boundscheck, false))
        rvalue = construct_similar($templateFSA, $tuple_expr)
        $(Expr(:boundscheck,:pop))
        rvalue
    end
end

@generated function map{F<:FixedArray}(func, arg1::F)
    unrolled_map_expr2(:func, arg1, size(F), (arg1,), (:arg1,))
end
@generated function map{F1<:FixedArray, F2<:FixedArray}(func, arg1::F1, arg2::F2)
    unrolled_map_expr2(:func, arg1, size(F1), (arg1,arg2), (:arg1,:arg2))
end
@generated function map{F<:FixedArray}(func, arg1::F, arg2::Union{Number,AbstractArray})
    unrolled_map_expr2(:func, arg1, size(F), (arg1,arg2), (:arg1,:arg2))
end
@generated function map{F<:FixedArray}(func, arg1::Union{Number,AbstractArray}, arg2::F)
    unrolled_map_expr2(:func, arg2, size(F), (arg1,arg2), (:arg1,:arg2))
end

