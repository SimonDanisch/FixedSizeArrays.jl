
# helper function for creating names.
tuple_to_string(t::(Any...), sep::AbstractString) = foldl((v0,d)-> string(v0)*sep*string(d), t)

# Different base names for different dimensionalities...this is definitely a trade off between commonly used names and homogenity.
vecname(sz::(Integer, Integer, Integer...), mutable::Bool, basename::AbstractString="FSArray") = symbol((mutable?"M":"")*basename*tuple_to_string(sz, "x"))
vecname(sz::(Integer, Integer), mutable::Bool, basename::AbstractString="Matrix")              = symbol((mutable?"M":"")*basename*tuple_to_string(sz, "x"))
vecname(sz::(Integer,), mutable::Bool, basename::AbstractString="Vec")                         = symbol((mutable?"M":"")*basename*string(first(sz)))


function fixedarray_type_expr(typename::Symbol, SIZE::(Integer...), mutable::Bool, fields::Function=fieldname)
    len = prod(SIZE)
    fields_names        = [fields(i) for i=1:len]
    fields_types        = [Expr(:(::), fields(i), :T) for i=1:len]
    NDim                = length(SIZE)
    
    typ_expr = mutable ? quote
        type $(typename){T} <: MutableFixedArray{T, $NDim, $SIZE}
            $(fields_types...)
        end
    end : quote 
        immutable $(typename){T} <: FixedArray{T, $NDim, $SIZE}
            $(fields_types...)
        end
    end
end
# maps the dimension of vectors always to (length,)
flatten_dimension(SIZE) = SIZE
flatten_dimension(SIZE::(Int, Int)) = any(x->x==1, SIZE) ? prod(SIZE) : SIZE

function gen_fixedsizevector_type(SIZE::(Integer...), mutable::Bool)
    typename = vecname(SIZE, mutable)
    # if already exists, there's nothing to be done here
    isdefined(FixedSizeArrays, typename) && return typename
    expr = fixedarray_type_expr(typename, SIZE, mutable)
    #Call the type into existence
    eval(FixedSizeArrays, expr)
    #return name of type
    typename
end

#General constructor for arbitrary fixedsizearrays
stagedfunction call{T, NDim, SIZE}(t::Type{FixedArray{T, NDim, SIZE}}, data::T...)
    N = length(data)
    @assert prod(SIZE) == N "not the right dimension"
    !t.abstract && return :(t(data...)) # return if type is known and not abstract
    typename = gen_fixedsizevector_type(SIZE, t.mutable)
    :($typename(data...))
end

nvec{T <: AbstractArray}(x::T)        = FixedArray(x)
nvec{T}(x::T...)                      = FixedArray{T, 1, (length(x),)}(x)
nvec{T}(SIZE::(Integer...,), x::T...) = FixedArray{T, length(SIZE), SIZE}(x)


macro gen_fixed_size_vector(basename, fields, N, mutable)
    # bring the sugar back
    fields = Symbol[elem.args[1] for elem in fields.args]
    fieldfunction = i->fields[i]
    N      = (N.args[1]):(N.args[2])
    expr   = [begin
        typename = vecname((i,), mutable, basename)
        fixedarray_type_expr(typename, (i,), mutable, fieldfunction) 
    end for i in N]
    esc(Expr(:block, expr...))
end


function gen_fixed_size_matrix(M, N, mutable)
    # bring the sugar back
    expr = []
    for i in M, j in N
        typename = vecname((i,j), mutable)
        push!(expr, fixedarray_type_expr(typename, (i,j), mutable))
        push!(expr, Expr(:export, typename))
    end

    eval(FixedSizeArrays, Expr(:block, expr...))
end

immutable RandFunc{T} <: Func{1} 
eltype::Type{T}
end
call{T}(::RandFunc{T}, x) = rand(T)

rand{FSA <: FixedArray}(x::Type{FSA}) = map(RandFunc(eltype(FSA)), FSA)


immutable ZeroFunc{T} <: Func{1}
    typ::Type{T}
end
call{T}(::ZeroFunc{T}, x) = zero(T)
zero{FSA <: FixedArray}(::Type{FSA}) = map(ZeroFunc(eltype(FSA)), FSA)

immutable EyeFunc{NDim} <: Func{1}
    size::NTuple{NDim, Int}
    eltype::DataType
end
function call{T}(ef::EyeFunc{T}, x)
    i,j = ind2sub(ef.size, x)
    i==j?one(ef.eltype) : zero(ef.eltype)
end
eye{FSA <: FixedArray}(::Type{FSA}) = map(EyeFunc(size(FSA), eltype(FSA)), FSA)

function row{T, N}(v::FixedVector{T, N})
    FixedMatrix{T, 1, N}(v...)
end
column{T, N}(v::FixedVector{T, N}) = v
   
#=
    Base.call{AT <: Array}(::Type{$(typename)}, A::AT) = $(typename)($([:(A[$i]) for i=1:len]...))
    $(typename)($(fields_names...)) = $(typename)(promote($(fields_names...))...)
    type $(mutabletypename){T} <: MutableFixedArray{T, $NDim, $SIZE}
        $(fields_types...)
        #$(mutabletypename)($(fields_names...)) = new($(fields_names...))
    end
    Base.call{AT <: Array}(::Type{$(mutabletypename)}, A::AT) = $(mutabletypename)($([:(A[$i]) for i=1:len]...))
    $(mutabletypename)($(fields_names...)) = $(mutabletypename)(promote($(fields_names...))...)
end
=#