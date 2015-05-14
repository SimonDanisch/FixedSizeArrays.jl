#simple reduce function
function reduce(f::Func{2}, a::FixedArray)
    alength = length(a)
    alength == 1 && return a[1]
    s = f(a[1], a[2])
    for i=3:alength
        s = f(s, a[i])
    end
    s
end

# either gives an expression to index into the fixedsizearray, or just inserts the constant
accessor{T <: FixedArray}(arg::Type{T}, i::Integer, name::Symbol) = :($name[$i])
accessor{T                      }(arg::Type{T}, i::Integer, name::Symbol) = :($name)


#To simplify things, map_expression assumes that the argument names of the staged functions are f for the callable and a,b,c,.. for the other arguments.
function map_expression{F <: Func, FSA <: FixedArray}(f::Type{F}, fsa::Type{FSA}, args...)
    Cardinality = first(super(F).parameters)
    @assert length(args) == Cardinality "cardinality of callable doesn't match arguments. Callable: $Cardinality, args: $(length(args))"
    argnames = [:a, :b, :c, :d, :e, :f, :g] # must be the same as the arguments of map
    #call f for every index in FSA, with all the args.
    expr = ntuple(length(FSA)) do i 
        quote
            f( $(ntuple(j-> accessor(args[j], i, argnames[j]), Cardinality)...)) 
        end
    end
    # the type FSA from map can't be used, because it contains the parameter already (e.g. FSA{Int}), which results in a conversion.
    # map should return a fixedsize array with the return type of f, so the parameter has to be stripped of FSA.
    type_name = :($(symbol(name(fsa)))) # custom types are not known in the Module FixedSizeArray, but in Main
    #eval(type_name)
    :(FSA($(expr...)))
end

@generated function map{FSA <: FixedArray, F <: Func{1}}(f::Union(Type{F}, F), a::Type{FSA})
    :(FSA($([:(f($i)) for i=1:length(FSA)]...)))
end
@generated function map{FSA <: FixedArray, F <: Func{2}}(f::Union(Type{F}, F), a::Type{FSA})
    #@assert ndims(FSA) == 2
    :(FSA($([:(f($i, $j)) for i=1:size(FSA,1), j=1:size(FSA,2)]...)))
end

@generated function map{FSA <: FixedArray}(f::Func{1}, a::FSA)
    map_expression(f, a, a)
end
@generated function map{FSA <: FixedArray}(f::Func{2}, a::FSA, b::FSA)
    map_expression(f, a, a, b)
end
@generated function map{FSA <: FixedArray, REAL <: Real}(f::Func{2}, a::FSA, b::REAL)
    map_expression(f, a, a, b)
end
@generated function map{FSA <: FixedArray, REAL <: Real}(f::Func{2}, a::REAL, b::FSA)
    map_expression(f, b, a, b)
end



@generated function product(f, it...)
    variables       = ntuple(fieldname, length(it))
    iterator_access = Expr(:block, [:($(variables[i]) = it[$i]) for i=1:length(it)]...)
    for_expression  = Expr(:for, 
        iterator_access, 
        :(f($(variables...)))
    )
    quote
    $for_expression
    nothing
    end
end

