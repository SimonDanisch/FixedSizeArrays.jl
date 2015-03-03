reduce{F <: Func{2}, T}(f::Type{F}, a::FixedVector{T, 1})   = a[1]
reduce{F <: Func{2}, T}(f::Type{F}, a::FixedMatrix{T, 1,1}) = a[1]
stagedfunction reduce(f::Func{2}, a::FixedVector)
    alength = length(a)
    @assert alength >= 2 "type $a not long enough for reduce (needs at least two elements)"
    quote 
        s = f(a[1], a[2])
        $([:(s = f(s, a[$i])) for i=3:alength]...)
    end
end

# either gives an expression to index into the fixedsizearray, or just inserts the constant
accessor{T <: FixedArray}(arg::Type{T}, i::Integer, name::Symbol) = :($name[$i])
accessor{T                      }(arg::Type{T}, i::Integer, name::Symbol) = :($name)


#To simplify things, map_expression assumes that the argument names of the staged functions are f for the callable and a,b,c,.. for the other arguments.
function map_expression{F <: Func, FSA <: FixedArray}(f::Type{F}, fsa::Type{FSA}, args...)
    FN = first(super(F).parameters)
    @assert length(args) == FN "cardinality of callable doesn't match arguments. Callable: $FN, args: $(length(args))"
    argnames = [:a, :b, :c, :d, :e, :f, :g] # must be the same as the arguments of map
    #call f for every index in FSA, with all the args.
    expr = ntuple(length(FSA)) do i 
        quote
            f( $(ntuple(j-> accessor(args[j], i, argnames[j]), FN)...)) 
        end
    end
    # the type FSA from map can't be used, because it contains the parameter already (e.g. FSA{Int}), which results in a conversion.
    # map should return a fixedsize array with the return type of f, so the parameter has to be stripped of FSA.
    type_name = :(Main.$(symbol(name(fsa)))) # custom types are not known in the Module FixedSizeArray, but in Main
    eval(type_name)
    :($type_name($(expr...)))
end



stagedfunction map{FSA <: FixedArray}(f::Func{1}, a::FSA)
    map_expression(f, a, a)
end
stagedfunction map{FSA <: FixedArray}(f::Func{2}, a::FSA, b::FSA)
    map_expression(f, a, a, b)
end
stagedfunction map{FSA <: FixedArray, REAL <: Real}(f::Func{2}, a::FSA, b::REAL)
    map_expression(f, a, a, b)
end
stagedfunction map{FSA <: FixedArray, REAL <: Real}(f::Func{2}, a::REAL, b::FSA)
    map_expression(f, b, a, b)
end



stagedfunction product(f, it...)
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

map{FSA <: FixedArrayWrapper}(f::Func{1}, a::FSA)                        = without_params(FSA)(map(f, a.(1)))
map{FSA <: FixedArrayWrapper}(f::Func{2}, a::FSA, b::FSA)                = without_params(FSA)(map(f, a.(1), b.(1)))
map{FSA <: FixedArrayWrapper, REAL <: Real}(f::Func{2}, a::FSA, b::REAL) = without_params(FSA)(map(f, a.(1), b))
map{FSA <: FixedArrayWrapper, REAL <: Real}(f::Func{2}, a::REAL, b::FSA) = without_params(FSA)(map(f, a, b.(1)))
