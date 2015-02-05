
reduce{F <: Func{2}, T}(f::Type{F}, a::AbstractFixedVector{T, 1}) = a

stagedfunction reduce(f::Func{2}, a::AbstractFixedVector)
    alength = length(a)
    @assert alength >= 2 "type $a not long enough for reduce (needs at least two elements)"
    quote 
        s = f(a[1], a[2])
        $([:(s = f(s, a[$i])) for i=3:alength]...)
    end
end

# either gives an expression to index into the fixedsizearray, or just inserts the constant
accessor{T <: AbstractFixedArray}(arg::Type{T}, i::Integer, name::Symbol) = :($name[$i])
accessor{T                      }(arg::Type{T}, i::Integer, name::Symbol) = :($name)

name(typ::DataType) = typ.name.name 

#To simplify things, map_expression assumes that the argument names of the staged functions are f for the callable and a,b,c,.. for the other arguments.
function map_expression{F <: Func, FSA <: AbstractFixedArray}(f::Type{F}, fsa::Type{FSA}, args...)
    FN = first(super(F).parameters)
    @assert length(args) == FN "cardinality of callable doesn't match arguments. Callable: $FN, args: $(length(args))"
    argnames = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j]
    #call f for every index in FSA, with all the args.
    expr = ntuple(i -> :(
        f( $(ntuple(j-> accessor(args[j], i, argnames[j]), FN)...)) 
    ),length(FSA))
    # the type FSA from map can't be used, because it contains the parameter already (e.g. FSA{Int}), which results in a conversion.
    # map should return a fixedsize array with the return type of f, so the parameter has to be stripped of FSA.
    type_name = :(Main.$(symbol(string(fsa.name.name)))) # custom types are not known in the Module FixedSizeArray, but in Main
    eval(type_name)
    :($type_name($(expr...)))
end


stagedfunction map{FSA <: AbstractFixedArray}(f::Func{1}, a::FSA)
    map_expression(f, a, a)
end
stagedfunction map{FSA <: AbstractFixedArray}(f::Func{2}, a::FSA, b::FSA)
    map_expression(f, a, a, b)
end
stagedfunction map{FSA <: AbstractFixedArray, REAL <: Real}(f::Func{2}, a::FSA, b::REAL)
    map_expression(f, a, a, b)
end
stagedfunction map{FSA <: AbstractFixedArray, REAL <: Real}(f::Func{2}, a::REAL, b::FSA)
    map_expression(f, b, a, b)
end



stagedfunction row{T, M, N}(A::AbstractFixedMatrix{T, M, N}, i)
  return :(AbstractFixedVector{T, M}($(ntuple(j->:(A[i,$j]), M)...)))
end
stagedfunction col{T, M, N}(A::AbstractFixedMatrix{T, M, N}, i)
  return :(AbstractFixedVector{T, N}($(ntuple(j->:(A[$j, i]), N)...)))
end
# Matrix
stagedfunction (*){T, M, N, K}(a::AbstractFixedMatrix{T, M, N}, b::AbstractFixedMatrix{T, N, K})
    :(AbstractFixedMatrix{$T, $M, $K}( 
         $([foldl((v0, prod_expr) -> :($v0 + $prod_expr), [:(a.($(sub2ind((M,N),i,k)))*b.($(sub2ind((N,K),k,j)))) for k=1:N]) for i=1:M, j=1:K]...)
    ))
end