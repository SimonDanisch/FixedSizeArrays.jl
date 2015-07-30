@inline map{FSA <: FixedArray}(f::Func{1}, a::FSA)          = FSA(map(f, a.(1)))
@inline map{FSA <: FixedArray}(f::Func{2}, a::FSA, b::FSA)  = FSA(map(f, a.(1), b.(1)))

inner_map(len::Int, inner) = ntuple(i -> quote f($(inner(i))) end, len)

function reduce{FSA <: FixedArray}(f::Func{2}, a::FSA)
    red = f(a[1], a[2])
    for i=3:length(a)
        red = f(red, a[i])
    end
    red
end

@generated function map{FSA <: FixedMatrix}(f::Func{1}, a::Union(FSA, Type{FSA}))
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:R
        push!(exprs, quote tuple($(inner_map(C, x->:($(sub2ind((R,C), i,x))))...)) end)
    end
    quote 
        Main.Mat{$R, $C, $T}(tuple($(exprs...)))
    end
end

@generated function map{FSA <: FixedMatrix}(f::Func{2}, a::FSA, b::Number)
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:R
        push!(exprs, quote tuple($(ntuple(j -> quote f(a.(1)[$i][$j], b) end, C)...)) end)
    end
    quote 
        Main.Mat{$R, $C, $T}(tuple($(exprs...)))
    end
end
@generated function map{FSA <: FixedMatrix}(f::Func{2}, a::Number, b::FSA)
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:R
        push!(exprs, quote tuple($(ntuple(j -> quote f(a, b.(1)[$i][$j]) end, C)...)) end)
    end
    quote 
        Main.Mat{$R, $C, $T}(tuple($(exprs...)))
    end
end
@generated function map{FSA <: FixedArray}(f::Func{2}, a::FSA, b::Number)
    exprs = ntuple(i -> :(f(a[$i], b)), length(a))
    :(FSA($(exprs...)))
end
@generated function map{FSA <: FixedArray}(f::Func{2}, a::Number, b::FSA)
    exprs = ntuple(i -> :(f(a, b[$i])), length(b))
    :(FSA($(exprs...)))
end


@generated function map{FSA <: FixedArray}(f::Func{1}, a::Type{FSA})
    :(FSA($([:(f($i)) for i=1:length(FSA)]...)))
end

@generated function map{FSA <: FixedArray, F <: Func{1}}(f::Union(Type{F}, F), a::Type{FSA})
    :(FSA($([:(f($i)) for i=1:length(FSA)]...)))
end
@generated function map{FSA <: FixedArray, F <: Func{2}}(f::Union(Type{F}, F), a::Type{FSA})
    :(FSA($([:(f($i, $j)) for i=1:size(FSA,1), j=1:size(FSA,2)]...)))
end

