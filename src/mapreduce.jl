map{FSA <: FixedVector}(f::Func{1}, a::FSA)         = FSA(map(f, a.(1)))
map{FSA <: FixedVector}(f::Func{2}, a::FSA, b::FSA) = FSA(map(f, a.(1), b.(1)))

@generated function map{FSA <: FixedVectorNoTuple}(f::Func{2}, a::FSA, b::FSA)
    :(FSA($(ntuple( i->:(f(a.($i), b.($i))), length(a))...)))
end

map{R, C, T}(f::Func{1}, a::Mat{R,C,T}) = Mat(ntuple(c->map(f, a.(1)[c]), Val{C}))

map{R, C, T, T2}(f::Func{2}, a::Mat{R,C,T}, b::Mat{R,C,T2}) = Mat(ntuple(c->map(f, a.(1)[c], b.(1)[c]), Val{C}))


inner_map(len::Int, inner) = ntuple(i -> quote f($(inner(i))) end, len)

function reduce{FSA <: FixedArray}(f::Func{2}, a::FSA)
    red = f(a[1], a[2])
    @inbounds for i=3:length(a)
        red = f(red, a[i])
    end
    red
end
function Base.reduce{R,C,T}(f::Base.Func{2}, a::Mat{R,C,T})
    red = reduce(f, a.(1)[1])
    @inbounds for i=2:C
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



@generated function map{R, C, T}(f::Func{1}, a::Mat{R, C, T})
    exprs = [:(map(f, a.(1)[$i])) for i=1:C]
    :(Mat(tuple($(exprs...))))
end
@generated function map{FSA <: FixedMatrix}(f::Func{1}, a::Type{FSA})
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:C
        push!(exprs, quote tuple($(inner_map(R, x->:($(sub2ind((R,C), i,x))))...)) end)
    end
    quote
        Mat(tuple($(exprs...)))
    end
end

@generated function map{FSA <: FixedMatrix}(f::Func{2}, a::FSA, b::Number)
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:C
        push!(exprs, quote tuple($(ntuple(j -> quote f(a.(1)[$i][$j], b) end, R)...)) end)
    end
    quote
        Mat{$R, $C, $T}(tuple($(exprs...)))
    end
end
@generated function map{FSA <: FixedMatrix}(f::Func{2}, a::Number, b::FSA)
    exprs = []
    R, C = size(FSA)
    T = eltype(FSA)
    for i=1:C
        push!(exprs, quote tuple($(ntuple(j -> quote f(a, b.(1)[$i][$j]) end, R)...)) end)
    end
    quote
        Mat{$R, $C, $T}(tuple($(exprs...)))
    end
end
@generated function map{FSA <: FixedArray}(f::Func{2}, a::FSA, b::Number)
    exprs = ntuple(i -> :(f(a[$i], b)), length(a))
    :($FSA($(exprs...)))
end
@generated function map{FSA <: FixedArray}(f::Func{2}, a::Number, b::FSA)
    exprs = ntuple(i -> :(f(a, b[$i])), length(b))
    :($FSA(tuple($(exprs...))))
end


@generated function map{FSA <: FixedArray}(f::Func{1}, a::Type{FSA})
    :($FSA(tuple($([:(f($i)) for i=1:length(FSA)]...))))
end

@generated function map{FSA <: FixedArray, F <: Func{1}}(f::Union(Type{F}, F), a::Type{FSA})
    :($FSAa(tuple($([:(f($i)) for i=1:length(FSA)]...))))
end
@generated function map{FSA <: FixedArray, F <: Func{2}}(f::Union(Type{F}, F), a::Type{FSA})
    :($FSA(tuple($([:(f($i, $j)) for i=1:size(FSA,1), j=1:size(FSA,2)]...))))
end
