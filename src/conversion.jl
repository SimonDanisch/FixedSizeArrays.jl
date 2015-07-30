
function convert{DA <: DenseArray, FSA <: FixedArray}(::Type{DA}, b::FSA)
    elt = eltype(FSA)
    ptr = pointer_from_objref(b)
    sz = size(FSA)
    pointer_to_array(Ptr{elt}(ptr), sz)
end
convert{T,N, FSA <: FixedArray}(t::Type{Array{T, N}}, b::FSA) =
    convert(t, convert(Array, b))

function convert{T1 <: FixedArray, T2 <: FixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(b) % sizeof(a) == 0 "Type $a, with size: ($(sizeof(a))) doesn't fit into the array b: $(length(b)) x $(sizeof(eltype(b)))"
    eltype(T1) != eltype(T2) && return map(T1, b)
    reinterpret(T1, b, (div(sizeof(b), sizeof(a)),))
end

function convert{FSAA <: FixedArray}(a::Type{FSAA}, b::FixedArray)
	typeof(b) == a && return b
    map(IndexFunctor(b), a)
end
