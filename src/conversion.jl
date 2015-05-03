# dealing with an abstract FixedArrayType, this code is terrible
@generated function convert{SZ, T, ND}(a::Type{FixedArray{T, ND, SZ}}, v::Array{T, ND})
    returntype = gen_fixedsizevector_type(SZ, false)
    quote
        length(v) != length(a) && throw(DimensionMismatch("Lenght of Array: $(length(v)) is not equal to length of FixedSizeArray: $(length(FSA))"))
        unsafe_load(Ptr{$returntype{T}}(pointer(v)))
    end
end
#did I say terrible?
@generated function convert{FSA <: FixedArray, T, ND}(::Type{FSA}, v::Array{T, ND})
    symbol(FSA.name.name) == :FixedArray && return :(convert(FixedArray{T, ND, size(v)}, v))
    conversion = isleaftype(FSA) ? quote # differentiate between FixedArrays which come with an element type and abstract one without
        ptr_typ, converted, len = Ptr{FSA}, convert(Array{eltype(FSA)}, v), length(FSA)
    end : quote
        ptr_typ, converted, len = Ptr{FSA{T}}, v, length(v)
    end
    quote
        $conversion
        length(v) != len && throw(DimensionMismatch("Lenght of Array: $(length(v)) is not equal to length of FixedSizeArray: $(length(FSA))"))
        unsafe_load(ptr_typ(pointer(converted)))
    end
end

function convert{DA <: DenseArray, FSA <: FixedArray}(::Type{DA}, b::FSA)
    elt = eltype(FSA)
    ptr = pointer_from_objref(b)
    sz = size(FSA)
    pointer_to_array(Ptr{elt}(ptr), sz)
end
function convert{T,N, FSA <: FixedArray}(t::Type{Array{T, N}}, b::FSA)
    convert(t, convert(Array, b))
end

function convert{T1 <: FixedArray, T2 <: FixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(b) % sizeof(a) == 0 "Type $a, with size: ($(sizeof(a))) doesn't fit into the array b: $(length(b)) x $(sizeof(eltype(b)))"
    eltype(T1) != eltype(T2) && return map(T1, b)
    reinterpret(T1, b, (div(sizeof(b), sizeof(a)),))
end

function convert{FSAA <: FixedArray}(a::Type{FSAA}, b::FixedArray)
	typeof(b) == a && return b
    map(IndexFunctor(b), a)
end
