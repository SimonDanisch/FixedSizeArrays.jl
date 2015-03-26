stagedfunction convert{FSA <: FixedArray, T, ND}(::Type{FSA}, v::Array{T, ND})
    conversion = isleaftype(FSA) ? quote # differentiate between FixedArrays which come with an element type and abstract one without
        ptr_typ, converted = Ptr{FSA}, convert(Array{eltype(FSA)}, v)  
    end : quote
        ptr_typ, converted = Ptr{FSA{T}}, v
    end
    quote
        length(v) != $(length(FSA)) && throw(DimensionMismatch("Lenght of Array: $(length(v)) is not equal to length of FixedSizeArray: $(length(FSA))"))
        $conversion
        unsafe_load(ptr_typ(pointer(converted)))
    end
end

function convert{DA <: DenseArray, FSA <: FixedArray}(::Type{DA}, b::FSA)
    pointer_to_array(Ptr{eltype(FSA)}(pointer_from_objref(b)), size(FSA))
end
function convert{T,N, FSA <: FixedArray}(t::Type{Array{T, N}}, b::FSA)
    convert(t, convert(Array, b))
end

function convert{T1 <: FixedArray, T2 <: FixedArray}(a::Type{T1}, b::Array{T2})
    @assert sizeof(b) % sizeof(a) == 0 "Type $a, with size: ($(sizeof(a))) doesn't fit into the array b: $(length(b)) x $(sizeof(eltype(b)))"
    reinterpret(a, b, (div(sizeof(b), sizeof(a)),))
end

function convert{FSAA <: FixedArray}(a::Type{FSAA}, b::FixedArray)
	typeof(b) == a && return b
    map(IndexFunctor(b), a)
end
