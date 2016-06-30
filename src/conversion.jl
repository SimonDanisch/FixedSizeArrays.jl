function convert{DA <: Array, FSA <: FixedArray}(::Type{DA}, b::FSA)
    result = Array(eltype(b), size(b)...)
    @inbounds for i=1:length(b)
        result[i] = b[i]
    end
    result
end


convert{FA<:FixedArray}(a::Type{FA}, b::FA) = b
convert{FA<:Mat}(a::Type{FA}, b::FA) = b

function convert{R, C, R2, C2, T}(MT::Type{Mat{R, C, T}}, b::Mat{R2, C2, T})
    MT(ntuple(c -> ntuple(r -> b[r,c], Val{R}), Val{C}))
end

#conversion
convert{N}(T::Type{NTuple{N}}, x::Real) = ntuple(ConstFunctor(eltype(T)(x)), Val{N})
convert{T <: Tuple, FSA <: FixedArray}(::Type{T}, x::FSA) = convert(T, Tuple(x))

function convert{FSA1 <: FixedArray}(t::Type{FSA1}, f::FixedArray)
    map(ConversionIndexFunctor(f, eltype_or(FSA1, eltype(typeof(f)))), FSA1)
end

convert{T <: FixedArray}(::Type{T}, x) = T(x)
convert{T <: FixedArray}(::Type{T}, x...) = T(x)
