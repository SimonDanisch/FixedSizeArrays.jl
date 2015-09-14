function convert{DA <: Array, FSA <: FixedArray}(::Type{DA}, b::FSA)
    result = Array(eltype(b), size(b)...)
    @inbounds for i=1:length(b)
        result[i] = b[i]
    end
    result
end

convert{R, C, R2, C2, T}(a::Type{Mat{R, C, T}}, b::Mat{R2, C2, T}) =
    Mat(ntuple(c -> ntuple(r -> b[r,c], Val{R}), Val{C}))

#conversion
convert{T <: Tuple}(::Type{T}, x::Real)                   = ntuple(ConstFunctor(eltype(T)(x)), Val{length(T.parameters),})
convert{T <: Tuple, FSA <: FixedArray}(::Type{T}, x::FSA) = map(eltype(T), x.(1)[1:length(T.parameters)])
convert{T <: FixedArray}(t::Type{T}, f::T)                = f
convert{FSA1 <: FixedArray}(t::Type{FSA1}, f::FixedArray) =
    map(ConversionIndexFunctor(f, eltype_or(FSA1, eltype(typeof(f)))), FSA1)

convert{T <: FixedArray}(::Type{T}, x) = T(x)
convert{T <: FixedArray}(::Type{T}, x...) = T(x)
