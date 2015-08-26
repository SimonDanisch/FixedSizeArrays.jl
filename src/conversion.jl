
function convert{DA <: Array, FSA <: FixedArray}(::Type{DA}, b::FSA)
    result = Array(eltype(b), size(b)...)
    for i=1:length(b)
        result[i] = b[i]
    end
    result
end
function convert{T, N, FSA <: FixedArray}(::Type{Array{T, N}}, b::FSA)
    result = Array(eltype(b), size(b)...)
    for i=1:length(b)
        result[i] = T(b[i])
    end
    result
end

convert{R, C, R2, C2, T}(a::Type{Mat{R, C, T}}, b::Mat{R2, C2, T}) =
    Mat(ntuple(c -> ntuple(r -> b[r,c], Val{R}), Val{C}))
