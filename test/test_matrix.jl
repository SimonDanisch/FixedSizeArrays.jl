
abstract AbstractFixedArray{T, N, SZ}
abstract AbstractFixedArrayWrapper{T <: AbstractFixedArray}
immutable Vec{T, N, FS} <: AbstractFixedArrayWrapper{FS}
    val::FS
end