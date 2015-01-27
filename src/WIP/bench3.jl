immutable A{T, SZ}
    ptr::Ptr{T}
end
immutable Vec4{T}
    x::T
    y::T
    z::T
    w::T
end
function A4{T}(a::T,b::T,c::T,d::T)
    ptr = convert(Ptr{T}, c_calloc(4, sizeof(T)))
    unsafe_store!(ptr, a, 1)
    unsafe_store!(ptr, b, 2)
    unsafe_store!(ptr, c, 3)
    unsafe_store!(ptr, d, 4)
    A{T, (4,)}(ptr)
end 
Base.length{T, SZ}(a::A{T, SZ}) = prod(SZ)
immutable Field{S} end

function Base.getindex{T, SZ}(a::A{T, SZ}, i::Field{:x})
    unsafe_load(a.ptr, 1)
end
function Base.getindex{T, SZ}(a::A{T, SZ}, i::Field{:y})
    unsafe_load(a.ptr, 2)
end
function Base.getindex{T, SZ}(a::A{T, SZ}, i::Field{:z})
    unsafe_load(a.ptr, 3)
end
function Base.getindex{T, SZ}(a::A{T, SZ}, i::Field{:w})
    unsafe_load(a.ptr, 4)
end


function Base.getindex{T}(a::Vec4{T}, i::Integer)
    getfield(a, i)
end
const a = A4(1.0,2.0,3.0,4.0)
const b = Vec4(1.0,2.0,3.0,4.0)


@inline function +(a::A, b::A)
    A4(a[Field{:x}()]+b[Field{:x}()],
    a[Field{:y}()]+b[Field{:y}()],
    a[Field{:z}()]+b[Field{:z}()],
    a[Field{:w}()]+b[Field{:w}()],) 
end
function +(a::Vec4, b::Vec4)
    Vec4(a.x+b.x, 
        a.x+b.x,
        a.x+b.x,
        a.x+b.x)
end

@inline function Base.map{T, SZ}(a::A{T, SZ}, i)
    [1,2,3,4]
end
@inline function Base.map(a::Vec4, i)
    Vec4(2,2,3,4)
end
function test(n, r)
    gc_disable()
    for i=1:n
        map(r,1)
    end
    gc_enable()

end
println(code_lowered(c_calloc, (Int, Int)))
test(1, a)
@time test(10^7, a)
@time test(10^7, a)
@time test(10^7, a)
@time test(10^7, a)
test(1,b)
@time test(10^7, b)
@time test(10^7, b)
@time test(10^7, b)
@time test(10^7, b)