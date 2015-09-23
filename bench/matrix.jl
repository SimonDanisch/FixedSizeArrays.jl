function parallelmap{S,T}(f, A::AbstractArray{S}, B::AbstractArray{T})
    F = similar(A, promote_type(S,T), promote_shape(size(A),size(B)))
    @parallel for i in eachindex(A,B)
        @inbounds F[i] = f(A[i], B[i])
    end
    return F
end

function simplemap{S,T}(f, A::AbstractArray{S}, B::AbstractArray{T})
    F = similar(A, promote_type(S,T), promote_shape(size(A),size(B)))
    for i in eachindex(A,B)
        @inbounds F[i] = f(A[i], B[i])
    end
    return F
end
function test()
    const a = rand(10^7)
    const b = rand(10^7)


    c = parallelmap(+, a, b)
    println("parallel:")
    @time parallelmap(+, a, b)
    d = parallelmap(Base.AddFun(), a, b)
    println("functor:")
    @time parallelmap(Base.AddFun(), a, b)


    e = simplemap(+, a, b)

    println("simplemap:")
    @time simplemap(+, a, b)
    f = simplemap(Base.AddFun(), a, b)
    println("functor:")
    @time simplemap(Base.AddFun(), a, b)

    g = map(+, a, b)
    println("map:")
    @time map(+, a, b)
    h = map(Base.AddFun(), a, b)
    println("functor:")
    @time map(Base.AddFun(), a, b)

    println("plus")
    i = a+b
    @time a+b

    @assert c==d && c==e && c==f && c==g && c==h && c==i
end
test()

