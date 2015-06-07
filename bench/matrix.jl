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

const a = rand(10^7)
const b = rand(10^7)


@time c = parallelmap(+, a, b)
@time parallelmap(+, a, b)
@time d = parallelmap(Base.AddFun(), a, b)
@time parallelmap(Base.AddFun(), a, b)


@time e = simplemap(+, a, b)
@time simplemap(+, a, b)
@time f = simplemap(Base.AddFun(), a, b)
@time simplemap(Base.AddFun(), a, b)

@time g = map(+, a, b)
@time map(+, a, b)
@time h = map(Base.AddFun(), a, b)
@time map(Base.AddFun(), a, b)

@assert c==d && c==e && c==f && c==g && c==h


