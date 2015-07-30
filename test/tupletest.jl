importall Base
import Base.Func
#=
typealias FixedMatrix{T, M, N} NTuple{M, NTuple{N, T}}
=#
immutable FixedMatrix{T, M, N}
    x::NTuple{M, NTuple{N, T}}
end
typealias FixedVector{N, T} NTuple{N,T}
@inline Base.getindex{T, M, N}(x::FixedMatrix{T, M, N}, i::Integer, j::Union(Integer, Range)) = x.x[i][j]

@inline row(x, i::Integer) = x[i]

@generated function column{T, Row, Column}(A::FixedMatrix{T, Row, Column}, j::Integer)
    fields = ntuple(i-> :(A[$i][j]), Row)
    :(tuple($(fields...)))
end
@inline map_extended(F, a::NTuple, b::NTuple...) = map(F, a, b...)


@generated function map_extended(f::Func{2}, a::NTuple, b::Number)
    exprs = ntuple(i -> :(f(a[$i], b)), length(a.parameters))
    :(tuple($(exprs...)))
end

@generated function map_extended(f::Func{2}, a::Number, b::NTuple)
    exprs = ntuple(i -> :(f(a, b[$i])), length(b.parameters))
    :(tuple($(exprs...)))
end
@generated function map_extended{FSA <: NTuple}(f::Func{1}, a::Type{FSA})
    :(tuple($([:(f($i)) for i=1:length(FSA.parameters)]...)))
end

@inline .*{N}(a::NTuple{N}, b::NTuple{N}) = map(Base.MulFun(), a, b)

@inline sum{T, N}(a::NTuple{T, N}) = reduce(Base.AddFun(), a)
@inline dot{T, N}(a::NTuple{T, N}, b::NTuple{T, N}) = sum(a.*b)


@generated function (*){T, M, N, K}(a::FixedMatrix{T, M, N}, b::FixedMatrix{T, N, K})
	expr = []
	for i=1:M 
		rowt = [:(+($(ntuple(k->:(a.x[$k][$i]*b.x[$k][$j]), N)...))) for j=1:K]
		push!(expr, :(tuple($(rowt...))))
	end
	:(FixedMatrix(tuple($(expr...))))
end

function test(a)
    r = a
    for i=1:10^6
        r = r*a
    end
    r
end

function test1()
    const a = FixedMatrix((
        (1,2,3,4),
        (1,2,3,4),
        (1,2,3,4),
        (1,2,3,4),
    ))

    const b = [
        1 2 3 4;
        1 2 3 4;
        1 2 3 4;
        1 2 3 4;
    ]
    @time test(a)
    @time test(a)

    @time test(b)
    @time test(b)
end

test1()
