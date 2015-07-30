using Base.Test
importall Base
abstract FixedArray{ND, SZ, T}

immutable Mat{R, C, L, T} <: FixedArray{2, Tuple{R, C}, T}
	_::NTuple{L, T}
end
@inline getindex{R, C, L, T}(a::Mat{R, C, L, T}, i, j) = a.(1)[((j-1)*R+(i-1))+1]
@inline getindex{N, M, L, T}(a::Mat{N, M, L, T}, i)    = a.(1)[i]
@inline Base.length{R, C, L, T}(::Mat{R, C, L, T}) 	   = L


immutable IndexFunc2{T} <: Base.Func{1}
	arg::T
	i::Int
end
@inline Base.call{T}(a::IndexFunc2{T}, j) = a.arg[j, a.i]
@inline column{R, C, L, T}(a::Mat{R, C, L, T}, i) = ntuple(IndexFunc2(a,i), Val{C})::NTuple{C, T}
immutable IndexFunc{T} <: Base.Func{1}
	arg::T
	i::Int
end
Base.call{T}(a::IndexFunc{T}, j) = a.arg[a.i,j] 
@inline row{R, C, L, T}(a::Mat{R, C, L, T}, j) = ntuple(IndexFunc(a,j), Val{R})::NTuple{R, T}


immutable Mat2{R, C, T} <: FixedArray{2, Tuple{R, C}, T}
	_::NTuple{C, NTuple{R, T}}
end
@inline getindex{N, M, T}(a::Mat2{N, M, T}, i, j) = a.(1)[j][i]
@inline getindex{N, M, T}(a::Mat2{N, M, T}, i) 	  = a[ind2sub((N,M), i)...]

@inline column{R, C, T}(a::Mat2{R, C, T}, i) = a.(1)[i]
@inline row{R, C, T}(a::Mat2{R, C, T}, j) 	 = ntuple(IndexFunc(a,j), Val{R})::NTuple{R, T}
@inline Base.length{R, C, T}(::Mat2{R, C, T}) = R*C



@generated function (*){T, M, N, K}(a::Mat2{M, N, T}, b::Mat2{N, K, T})
    expr = []
    for i=1:M 
        rowt = [:(+($(ntuple(k->:(a.(1)[$k][$i]*b.(1)[$k][$j]), N)...))) for j=1:K]
        push!(expr, :(tuple($(rowt...))))
    end
    :(Mat2(tuple($(expr...))))
end
@generated function (*){T, M, N, K, Len}(a::Mat{M, N, Len, T}, b::Mat{N, K, Len, T})
    :(Mat{$M, $K, $(M*K), $T}(tuple(
         $([:(+($(ntuple(k->:(

         	a.(1)[$(sub2ind((M,N), k, i))]*b.(1)[$(sub2ind((N,K), k, j))]

         ), N)...)))for i=1:M, j=1:K]...)
    )))
end

function Base.reduce{R,C,T}(f::Base.Func{2}, a::Mat2{R,C,T})
    red = reduce(f, a.(1)[1])
    for i=2:C
        red = f(red, reduce(f, a.(1)[i]))
    end
    red
end
function Base.reduce{FSA <: FixedArray}(f::Base.Func{2}, a::FSA)
    red = f(a[1], a[2])
    for i=3:length(a)
        red = f(red, a[i])
    end
    red
end
sum{FSA <: FixedArray}(a::FSA) = reduce(Base.AddFun(), a)

.*(a::NTuple, b::NTuple) = map(Base.MulFun(), a, b)
dot{T}(a::NTuple{4,T}, b::NTuple{4,T}) = sum(a.*b)

function test(n, a, b)
	r = a*b
	for i=1:n
		r = r*a*b
	end
	r
end

function test_row(n, a)
	r = row(a,1)
	for i=1:n
		for j=1:4
			r = row(a,j)
		end
	end
	r
end

function test_column(n, a)
	r = column(a,1)
	for i=1:n
		for j=1:4
			r = column(a,j)
		end
	end
	r
end
function test_sum(n, a)

	r = sum(a)
	for i=1:n
		for j=1:4
			r = sum(a)+j
		end
	end
	r
end

function test_dot(n, a)
	r = dot(row(a,1), column(a,1))
	for i=1:n
		for j=1:4
			r += dot(row(a,j), column(a,5-j))
		end
	end
	r
end
function test()
	aj = [1 1 1 1; 2 2 2 2; 3 3 3 3; 4 4 4 4]

	const a = Mat{4, 4, 16, Int}((1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
	const b = Mat{4, 4, 16, Int}((1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))

	for i=1:length(aj)
		@test aj[i] == a[i]
	end
	for i=1:4, j=1:4
		@test aj[i,j] == a[i,j]
	end

	@time test(10^6, a, b)
	@time test(10^6, a, b)

	@time test_column(10^6, a)
	@time test_column(10^6, a)

	@time test_row(10^6, a)
	@time test_row(10^6, a)

	sum1 = test_sum(10^6, a)
	@time test_sum(10^6, a)
	@time test_sum(10^6, a)
	println("dot:")
	@time test_dot(10^6, a)
	@time test_dot(10^6, a)

	const a2 = Mat2((
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4)
	))
	const b2 = Mat2((
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4),
		(1,2,3,4)
	))
	for i=1:length(aj)
		@test aj[i] == a2[i]
	end
	for i=1:4, j=1:4
		@test aj[i,j] == a2[i,j]
	end
	row(a2, 1)

	@time test(10^6, a2, b2)
	@time test(10^6, a2, b2)


	@time test_column(10^6, a2)
	@time test_column(10^6, a2)

	@time test_row(10^6, a2)
	@time test_row(10^6, a2)

	sum2 = test_sum(10^6, a2)
	@time test_sum(10^6, a2)
	@time test_sum(10^6, a2)

	@test sum2 == sum1

	println("dot:")
	@time test_dot(10^6, a2)
	@time test_dot(10^6, a2)
end
@inbounds test()