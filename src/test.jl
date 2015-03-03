immutable Vec3{T} <: DenseArray{T, 1}
	x::T
	y::T
	z::T
end
Vec3(a) 												= Vec3(a,a,a)
Base.zero{T}(a::Type{Vec3{T}}) 							= Vec3(zero(T));
Base.length(a::Vec3) 							  		= 3;
Base.getindex(A::Vec3, i::Integer) 				  		= A.(i);
Base.size(A::Vec3)						  		  		= (3,);
Base.size(A::Vec3, i::Integer)			  		  		= (3,)[i];
Base.ndims(A::Vec3)						  		  		= 1;
Base.eltype{T}(A::Vec3{T})						  		= T;

Base.start(A::Vec3)            					  		= 1;
Base.next (A::Vec3, state::Integer) 			  		= (A[state], state+1);
Base.done (A::Vec3, state::Integer) 		      		= length(A) < state;

#with cheated similar
Base.similar{T}(a::Vec3, t::Type{T}, z::(Integer...)) 	= Array(t,z)

#without similiar
try
	println("[Vec3(1,3,4); Vec3(1,2,3)]: ", [Vec3(1,3,4); Vec3(1,2,3)]) # #MethodError(similar,([1,3,4],Int64,(6,)))
catch e
	println("ERROR: [Vec3(1,3,4); Vec3(1,2,3)]")
	println(e)
end
println("###############################")
try
println("[Vec3(1,3,4), Vec3(1,2,3)]: ", [Vec3(1,3,4), Vec3(1,2,3)]) #[a,b] concatenation is deprecated; use [a;b] instead + MethodError(similar,([1,3,4],Int64,(6,)))
catch e
	println("ERROR: [Vec3(1,3,4), Vec3(1,2,3)]")
	println(e)
end
println("###############################")

try
@show [Vec3(1,3,4)  Vec3(1,2,3)] #MethodError(similar,([1,3,4],Int64,(3,2)))
catch e
	println("[Vec3(1,3,4)  Vec3(1,2,3)]")
	println(e)
end
println("###############################")

try

@show a = Vec3[Vec3(1,3,4); Vec3(1,2,3)] #MethodError(similar,([1,3,4],Vec3{T},(6,)))
catch e
	println("Vec3[Vec3(1,3,4); Vec3(1,2,3)]")
	println(e)
end
println("###############################")

try
@show a = Vec3[Vec3(1,3,4), Vec3(1,2,3)] #Vec3[[1,3,4],[1,2,3]]
catch e
	println("Vec3[Vec3(1,3,4), Vec3(1,2,3)]")
	println(e)
end

println("###############################")

try
@show a = Vec3[Vec3(1,3,4) Vec3(1,2,3)] #MethodError(similar,([1,3,4],Vec3{T},(3,2)))
catch e
	println("Vec3[Vec3(1,3,4) Vec3(1,2,3)]")
	println(e)
end
#=



@show a = [Vec3(1,3,4); Vec3(1,2,3)] # ERROR: no setindex! (obviously, as this is the whole point of similar)
@show a = [Vec3(1,3,4), Vec3(1,2,3)] # ERROR: no setindex!
@show a = [Vec3(1,3,4)  Vec3(1,2,3)] # ERROR: BoundsError at abstractarray.jl:562

@show a = Vec3[Vec3(1,3,4); Vec3(1,2,3)] # ERROR: `zero` has no method matching zero(::Type{Vec3{T}})
@show a = Vec3[Vec3(1,3,4), Vec3(1,2,3)] # Array{Vec3{T},1} 
@show a = Vec3[Vec3(1,3,4)  Vec3(1,2,3)] # ERROR:`zero` has no method matching zero(::Type{Vec3{T}})


#with similar, but immutable



@show a = [Vec3(1,3,4); Vec3(1,2,3)] # -> [1,3,4,1,2,3]
@show a = [Vec3(1,3,4), Vec3(1,2,3)] # -> [1,3,4,1,2,3] + deprecation warning
@show a = [Vec3(1,3,4)  Vec3(1,2,3)] # [1 1
 									 #  3 2
 									 #  4 3]

@show a = Vec3[Vec3(1,3,4); Vec3(1,2,3)] # has no method matching convert(::Type{Vec3{T}}, ::Int64)
@show a = Vec3[Vec3(1,3,4), Vec3(1,2,3)] # Vec3[[1,3,4],[1,2,3]]
@show a = Vec3[Vec3(1,3,4)  Vec3(1,2,3)] # error has no method matching convert(::Type{Vec3{T}}, ::Int64)

stagedfunction map{T <: FixedArray}(f, A::Type{T})
	quote
		A($([:(f($i)) for i=1:length(T)]...))
	end
end



stagedfunction products(f, it...)
	variables 		= ntuple(fieldname, length(it)) #-> [:i_1, :i_2, ...]
	iterator_access = Expr(:block, [:($(variables[i]) = it[$i])	for i=1:length(it)]...) # -> for i_1 in it[1], i_2 in it[2]...
	for_expression  = Expr(:for, 
		iterator_access, 
		:(f($(variables...))) # -> body of the nested for loop f(i_1, i_2, ...)
	)
	quote
	$for_expression
	nothing
	end
end

immutable IndexFunctor2{T}
	A::T
	inds::T
end
call(f::IndexFunctor2, i) = f.A[inds[i]]


@show map(IndexFunctor2(rand(5,5), ), FixedArray{Float64, 2, (5,5)})
=#