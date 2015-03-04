tic()
using FixedSizeArrays
tuple_to_string(t::(), sep)            			= ""
tuple_to_string(t::(Any,), sep)        			= "$(t[1])"
tuple_to_string(t::(Any, Any...), sep) 			= "$(t[1])$sep" * tuple_to_string(Base.tail(t), sep)
vec_name(sz::(Integer, Integer...), base_name) 	= symbol(base_name * tuple_to_string(sz, 'x'))
vec_name(sz::(Integer,), base_name)            	= symbol(base_name * string(first(sz)))

function fixedarray_type_expr(base_name, fields, SIZE::(Integer...))
	@assert prod(SIZE) <= length(fields)
    fields      = [Expr(:(::), fields[i], :T) for i=1:prod(SIZE)]
    NDim        = length(SIZE)
    typename    = vec_name(SIZE, base_name)
    quote
        immutable $(typename){T} <: FixedArray{T, $NDim, $SIZE}
            $(fields...)
        end
    end
end
function genfsa(name, fields, N)
	expr = [fixedarray_type_expr(name, fields, (i,)) for i in N]
	eval(Expr(:block, expr...))
end
genfsa("Vector", 	[:x,:y,:z,:w], 				1:4)
genfsa("Point",  	[:x,:y,:z,:w], 				1:4)
genfsa("Normal",  	[:x,:y,:z,:w], 				1:4)
genfsa("UV",		[:u,:v], 					2:2)
genfsa("UVW",  		[:u,:v,:w], 				3:3)
genfsa("Face",   	[:a,:b,:c,:d,:e,:f,:g,:h], 	3:8)
toc()
@show Vector4(1,2,3,4) + Vector4(1,2,3,4)

