tuple_to_string(t::(), sep)            = ""
tuple_to_string(t::(Any,), sep)        = "$(t[1])"
tuple_to_string(t::(Any, Any...), sep) = "$(t[1])$sep" * tuple_to_string(Base.tail(t), sep)
vec_name(sz::(Integer, Integer...))    = symbol("NVec" * tuple_to_string(sz, 'x'))
vec_name(sz::(Integer,))               = symbol("NVec" * string(first(sz)))

gen_fixedsizevector_type(name::DataType, T::Symbol, N::Int) = gen_fixedsizevector_type(symbol(string(name.name.name)), T, N)

# maps the dimension of vectors always to (length,)
flatten_dimension(SIZE) = SIZE
flatten_dimension(SIZE::(Int, Int)) = any(x->x==1, SIZE) ? prod(SIZE) : SIZE

function gen_fixedsizevector_type(SIZE::(Integer...))
    # Make sure, that vectors (1,x) OR (x,1) have the same type => (x,)
    SIZE 		= flatten_dimension(SIZE)
    fields      = [Expr(:(::), fieldname(i), :T) for i=1:prod(SIZE)]
    NDim        = length(SIZE)
    typename    = vec_name(SIZE)
    #only eval if not already defined
    !isdefined(Main, typename) && eval(Main, quote
        immutable $(typename){T} <: FixedArray{T, $NDim, $SIZE}
            $(fields...)
        end
    end)
    :(Main.$typename)
end


call{T, NDim, SIZE}(t::Type{FixedArray{T, NDim, SIZE}}, data::T...) = t(data)
stagedfunction call{T, NDim, SIZE}(t::Type{FixedArray{T, NDim, SIZE}}, data)
    N = length(data)
    @assert prod(SIZE) == N "not the right dimension"
    typename = gen_fixedsizevector_type(SIZE)
    :($(typename)(data...))
end
call{FS <: FixedArray, T, N}(::Type{FS}, a::Array{T, N}) = FixedArray{T, N, size(a)}(a...) 
nvec{T, N}(x::Array{T,N})             = FixedArray(x)
nvec{T}(x::T...)                      = FixedArray{T, 1, (length(x),)}(x)
nvec{T}(SIZE::(Integer...,), x::T...) = FixedArray{T, length(SIZE), SIZE}(x)


