abstract FixedSizeVector{T, N}

gen_fixedsizevector_type(name::DataType, T::Symbol, N::Int) = gen_fixedsizevector_type(symbol(string(name.name.name)), T, N)
function gen_fixedsizevector_type(name::Symbol, T::Symbol, N::Int)
    fields = [Expr(:(::), symbol("I_$i"), T) for i = 1:N]
    typename = symbol("FS$name")
    eval(quote
        immutable $(typename){$T} <: $(name){$T}
            $(fields...)
        end
    end)
    typename
end

stagedfunction Base.call{T <: FixedSizeVector, ET}(t::Type{T}, data::ET...)
    Tsuper, Nsuper = super(T).parameters
    N = length(data)
    @assert Nsuper == N "not the right dimension"
    typename = gen_fixedsizevector_type(T, symbol(string(Tsuper.name)), N)
    original_typename = t.name.name
    :($(original_typename)($(typename)(data...)))
end

immutable LOL{T <: FixedSizeVector{Real, 3}} <: FixedSizeVector{Real, 3}
	data::T
end

@show LOL(1,2,3)