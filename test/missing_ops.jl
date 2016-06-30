#=
	This file implements FSA operations on other datatypes
	which are not defined in the respective package.
	This way, the benchmark functions stay simple
=#
import Base: +, -, .*, ./, *, /
ops = (:(+), :(-), :(.*), :(./), :(*), :(/))
for N=1:10, op in ops
	@eval begin
		@inline function $(op){T}(a::NTuple{$N, T}, b::NTuple{$N, T})
			@inbounds return tuple($([:($(op)(a[$i], b[$i])) for i=1:N]...))
		end
		@inline function $(op){T}(a::NTuple{$N, T}, b)
			@inbounds return tuple($([:($(op)(a[$i], b)) for i=1:N]...))
		end
		@inline function $(op){T}(a, b::NTuple{$N, T})
			@inbounds return tuple($([:($(op)(a, b[$i])) for i=1:N]...))
		end

		@inline function $(op){T<:VecElement}(a::NTuple{$N, T}, b::NTuple{$N, T})
			@inbounds return tuple($([:(VecElement($(op)(a[$i].value, b[$i].value))) for i=1:N]...))
		end
		@inline function $(op){T<:VecElement}(a::NTuple{$N, T}, b)
			@inbounds return tuple($([:(VecElement($(op)(a[$i].value, b))) for i=1:N]...))
		end
		@inline function $(op){T<:VecElement}(a, b::NTuple{$N, T})
			@inbounds return tuple($([:(VecElement($(op)(a, b[$i].value))) for i=1:N]...))
		end
	end
end
unary_ops = (:(+), :(-))
for N=1:10, op in ops
	@eval begin
		@inline function $(op){T}(a::NTuple{$N, T})
			@inbounds return tuple($([:($(op)(a[$i])) for i=1:N]...))::NTuple{$N, T}
		end
		@inline function $(op){T<:VecElement}(a::NTuple{$N, T})
			@inbounds return tuple($([:(VecElement($(op)(a[$i].value))) for i=1:N]...))::NTuple{$N, T}
		end
	end
end
Base.zero{T}(::VecElement{T}) = VecElement(zero(T))
function +(a::VecElement, b::VecElement)
	VecElement(a.value + b.value)
end
function *(a::VecElement, b::VecElement)
	VecElement(a.value * b.value)
end
function Base.sum(a::Tuple)
	reduce(+, a)
end

@inline function (*)(a::Vector, b::Vector)
	a .* b
end

@inline function (*)(a::FixedSizeArrays.FixedVector, b::FixedSizeArrays.FixedVector)
	a .* b
end
@inline function (./)(a::SIMD.Vec, b::SIMD.Vec)
	a / b
end
@inline function (./){T}(a::SIMD.Vec, b::T)
	a / b
end


@inline function (/){N, T<:Integer, T2<:Integer}(a::SIMD.Vec{N, T}, b::SIMD.Vec{N, T2})
	TT = Base.promote_op(/, T, T2)
	SIMD.Vec{N, TT}(a) / SIMD.Vec{N, TT}(b)
end

@inline function (/){N,T1<:Integer,T2<:Integer}(a::SIMD.Vec{N, T1}, b::T2)
	TT =  Base.promote_op(/, T1, T2)
	SIMD.Vec{N, TT}(a) / TT(b)
end





