immutable FixedUnitRange{From, To} end
immutable FixedStepRange{From, Step, To} end
immutable LOL
    x::Float32
end
Base.endof(a::LOL) = 1
Base.getindex{From, To}(a::LOL, b::FixedUnitRange{From, To}) = println(From, " ", To)
Base.getindex(a::LOL, b::UnitRange)      = a[FixedUnitRange{first(b), last(b)}()]
a(x::Colon) = println(x)
a(:)
LOL(1f0)[:]