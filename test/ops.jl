using FixedSizeArrays


square(x) = x*x

const unaryOps = (-, ~, conj, abs, 
                  sin, cos, tan, sinh, cosh, tanh, 
                  asin, acos, atan, asinh, acosh, atanh,
                  sec, csc, cot, asec, acsc, acot,
                  sech, csch, coth, asech, acsch, acoth,
                  sinc, cosc, cosd, cotd, cscd, secd,
                  sind, tand, acosd, acotd, acscd, asecd,
                  asind, atand, radians2degrees, degrees2radians,
                  log, log2, log10, log1p, exponent, exp,
                  exp2, expm1, cbrt, sqrt, square, erf, 
                  erfc, erfcx, erfi, dawson, ceil, floor,
                  trunc, round, significand, lgamma, hypot,
                  gamma, lfact, frexp, modf, airy, airyai,
                  airyprime, airyaiprime, airybi, airybiprime,
                  besselj0, besselj1, bessely0, bessely1,
                  eta, zeta, digamma)

# vec-vec and vec-scalar
const binaryOps = (.+, .-,.*, ./, .\, .^,*,/,
                   .==, .!=, .<, .<=, .>, .>=, +, -,
                   min, max,
                   div, fld, rem, mod, mod1, cmp,
                   atan2, besselj, bessely, hankelh1, hankelh2, 
                   besseli, besselk, beta, lbeta)

const reductions = (sum, prod, minimum,(maximum)
testresult = Dict{Function, Any}()

function Base.filter(fun, x...)
      result = []
      @assert length(unique(map(length, x))) == 1 "all iterables need to have the same length. Lengths given: $(map(length, x))"
      for i=1:length(x[1])
            args = ntuple(j-> x[j][i], length(x))
            fun(args...) && push!(result, args)
      end
      result
end

@show filter(!=, [2,3,4,5], [2,1,4,7])

Base.call{FS <: AbstractFixedSizeArray, T, N}(::Type{FS}, a::Array{T, N}) = AbstractFixedSizeArray{T, N, size(a)}(a...) 
Base.call{FS <: AbstractFixedSizeArray, T, N}(::Type{FS}, a::Array{T, N}) = AbstractFixedSizeArray{T, N, size(a)}(a...) 



nvec{T, N}(x::Array{T,N}) = AbstractFixedSizeArray(x)
     

function testunaray()
      baseline = Any[rand(2), rand(3), rand(3,7), rand(7,2), rand(Float32, 4, 1), rand(Float32, 4), rand(Float64, 4,4)]
      test     = map(nvec, baseline)
      for op in unaryOps
            baeline_result = map(op, baseline)
            test_result    = map(op, test)

            if all(map(==, baeline_result, test_result))
                  testresult[op] = "passed"
            else
                  testresult[op] = ["didn't pass: ", filter(!=, baseline_result, test_result)]
            end
      end
            
end