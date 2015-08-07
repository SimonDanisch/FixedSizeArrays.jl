using FixedSizeArrays


square(x) = x*x

const unaryOps = (-, ~, conj, abs, 
                  sin, cos, tan, sinh, cosh, tanh, 
                  asin, acos, atan, asinh, acosh, atanh,
                  sec, csc, cot, asec, acsc, acot,
                  sech, csch, coth, asech, acsch, acoth,
                  sinc, cosc, cosd, cotd, cscd, secd,
                  sind, tand, acosd, acotd, acscd, asecd,
                  asind, atand, rad2deg, deg2rad,
                  log, log2, log10, log1p, exponent, exp,
                  exp2, expm1, cbrt, sqrt, erf, 
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




facts("mapping operators") do 
    context("unary: ") do 
        baseline = Any[rand(2), rand(3), rand(3,7), rand(7,2), rand(Float32, 4, 1), rand(Float32, 4), rand(Float64, 4,4)]
        test     = map(Vec, baseline)
        for op in unaryOps
            for t in test
                if typeof(op(t[1])) == eltype(t) # map does not handle operators that change the type yet.
                    v = op(t)

                    for i=1:length(v)
                        @fact v[i] => op(t[i])
                    end 
                end
            end
        end
    end

    context("binary: ") do
        baseline = Any[rand(2), rand(3), rand(3,7), rand(7,2), rand(Float32, 4, 1), rand(Float32, 4), rand(Float64, 4,4)]
        baseline2 = Any[rand(2), rand(3), rand(3,7), rand(7,2), rand(Float32, 4, 1), rand(Float32, 4), rand(Float64, 4,4)]
        test1     = map(Vec, baseline)
        test2     = map(Vec, baseline2)
        for op in binaryOps
            for i=1:length(test)
                v1 = test1[i]
                v2 = test2[i]
                r = op(v1, v2)
                for j=1:length(v)
                    @fact r[j] => op(v1[j], v2[j])
                end 
            end
        end
    end
end