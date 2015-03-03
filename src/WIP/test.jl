#=
using Romeo, GLAbstraction
obj = visualize(Texture(joinpath(homedir(),"Desktop/test.jpg")))
push!(Romeo.RENDER_LIST, obj)
while Romeo.window.inputs[:open].value
  Romeo.renderloop(Romeo.window)
end
GLFW.Terminate()
=#



function generate_arrays(maxSz::Integer, result)
    
        
    # operations
    const unaryOps = (:-, :~, :conj, :abs, 
                      :sin, :cos, :tan, :sinh, :cosh, :tanh, 
                      :asin, :acos, :atan, :asinh, :acosh, :atanh,
                      :sec, :csc, :cot, :asec, :acsc, :acot,
                      :sech, :csch, :coth, :asech, :acsch, :acoth,
                      :sinc, :cosc, :cosd, :cotd, :cscd, :secd,
                      :sind, :tand, :acosd, :acotd, :acscd, :asecd,
                      :asind, :atand, :radians2degrees, :degrees2radians,
                      :log, :log2, :log10, :log1p, :exponent, :exp,
                      :exp2, :expm1, :cbrt, :sqrt, :square, :erf, 
                      :erfc, :erfcx, :erfi, :dawson, :ceil, :floor,
                      :trunc, :round, :significand, :lgamma, :hypot,
                      :gamma, :lfact, :frexp, :modf, :airy, :airyai,
                      :airyprime, :airyaiprime, :airybi, :airybiprime,
                      :besselj0, :besselj1, :bessely0, :bessely1,
                      :eta, :zeta, :digamma)

    # vec-vec and vec-scalar
    const binaryOps = (:.+, :.-,:.*, :./, :.\, :.^, 
                       :.==, :.!=, :.<, :.<=, :.>, :.>=,
                       :min, :max,
                       :div, :fld, :rem, :mod, :mod1, :cmp,
                       :atan2, :besselj, :bessely, :hankelh1, :hankelh2, 
                       :besseli, :besselk, :beta, :lbeta)
    
    # vec-vec only
    const binaryOps2 = (:+,:-)

    const reductions = ((:sum,:+),(:prod,:*),(:minimum,:min),(:maximum,:max))

    # expression functions
    vecTyp(n)       = symbol(string("Vector",n))
    vecTypT(n)      = Expr(:curly, vecTyp(n), :T)
    matTyp(r,c)     = symbol(string("Matrix",r,"x",c))
    matTypT(r,c)    = Expr(:curly, matTyp(r,c,), :T)
    elt(i)          = symbol(string("e",i))
    col(i)          = symbol(string("c",i))
    mem(s,e)        = Expr(:.,s,Expr(:quote,e))
    velt(v,i)       = mem(v,elt(i))
    melt(m,i,j)     = mem(mem(m,col(j)),elt(i))

    # vector types
    for sz = 1:maxSz
        local Typ = vecTyp(sz)
        local TypT = vecTypT(sz)

        # the body of the type definition
        local defn = :(immutable $TypT <: FixedVector{T, $sz} end)

        # the members of the type
        for i = 1:sz
            local e = elt(i)
            push!(defn.args[3].args, :($e::T))
        end

        # instantiate the type definition
        println(result, defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::T) = $TypT())
        for i = 1:sz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::T))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        println(result, ctorn)
        println(result, ctor1)

        # construct or convert from other vector types
        println(result, :($Typ(a::AbstractVector) = $Typ(ntuple($sz, x-> a[x])...)))
        println(result, :(convert{T}(::Type{$TypT}, x::AbstractVector) = $Typ(x)))

        # convert to Array
        println(result, quote
            function convert{T}(::Type{Vector{T}}, v::$TypT)
                a = Array(T,$sz)
                for i = 1:$sz
                    a[i] = v[i]
                end
                a
            end
        end)
            

        # helper for defining maps
        mapBody(f,j) = begin
            mp = :($Typ())
            for i = 1:sz
                local ff = copy(f)
                ff.args[j] = mem(:v,elt(i))
                push!(mp.args, ff)
            end
            mp
        end

        for op = unaryOps
            local bdy = mapBody(:($op(x)),2)
            println(result, :($op(v::$Typ) = $bdy))
        end

        for op = binaryOps
            local bdy = :($Typ())
            for i = 1:sz
                local mem1 = mem(:v1,elt(i))
                local mem2 = mem(:v2,elt(i))
                push!(bdy.args, Expr(:call,op,mem1,mem2))
            end
            #println(result, :($op(v1::$Typ,v2::$Typ) = $bdy))

            bdy = mapBody(:($op(s,x)),3)
            if op == :.^ # special version for MathConst{:e}
                #println(result, :($op(s::MathConst{:e},m::$Typ) = $bdy))
            end
            if op == :min || op == :max
                #println(result, :($op{T2<:Real}(s::T2,v::$Typ) = $bdy))
            else
                #println(result, :($op(s::Number,v::$Typ) = $bdy))
            end

            bdy = mapBody(:($op(x,s)),2)
            if op == :min || op == :max
                #println(result, :($op{T2<:Real}(v::$Typ,s::T2) = $bdy))
            else
                #println(result, :($op(v::$Typ,s::Number) = $bdy))
            end
        end
        
        for op = binaryOps2
            local bdy = :($Typ())
            for i = 1:sz
                local mem1 = mem(:v1,elt(i))
                local mem2 = mem(:v2,elt(i))
                push!(bdy.args, Expr(:call,op,mem1,mem2))
            end
            #println(result, :($op(v1::$Typ,v2::$Typ) = $bdy))
        end

        for pr = reductions
            local bdy = Expr(:call,pr[2])
            for i = 1:sz
                push!(bdy.args, mem(:v,elt(i)))
            end
            local meth = pr[1]
            println(result, :($meth(v::$Typ) = $bdy))
        end

        # convert to column matrix
        local colMatT = matTypT(sz,1)
        println(result, :(column{T}(v::$TypT) = $colMatT(v)))

        # convert to row matrix
        local rowMat = Expr(:call,matTyp(1,sz))
        for i = 1:sz
            local val = mem(:v,elt(i))
            push!(rowMat.args, :(Vector1{T}($val)))
        end
        println(result, :(row{T}(v::$TypT) = $rowMat))

        # vector norms
        println(result, :(norm{T}(v::$TypT) = sqrt(dot(v,v))))
        println(result, :(norm{T}(v::$TypT,p::Number) = begin
            if p == 1
                sum(abs(v))
            elseif p == 2
                norm(v)
            elseif p == Inf
                max(abs(v))
            else
                norm(copy(v),p)
            end
        end))

        # standard basis vectors
        local bdy = :($TypT())
        for j = 1:sz
            push!(bdy.args, :(i==$j?one(T):zero(T)))
        end
        println(result, :(unit{T}(::Type{$TypT}, i::Integer) = $bdy))

        # diagonal matrix
        sqMatTypT = matTypT(sz,sz)
        bdy = :($sqMatTypT())
        for i = 1:sz
            push!(bdy.args, :(v[$i].*unit($TypT,$i)))
        end
        println(result, :(diagm{T}(v::$TypT) = $bdy))

        # elementwise type conversion
        bdy = mapBody(:(convert(T,x)),3)
        println(result, :(convert{T}(::Type{$TypT}, v::$Typ) = $bdy))

        # some one-liners
        println(result, :(similar{T}(::$TypT, t::DataType, dims::Dims) = Array(t, dims)))
        println(result, :(zero{T}(::Type{$TypT}) = $Typ(zero(T))))
        println(result, :(dot{T}(v1::$TypT,v2::$TypT) = sum(v1.*conj(v2))))
        println(result, :(unit{T}(v::$TypT) = v/norm(v)))
    end

    # matrix types
    for rSz = 1:maxSz, cSz = 1:maxSz
        local Typ = matTyp(rSz,cSz)
        local TypT = matTypT(rSz,cSz)
        local ColTyp = vecTyp(rSz)
        local ColTypT = vecTypT(rSz)
        local RowTyp = vecTyp(cSz)
        local RowTypT = vecTypT(cSz)

        # the body of the type definition
        local defn = :(immutable $TypT <: FixedMatrix{T, $rSz, $cSz} end)

        # the members of the type
        for i = 1:cSz
            local c = col(i)
            push!(defn.args[3].args, :($c::$ColTypT))
        end

        # instantiate the type definition
        println(result, defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::$ColTypT) = $TypT())
        for i = 1:cSz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::$ColTypT))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        println(result, ctorn)
        println(result, ctor1)

        # construction from a scalar
        println(result, :($TypT(a::T) = $Typ($ColTyp(a))))

        # construct or convert from other matrix types
        println(result, :($Typ(a::AbstractMatrix) = $Typ(ntuple($cSz, c->
            $ColTyp(ntuple($rSz, r-> a[r,c])...))...)))
        println(result, :(onvert{T}(::Type{$TypT}, x::AbstractMatrix) = $Typ(x)))

        # convert to Array
        println(result, :( begin
            function convert{T}(::Type{Matrix{T}}, m::$TypT)
                a = Array(T,$rSz,$cSz)
                for i = 1:$rSz, j = 1:$cSz
                    a[i,j] = m[i,j]
                end
                a
            end
        end))


        # ctranspose
        local bdy = Expr(:call, matTypT(cSz,rSz))
        for i = 1:rSz
            local rw = :($RowTypT())
            for j = 1:cSz
                local val = mem(mem(:m,col(j)),elt(i))
                push!(rw.args, :(conj($val)))
            end
            push!(bdy.args, rw)
        end
        println(result, :(ctranspose{T}(m::$TypT) = $bdy))

        # helper for defining maps
        mapBody(f,k) = begin
            local bdy = :($Typ())
            for j = 1:cSz
                local cl = :($ColTyp())
                for i = 1:rSz
                    local ff = copy(f)
                    ff.args[k] = mem(mem(:m,col(j)),elt(i))
                    push!(cl.args, ff)
                end
                push!(bdy.args, cl)
            end
            bdy
        end

        for op = unaryOps
            local bdy = mapBody(:($op(x)),2)
            println(result, :($op(m::$Typ) = $bdy))
        end

        for op = binaryOps
            local bdy = :($Typ())
            for j = 1:cSz
                local cl = :($ColTyp())
                for i = 1:rSz
                    push!(cl.args, 
                          Expr(:call,op,
                               mem(mem(:m1,col(j)),elt(i)),
                               mem(mem(:m2,col(j)),elt(i))))
                end
                push!(bdy.args, cl)
            end
            #println(result, :($op(m1::$Typ,m2::$Typ) = $bdy))

            bdy = mapBody(:($op(s,x)),3)
            if op == :.^ # special version for MathConst{:e}
                #println(result, :($op(s::MathConst{:e},m::$Typ) = $bdy))
            end
            if op == :min || op == :max
                #println(result, :($op{T2<:Real}(s::T2,m::$Typ) = $bdy))
            else
                #println(result, :($op(s::Number,m::$Typ) = $bdy))
            end

            bdy = mapBody(:($op(x,s)),2)
            if op == :min || op == :max
               # println(result, :($op{T2<:Real}(m::$Typ,s::T2) = $bdy))
            else
                #println(result, :($op(m::$Typ,s::Number) = $bdy))
            end
        end
        
        for op = binaryOps2
            local bdy = :($Typ())
            for j = 1:cSz
                local cl = :($ColTyp())
                for i = 1:rSz
                    push!(cl.args, 
                          Expr(:call,op,
                               mem(mem(:m1,col(j)),elt(i)),
                               mem(mem(:m2,col(j)),elt(i))))
                end
                push!(bdy.args, cl)
            end
            println(result, :($op(m1::$Typ,m2::$Typ) = $bdy))
        end

        for pr = reductions
            local bdy = Expr(:call,pr[2])
            for i = 1:rSz, j = 1:cSz
                push!(bdy.args, mem(mem(:m,col(j)),elt(i)))
            end
            local meth = pr[1]
            println(result, :($meth(m::$Typ) = $bdy))
        end

        # vector-matrix multiplication
        bdy = :($RowTypT())
        for j = 1:cSz
            local e = :(+())
            for i = 1:rSz
                push!(e.args, 
                      Expr(:call, :*,
                           mem(:v,elt(i)),
                           mem(mem(:m,col(j)),elt(i))))
            end
            push!(bdy.args, e)
        end
        println(result, :(*{T}(v::$ColTypT,m::$TypT) = $bdy))

        # matrix-vector multiplication
        bdy = :($ColTypT())
        for i = 1:rSz
            local e = :(+())
            for j = 1:cSz
                push!(e.args, 
                      Expr(:call, :*,
                           melt(:m,i,j),
                           velt(:v,j)))
            end
            push!(bdy.args, e)
        end
        println(result, :(*{T}(m::$TypT,v::$RowTypT) = $bdy))

        # vector-matrix-vector multiplication
        bdy = :(+())
        for i = 1:rSz, j = 1:cSz
            push!(bdy.args, 
                  Expr(:call, :*, 
                       mem(:vl,elt(i)),
                       mem(mem(:m,col(j)),elt(i)),
                       mem(:vr,elt(j))))
        end
        println(result, :(*{T}(vl::$ColTypT,m::$TypT,vr::$RowTypT) = $bdy))

        # identity
        bdy = :($TypT())
        for j = 1:cSz
            push!(bdy.args, :(unit($ColTypT,$j)))
        end
        println(result, :(eye{T}(::Type{$TypT}) = $bdy))

        # matrix diagonal
        diagSz = min(rSz,cSz)
        diagTypT = vecTypT(diagSz)
        bdy = :($diagTypT())
        for i = 1:diagSz
            push!(bdy.args, mem(mem(:m,col(i)),elt(i)))
        end
        println(result, :(diag{T}(m::$TypT) = $bdy))

        # elementwise type conversion
        bdy = mapBody(:(convert(T,x)),3)
        println(result, :(convert{T}(::Type{$TypT}, m::$Typ) = $bdy))

        # some one-liners
        println(result, :(similar{T}(::$TypT, t::DataType, dims::Dims) = Array(t, dims)))
        println(result, :(zero{T}(::Type{$TypT}) = $Typ($ColTyp(zero(T)))))
    end

    # matrix-matrix multiplication
    for n = 1:maxSz, p = 1:maxSz, m = 1:maxSz
        local bdy = Expr(:call, matTypT(n,m))
        for j = 1:m
            local c = Expr(:call, vecTypT(n))
            for i = 1:n
                local e = :(+())
                for k = 1:p
                    push!(e.args,
                          Expr(:call, :*,
                               melt(:m1,i,k),
                               melt(:m2,k,j)))
                end
                push!(c.args, e)
            end
            push!(bdy.args, c)
        end
        local m1T = matTypT(n,p)
        local m2T = matTypT(p,m)
        println(result, :(*{T}(m1::$m1T,m2::$m2T) = $bdy))
    end

    # matrix determinant and inverse
    for sz in 1:maxSz
        local Typ = matTyp(sz,sz)
        local TypT = matTypT(sz,sz)

        println(result, :(det{T}(a::$TypT) = det(convert(Array{T,2},a))))
        println(result, :(inv{T}(a::$TypT) = $Typ(inv(convert(Array{T,2},a)))))
    end

    # cross products
    if maxSz >= 2
        println(result, :(cross(a::Vector2,b::Vector2) = a.e1*b.e2-a.e2*b.e1))
    end

    if maxSz >= 3
        println(result, :(cross(a::Vector3,b::Vector3) = Vector3(a.e2*b.e3-a.e3*b.e2, 
                                                     a.e3*b.e1-a.e1*b.e3, 
                                                     a.e1*b.e2-a.e2*b.e1)))
    end
end

result = IOBuffer(Uint8[], true, true)
println(result, readall(open("fixeddefinitions.jl")))
generate_arrays(4, result)
println(result, "end")
fd = open("FixedSizeArrays.jl", "w")
print(fd, takebuf_string(result))
