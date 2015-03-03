
function generate_arrays(maxSz::Integer, vecname::String, matname::String)
    # expression functions
    vecTyp(n)    = symbol(string(vecname,n))
    vecTypT(n)   = Expr(:curly, vecTyp(n), :T)
    matTyp(r,c)  = symbol(string(matname, r, "x", c))
    matTypT(r,c) = Expr(:curly, matTyp(r,c,), :T)
    elt(i)       = symbol(["x", "y", "z", "w"][i])
    col(i)       = symbol(string("c",i))
    mem(s,e)     = Expr(:.,s,Expr(:quote,e))
    velt(v,i)    = mem(v,elt(i))
    melt(m,i,j)  = mem(mem(m,col(j)),elt(i))

    # vector types
    for sz = 1:maxSz
        local Typ  = vecTyp(sz)
        local TypT = vecTypT(sz)

        # the body of the type definition
        local defn = :(immutable $TypT <: FixedVector{T, $sz} end)
        # the members of the type
        for i = 1:sz
            local e = elt(i)
            push!(defn.args[3].args, :($e::T))
        end
        # instantiate the type definition
        eval(defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::T) = $TypT())
        for i = 1:sz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::T))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        eval(ctorn)
        eval(ctor1)

        # construct or convert from other vector types
        @eval $Typ(a::AbstractVector) = $Typ(ntuple($sz, x-> a[x])...)
        @eval convert{T}(::Type{$TypT}, x::AbstractVector) = $Typ(x)

        # convert to Array
        @eval begin
            function convert{T}(::Type{Vector{T}}, v::$TypT)
                a = Array(T,$sz)
                for i = 1:$sz
                    a[i] = v[i]
                end
                a
            end
        end
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
        eval(defn)

        # unary and n-ary constructors
        ctorn = :($TypT() = $TypT())
        ctor1 = :($TypT(a::$ColTypT) = $TypT())
        for i = 1:cSz
            local arg = symbol(string("a",i))
            push!(ctorn.args[1].args, :($arg::$ColTypT))
            push!(ctorn.args[2].args, arg)
            push!(ctor1.args[2].args, :a)
        end
        eval(ctorn)
        eval(ctor1)

        # construction from a scalar
        @eval $TypT(a::T) = $Typ($ColTyp(a))

        # construct or convert from other matrix types
        @eval $Typ(a::AbstractMatrix) = $Typ(ntuple($cSz, c->
            $ColTyp(ntuple($rSz, r-> a[r,c])...))...)
        @eval convert{T}(::Type{$TypT}, x::AbstractMatrix) = $Typ(x)

        # convert to Array
        @eval begin
            function convert{T}(::Type{Matrix{T}}, m::$TypT)
                a = Array(T,$rSz,$cSz)
                for i = 1:$rSz, j = 1:$cSz
                    a[i,j] = m[i,j]
                end
                a
            end
        end

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
        @eval ctranspose{T}(m::$TypT) = $bdy

        
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
        @eval *{T}(v::$ColTypT,m::$TypT) = $bdy

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
        @eval *{T}(m::$TypT,v::$RowTypT) = $bdy

        # vector-matrix-vector multiplication
        bdy = :(+())
        for i = 1:rSz, j = 1:cSz
            push!(bdy.args, 
                  Expr(:call, :*, 
                       mem(:vl,elt(i)),
                       mem(mem(:m,col(j)),elt(i)),
                       mem(:vr,elt(j))))
        end
        @eval *{T}(vl::$ColTypT,m::$TypT,vr::$RowTypT) = $bdy

        # identity
        bdy = :($TypT())
        for j = 1:cSz
            push!(bdy.args, :(unit($ColTypT,$j)))
        end
        @eval eye{T}(::Type{$TypT}) = $bdy

        # matrix diagonal
        diagSz = min(rSz,cSz)
        diagTypT = vecTypT(diagSz)
        bdy = :($diagTypT())
        for i = 1:diagSz
            push!(bdy.args, mem(mem(:m,col(i)),elt(i)))
        end
        @eval diag{T}(m::$TypT) = $bdy
        
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
        @eval *{T}(m1::$m1T,m2::$m2T) = $bdy
    end

    # matrix determinant and inverse
    for sz in 1:maxSz
        local Typ = matTyp(sz,sz)
        local TypT = matTypT(sz,sz)
        @eval det{T}(a::$TypT) = det(convert(Array{T,2},a))
        @eval inv{T}(a::$TypT) = $Typ(inv(convert(Array{T,2},a)))
    end
end

generate_arrays(4, "Vector", "Matrix")