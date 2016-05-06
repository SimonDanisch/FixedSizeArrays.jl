### expm

function codeblock(stmts...)
    stmts_ = Any[]
    for s in stmts
        if !(s == nothing)
            push!(stmts_, s)
        end
    end
    isempty(stmts_) ? nothing : Expr(:block, stmts_...)
end


macro rotate3(A, c, s, i,j)
  l = setdiff(1:3, [i,j])[1]
  codeblock([
        quote
             $(Symbol("a$l$k")) = $(esc(A))[$l,$k]
             $(Symbol("a$i$k"))  = $(esc(A))[$i,$k] * $(esc(c)) + $(esc(A))[$j,$k] * $(esc(s))
             $(Symbol("a$j$k"))  = $(esc(A))[$j,$k] * $(esc(c)) - $(esc(A))[$i,$k] * $(esc(s))
        end
        for k in 1:3]...,        quote Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33))) end
    )
end

macro rotate3t(A, c, s, i,j)
  l = setdiff(1:3, [i,j])[1]
  codeblock([
        quote
             $(Symbol("a$k$l")) = $(esc(A))[$k,$l]
             $(Symbol("a$k$i")) = $(esc(A))[$k,$i] * $(esc(c)) + $(esc(A))[$k,$j] * $(esc(s))
             $(Symbol("a$k$j")) = $(esc(A))[$k,$j] * $(esc(c)) - $(esc(A))[$k,$i] * $(esc(s))
        end
        for k in 1:3]...,        quote Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33))) end
    )
end

macro dblrotate3(A, c, s, i,j)
    k = setdiff(1:3, [i,j])[1]
    quote
        $(Symbol("a$i$k")) = $(esc(A))[$i,$k] * $(esc(c)) + $(esc(A))[$j,$k] * $(esc(s))
        $(Symbol("a$j$k")) = $(esc(A))[$j,$k] * $(esc(c)) - $(esc(A))[$i,$k] * $(esc(s))
        $(Symbol("a$k$i")) = $(esc(A))[$k,$i] * $(esc(c)) + $(esc(A))[$k,$j] * $(esc(s))
        $(Symbol("a$k$j")) = $(esc(A))[$k,$j] * $(esc(c)) - $(esc(A))[$k,$i] * $(esc(s))

        $(Symbol("a$j$j")) = $(esc(A))[$j,$j] * abs2($(esc(c))) + $(esc(A))[$i,$i] * abs2($(esc(s))) - ($(esc(A))[$i,$j] + $(esc(A))[$j,$i]) * $(esc(s)) * $(esc(c))
        $(Symbol("a$j$i")) = $(esc(A))[$j,$i] * abs2($(esc(c))) - $(esc(A))[$i,$j] * abs2($(esc(s))) + ($(esc(A))[$j,$j] - $(esc(A))[$i,$i]) * $(esc(c)) * $(esc(s))
        $(Symbol("a$i$j")) = $(esc(A))[$i,$j] * abs2($(esc(c))) - $(esc(A))[$j,$i] * abs2($(esc(s))) + ($(esc(A))[$j,$j] - $(esc(A))[$i,$i]) * $(esc(c)) * $(esc(s))
        $(Symbol("a$i$i")) = $(esc(A))[$i,$i] * abs2($(esc(c))) + $(esc(A))[$j,$j] * abs2($(esc(s))) + ($(esc(A))[$i,$j] + $(esc(A))[$j,$i]) * $(esc(s)) * $(esc(c))

        $(Symbol("a$k$k")) = $(esc(A))[$k,$k]

        Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33)))
    end
end
function hessenberg3(A)
    c, s, _ = LinAlg.givensAlgorithm(A[2,1],A[3,1])
    @dblrotate3(A, c, s, 2, 3)
end


function reflect3{T}(x1,x2,x3, τ::Number, A::Mat{3, 3, T}) # apply reflector from left

            vA1 = τ'*(A[1, 1] +  x2'*A[2, 1] +  x3'*A[3, 1])
            a11 = A[1,1] - vA1
            a21 = A[2,1] - x2*vA1
            a31 = A[3,1] - x3*vA1

            vA2 = τ'*(A[1, 2] +  x2'*A[2, 2] +  x3'*A[3, 2])
            a12 = A[1,2] - vA2
            a22 = A[2,2] - x2*vA2
            a32 = A[3,2] - x3*vA2

            vA3 = τ'*(A[1, 3] +  x2'*A[2, 3] +  x3'*A[3, 3])
            a13 = A[1,3] - vA3
            a23 = A[2,3] - x2*vA3
            a33 = A[3,3] - x3*vA3
            Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33)))::Mat{3, 3, T}

end

function reflect3t{T}(x1,x2,x3, τ::Number, A::Mat{3, 3, T}) # apply reflector from left

            vA1 = τ'*(A[1, 1] +  x2'*A[1, 2] +  x3'*A[1, 3])
            a11 = A[1,1] - vA1
            a21 = A[1,2] - x2*vA1
            a31 = A[1,3] - x3*vA1

            vA2 = τ'*(A[2, 1] +  x2'*A[2, 2] +  x3'*A[2, 3])
            a12 = A[2,1] - vA2
            a22 = A[2,2] - x2*vA2
            a32 = A[2,3] - x3*vA2

            vA3 = τ'*(A[3, 1] +  x2'*A[3, 2] +  x3'*A[3, 3])
            a13 = A[3,1] - vA3
            a23 = A[3,2] - x2*vA3
            a33 = A[3,3] - x3*vA3
            Mat(((a11,a12,a13),(a21,a22,a23),(a31,a32,a33)))::Mat{3, 3, T}

end


#francis QR algorithm with single or double shift

function francisdbl{T<:Real}(AA::Mat{3, 3, T})
    A = hessenberg3(AA)
    EPS = eps(T)
    i = 0
    loc = 0 #location of zero
    la1::T = la2::T = la3::T = 0.
    for i in 1:200
        if abs(A[3,2]) <= EPS * (abs(A[2,2])+abs(A[3,3]))
            la1 = A[3,3]
            loc = 1
            break
        elseif abs(A[2,1]) <= EPS * (abs(A[1,1])+abs(A[2,2]))
            la1 = A[1,1]
            loc = 3
            break
        end
        l = 3
        # if real ev in sub 2x2 matrix do twice a raleigh shift
        if 0 < abs2(A[l,l] - A[2,2]) + 4*A[l,2]*A[2,l]
            mu = A[l,l]
            A = A - mu*(eye(typeof(A))::Mat{3, 3, T})
            c, s, _ = LinAlg.givensAlgorithm(A[1,1],A[2,1])::Tuple{T,T,T}
            A2 = @rotate3(A, c, s, 1, 2)
            c2, s2, _ = LinAlg.givensAlgorithm(A2[2,2],A2[3,2])::Tuple{T,T,T}
            A = @rotate3(A2, c2, s2, 2, 3)

            A2 = @rotate3t(A, c, s, 1, 2)
            A = @rotate3t(A2, c2, s2, 2, 3)

            c, s, _ = LinAlg.givensAlgorithm(A[1,1],A[2,1])::Tuple{T,T,T}
            A2 = @rotate3(A, c, s, 1, 2)
            c2, s2, _ = LinAlg.givensAlgorithm(A2[2,2],A2[3,2])::Tuple{T,T,T}
            A = @rotate3(A2, c2, s2, 2, 3)

            A2 = @rotate3t(A, c, s, 1, 2)
            A = @rotate3t(A2, c2, s2, 2, 3)
            A = A + mu*(eye(typeof(A))::Mat{3, 3, T})
        else #do an implicit francis double shift

            x = (A[1,1]*A[1,1] + A[1,2]*A[2,1] + A[3,1]*A[1,3]) -  A[1,2]*A[2,1] -
                (A[3,3] + A[2,2])*A[1,1]  + A[3,3] * A[2,2] - A[2,3] * A[3,2]
            y = A[2,1]*(A[1,1] - A[3,3])
            z = A[2,1]*A[3,2]

            ξ1 = x
            normu = abs2(x) + abs2(y) + abs2(z)
            if normu == zero(normu)
               g = zero(ξ1/normu)
            else
                normu = sqrt(normu)
                ν = copysign(normu, real(ξ1))
                ξ1 += ν
                x = -ν
                y /= ξ1
                z /= ξ1
                g = ξ1/ν
            end
            A = reflect3(x,y,z, g, A)
            A = reflect3t(x,y,z, g, A)
            A = hessenberg3(A)
        end
    end

    # check convergence
    if loc == 0
        converged = false
        loc = 1
    else
        converged = true
    end

    # compute eigenvalues of 2x2 submatrix
    spur = A[loc,loc] + A[2,2]
    det2 = A[loc,loc] * A[2,2] - A[2,loc] * A[loc,2]
    dis = abs2(A[loc,loc] - A[2,2]) + 4*A[loc,2]*A[2,loc]


    # return eigenvalue, discriminant of 2x2 matrix of end result and convergence status
    la1, 0.5*spur, dis/4, converged

end


# default version using conversion to and from Matrix type
expm{N,T}(m::Mat{N,N,T}) = Mat{N,N,T}(expm(convert(Matrix{T}, m)))

expm{T}(A::Mat{1, 1, T}) = Mat{1, 1, T}(((expm(A[1,1]),),))
function expm{T<:Complex}(A::Mat{2, 2, T})
 	a = A[1,1]
	b = A[1,2]
	c = A[2,1]
	d = A[2,2]

	z = sqrt((a-d)*(a-d) + 4.0*b*c )
	e = exp(a/2.0 + d/2.0 - z/2.0)
	f = exp(a/2.0 + d/2.0 + z/2.0)

    Mat{2, 2, T}(
        ( -(e*(a - d - z))/(2.0* z) + (f*(a - d + z))/(2.0* z), -((e * c)/z) + (f * c)/z),
        ( -((e * b)/z) + (f * b)/z,	-(e*(-a + d - z))/(2.0* z) + (f*(-a + d + z))/(2.0* z))
    )
end

# presumably better without resorting to complex numbers, but for now...
function expm{T<:Real}(A::Mat{2, 2, T})
 	a = A[1,1]
	b = A[1,2]
	c = A[2,1]
	d = A[2,2]

	z = sqrt(Complex((a-d)*(a-d) + 4.0*b*c))
	e = exp(a/2.0 + d/2.0 - z/2.0)
	f = exp(a/2.0 + d/2.0 + z/2.0)

    Mat{2, 2, T}(
        ( real(-(e*(a - d - z))/(2.0* z) + (f*(a - d + z))/(2.0* z)), real(-((e * c)/z) + (f * c)/z)),
        ( real(-((e * b)/z) + (f * b)/z), real(-(e*(-a + d - z))/(2.0* z) + (f*(-a + d + z))/(2.0* z)))
    )
end


function expm{T<:Real}(AA::Mat{3, 3, T})
    t = sqrt(sum(AA.^2))

    A = AA/t #normalize

    la1, x, dis, conv = francisdbl(A)
    if !conv # if convergence takes to long, call lapack
       return Mat{3,3,T}(expm(convert(Matrix{T}, AA)))
    end
    putzer(t, A, la1, x, dis)
end

function dexp(t, x, y)
    if abs2(x-y) > eps()
        2exp(t*(x+y)/2)*sinh(t*(x-y)/2)/((x-y))
    else
        exp(t*(x+y)/2)*(t+(t^3*(x-y)^2)/24)
    end
end

function ddexp(t, x, y, z)
    d, i = findmax([abs2(y-x), abs2(z-y),abs2(z-x)])
    if d < eps()
        exp(t*(x+y+z)/3)*t^2/2
    elseif i == 1
        (dexp(t, y, z) - dexp(t, x,z))/(y-x)
    elseif i==2
        (dexp(t, z, x) - dexp(t, y,x))/(z-y)
    elseif i==3
        (dexp(t, z, y) - dexp(t, x,y))/(z-x)
    end
end

function putzer{T}(t, A::Mat{3,3,T}, x, y, d)
    deltaabs = sqrt(abs(d))
    delta = sqrt(complex(d))
    r1 = real(exp(t * (y+delta)))
    la1 = (y+delta)
    la2 = (y-delta)
    la3 = x
    r2 = real(dexp(t, la1, la2))
    r3 = real(ddexp(t, la1, la2, la3))
    R = r2 - la2*r3
    r3 * A*A +  real(R - la1*r3)*A +  real(r1 - la1 * R)*eye(Mat{3,3,T})
end
