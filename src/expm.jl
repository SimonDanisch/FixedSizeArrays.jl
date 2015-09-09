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
             $(symbol("a",l,k)) = $(esc(A))[$l,$k]
             $(symbol("a",i,k))  = $(esc(A))[$i,$k] * $(esc(c)) + $(esc(A))[$j,$k] * $(esc(s))
             $(symbol("a",j,k))  = $(esc(A))[$j,$k] * $(esc(c)) - $(esc(A))[$i,$k] * $(esc(s))
        end
        for k in 1:3]...,        quote Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33))) end
    )
end

macro rotate3t(A, c, s, i,j)
  l = setdiff(1:3, [i,j])[1]
  codeblock([
        quote
             $(symbol("a",k,l)) = $(esc(A))[$k,$l]
             $(symbol("a",k,i)) = $(esc(A))[$k,$i] * $(esc(c)) + $(esc(A))[$k,$j] * $(esc(s))
             $(symbol("a",k,j)) = $(esc(A))[$k,$j] * $(esc(c)) - $(esc(A))[$k,$i] * $(esc(s))
        end
        for k in 1:3]...,        quote Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33))) end
    )
end

macro dblrotate3(A, c, s, i,j)
  k = setdiff(1:3, [i,j])[1]
  quote
             $(symbol("a",i,k)) = $(esc(A))[$i,$k] * $(esc(c)) + $(esc(A))[$j,$k] * $(esc(s))
             $(symbol("a",j,k)) = $(esc(A))[$j,$k] * $(esc(c)) - $(esc(A))[$i,$k] * $(esc(s))
             $(symbol("a",k,i)) = $(esc(A))[$k,$i] * $(esc(c)) + $(esc(A))[$k,$j] * $(esc(s))
             $(symbol("a",k,j)) = $(esc(A))[$k,$j] * $(esc(c)) - $(esc(A))[$k,$i] * $(esc(s))
             
             $(symbol("a",j,j)) = $(esc(A))[$j,$j] * abs2($(esc(c))) + $(esc(A))[$i,$i] * abs2($(esc(s))) - ($(esc(A))[$i,$j] + $(esc(A))[$j,$i]) * $(esc(s)) * $(esc(c))
             $(symbol("a",j,i)) = $(esc(A))[$j,$i] * abs2($(esc(c))) - $(esc(A))[$i,$j] * abs2($(esc(s))) + ($(esc(A))[$j,$j] - $(esc(A))[$i,$i]) * $(esc(c)) * $(esc(s))
             $(symbol("a",i,j)) = $(esc(A))[$i,$j] * abs2($(esc(c))) - $(esc(A))[$j,$i] * abs2($(esc(s))) + ($(esc(A))[$j,$j] - $(esc(A))[$i,$i]) * $(esc(c)) * $(esc(s))
             $(symbol("a",i,i)) = $(esc(A))[$i,$i] * abs2($(esc(c))) + $(esc(A))[$j,$j] * abs2($(esc(s))) + ($(esc(A))[$i,$j] + $(esc(A))[$j,$i]) * $(esc(s)) * $(esc(c))
             
             $(symbol("a",k,k)) = $(esc(A))[$k,$k]
             
            Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33))) 
   end
end
function hessenberg3(A)
    c, s, _ = LinAlg.givensAlgorithm(A[2,1],A[3,1])
    @dblrotate3(A, c, s, 2, 3)
end


function reflect3(x::AbstractVector, τ::Number, A) # apply reflector from left

            vA1 = τ'*(A[1, 1] +  x[2]'*A[2, 1] +  x[3]'*A[3, 1])
            a11 = A[1, 1] - vA1
            a21 = A[2, 1] - x[2]*vA1
            a31 = A[3, 1] - x[3]*vA1
            
            vA2 = τ'*(A[1, 2] +  x[2]'*A[2, 2] +  x[3]'*A[3, 2])
            a12 = A[1,2] - vA2
            a22 = A[2,2] - x[2]*vA2
            a32 = A[3,2] - x[3]*vA2
            
            vA3 = τ'*(A[1, 3] +  x[2]'*A[2, 3] +  x[3]'*A[3, 3])
            a13 = A[1,3] - vA3
            a23 = A[2,3] - x[2]*vA3
            a33 = A[3,3] - x[3]*vA3
            Mat(((a11,a21,a31),(a12,a22,a32),(a13,a23,a33)))

end

#francis QR algorithm with single or double shift

function francisdbl{T<:Real}(AA::Mat{3, 3, T})
    A = hessenberg3(AA)
    EPS = eps(T)
    i = 0
    loc = 0 #location of zero
    la1::T = la2::T = la3::T = 0.
    for i in 1:80
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
            v = [x,y,z]
            g = LinAlg.reflector!(v)
           
            A = reflect3(v, g, A)
            A = reflect3(v, g, A')'
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
    dis2 = abs2(A[loc,loc] - A[2,2]) + 4*A[loc,2]*A[2,loc]

    # return eigenvalues (if complex, return real and imaginary part as real numbers)
    if dis2 >= 100EPS # real, real, real
        la2 = 0.5*(spur + copysign(sqrt(dis2), spur))
        if abs(la2) >  eps()
            la3 = det2/la2
        else
            la3 = 0.5*(spur - copysign(sqrt(dis2), spur))
        end
    elseif dis2 >= -100EPS # treat as real
        la3 = la2 = 0.5*spur 
        dis2 = 0.
    else # real +  complex pair
        la2 = 0.5*spur #real part
        la3 = 0.5*sqrt(-dis2) #imaginary part
    end
    
    disc = dis2
    
    # return eigenvalus, discriminant of 2x2 matrix of end result and convergence status
    la1, la2, la3, disc, converged
   
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

# see https://en.wikipedia.org/w/index.php?title=Eigenvalue_algorithm&oldid=671391308#3.C3.973_matrices
# the eigenvalues satisfy eig3 <= eig2 <= eig1

function eigvalssym{T<:Real}(A::Mat{3, 3, T})
    p1 = abs2(A[1,2]) + abs2(A[1,3]) + abs2(A[2,3])
    if p1 == 0
        return A[1,1],  A[2,2], A[3,3]
    else
        q = (A[1,1] + A[2,2] + A[3,3])/3
        p2 = abs2(A[1,1] - q) + abs2(A[2,2] - q) + abs2(A[3,3] - q) + 2 * p1
        p = sqrt(p2 / 6)
        B = (1 / p) * (A - q * eye(Mat{3,3,T}))    
        d = det(B)
        r = d / 2
        
        phi = acos(clamp(r, -one(T), one(T))) / 3 # correcting rounding errors

        eig1 = q + 2 * p * cos(phi)
        eig3 = q + 2 * p * cos(phi + (2*pi/3))
        eig2 = 3 * q - eig1 - eig3   
    end
   return eig1, eig2, eig3
end
    

function expm{T<:Real}(AA::Mat{3, 3, T})
    t = sqrt(sum(AA.^2))
    A = AA/t #normalize

    x1, x2, x3, dis2, conv = francisdbl(A)
    if !conv # if convergence takes to long, call lapack
       return Mat{3,3,T}(expm(convert(Matrix{T}, AA)))
    end 
    dis = dis2 >= 0 ? 1 : -1
    if dis < 0 && abs2(x3) < eps(T)
        dis = 1
        x3 = x2
    end
 
    if dis >= 0
        if abs2(x1 - x2) < eps(T)
            putzer3(t, A, x3, x1)
        elseif abs2(x2 - x3) < eps(T)
            putzer3(t, A, x1, x2)
        elseif abs2(x1 - x3) < eps(T)
            putzer3(t, A, x2, x1)            
        else
            putzer3(t, A, x1, x2, x3)
        end
    else 
        putzer3(t, A, x1, x2 + im*x3)
    end
end    

function putzer3{T}(t, A::Mat{3,3,T}, la1, la2) #la2 == la3, maybe la1 == la2 == la3
        if abs(la1 - la2) < 10eps(T)  # la1 == la2 == la3
            r1 = exp(t* la1)
            r2 = t * exp(t* la1)
            r3 = 0.5 * t * t * exp(t * la1)
            R = r2 - la2 * r3

        else
            r1 = exp(t * la1)
            r2 = (expm1(t * la1) - expm1(t * la2)) / (la1 - la2) # la1 != la2
            r3 = (exp(t * la1) - exp(t * la2)*(1. + t*(la1-la2))) / abs2(la1-la2)
            R = r2 - la2*r3

       	end
        r3 * A*A +  (R - la1*r3)*A +  (r1 - la1 * R)*eye(Mat{3, 3, Float64})
end
function putzer3{T}(t, A::Mat{3,3,T}, la1, la2::Real, la3::Real) # three distinct real eigenvalues
        r1 = exp(t * la1)
        r2 = (expm1(t * la1) - expm1(t * la2))/ (la1 - la2) 
        r3 = -((la2 - la3)*exp(t*la1) + (la3 - la1)*exp(t*la2) + (la1 - la2)*exp(t*la3)) / ((la1 - la2)*(la2 - la3)*(la3 - la1))
        R = r2 - la2*r3

        r3 * A*A +  (R - la1*r3)*A +  (r1 - la1 * R)*eye(Mat{3,3,T})
        
end
function putzer3{T}(t, A::Mat{3,3,T}, la1, la2::Complex) # one real and a complex conjugate pair of eigenvalues, all different
        la3 = la2'
        r1 = exp(t * la1)
        r2 = (expm1(t * la1) - expm1(t * la2)) / (la1 - la2) # possibly complex
        r3 = real(-((la2 - la3)*exp(la1*t) + (la3 - la1)*exp(t*la2) + (la1 - la2)*exp(t*la3)) / ((la1 - la2)*(la2 - la3)*(la3 - la1))) # r3 is real, see short comm. by r.r. huilgol
        R = real(r2 - la2*r3) # also real
     
        r3 * A*A +  (R - la1*r3)*A +  (r1 - la1 * R)*eye(Mat{3,3,T})
   
end




