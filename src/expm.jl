### expm

#todo replace I3 by eye(Mat{3,3,T}). Only that is slooow
const I3 = Mat([1.0 0.0 0.0
                                        0.0 1.0 0.0
                                        0.0 0.0 1.0])

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
        B = (1 / p) * (A - q * I3)    
        d = det(B)
        r = d / 2
        
        phi = acos(clamp(r, -one(T), one(T))) / 3 # correcting rounding errors

        eig1 = q + 2 * p * cos(phi)
        eig3 = q + 2 * p * cos(phi + (2*pi/3))
        eig2 = 3 * q - eig1 - eig3   
    end
   return eig1, eig2, eig3
end
    

function expm{T<:Real}(A::Mat{3, 3, T})
    if A[1,2] != A[2,1] || A[1,3] != A[3,1] || A[2,3] != A[3,2]
        Mat{3,3,T}(expm(convert(Matrix{T}, A)))
    else
        x1, x2, x3 = eigvalssym(A)
        #x1, x2, x3 = eigvals(Matrix(A))
        if abs2(x1 - x2) < eps(T)
            putzer3(one(T), A, x3, x1)
        elseif abs2(x2 - x3) < eps(T)
            putzer3(one(T), A, x1, x2)
        elseif abs2(x1 - x3) < eps(T)
            putzer3(one(T), A, x2, x1)            
        else
            putzer3(one(T), A, x1, x2, x3)
        end
    end
end    

function putzer3(t, A, la1, la2) #la2 == la3, maybe la1 == la2 == la3
        if abs(la1 - la2) < eps()  # la1 == la2 == la3
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
function putzer3(t, A, la1, la2::Real, la3::Real) # three distinct real eigenvalues
        r1 = exp(t * la1)
        r2 = (expm1(t * la1) - expm1(t * la2))/ (la1 - la2) 
        r3 = -((la2 - la3)*exp(t*la1) + (la3 - la1)*exp(t*la2) + (la1 - la2)*exp(t*la3)) / ((la1 - la2)*(la2 - la3)*(la3 - la1))
        R = r2 - la2*r3

        r3 * A*A +  (R - la1*r3)*A +  (r1 - la1 * R)*I3
        
end
function putzer3(t, A, la1, la2::Complex) # one real and a complex conjugate pair of eigenvalues, all different
        la3 = la2'
        r1 = exp(t * la1)
        r2 = (expm1(t * la1) - expm1(t * la2)) / (la1 - la2) # possibly complex
        r3 = real(-((la2 - la3)*exp(la1*t) + (la3 - la1)*exp(t*la2) + (la1 - la2)*exp(t*la3)) / ((la1 - la2)*(la2 - la3)*(la3 - la1))) # r3 is real, see short comm. by r.r. huilgol
        R = real(r2 - la2*r3) # also real
     
        r3 * A*A +  (R - la1*r3)*A +  (r1 - la1 * R)*I3
   
end




