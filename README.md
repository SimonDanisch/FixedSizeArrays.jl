# FixedSizeArrays

[![Join the chat at https://gitter.im/SimonDanisch/FixedSizeArrays.jl](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/SimonDanisch/FixedSizeArrays.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl.svg?branch=master)](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl)
[![Coverage Status](https://coveralls.io/repos/SimonDanisch/FixedSizeArrays.jl/badge.svg?branch=master)](https://coveralls.io/r/SimonDanisch/FixedSizeArrays.jl?branch=master)
[![codecov.io](http://codecov.io/github/SimonDanisch/FixedSizeArrays.jl/coverage.svg?branch=master)](http://codecov.io/github/SimonDanisch/FixedSizeArrays.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/k6bqy1h4jk322cg6/branch/master?svg=true)](https://ci.appveyor.com/project/SimonDanisch/fixedsizearrays-jl/branch/master)

#### This package doesn't support 0.3 and it's not planned to update it to 0.6. Use [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl/pull/138) instead.

#### Packages that use FixedSizeArrays:
[GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl)

#### Usage and advantages:
FixedSizeArrays is giving any composite type array like behavior by inheriting from FixedSizeArrays.
So you can do something like this:
```Julia
immutable RGB{T} <: FixedVectorNoTuple{3, T}
    r::T
    g::T
    b::T
end
immutable Vec{N, T} <: FixedVector{N, T} # defined in FixedSizeArrays already
    _::NTuple{N, T}
end
Vec{3, Float32}(77) # constructor with 1 argument already defined
rand(Vec{3, Float64})+sin(Vec(0.,2.,2.)) # a lot of array functions are already defined
#There is also a matrix type
eye(Mat{3,3,Float32}) * rand(Vec{3, Float32}) # will also "just work"
a = Vec(1,2,3)[1:2] # returns (1,2)
```
Note that all of the above types are stack allocated and the speed of operations should be very fast!
If you find operations to be slow, please file a bug report!

FixedSizeArrays can be used in a lot of different ways.
You can define color types the same way, and arbitrary other point types like normals, vertices etc.
As they all inherit from FixedSizeArray, it's very easy to handle them in the same way.

For some more advantages, you can take a look at [MeshIO](https://github.com/JuliaIO/MeshIO.jl).

Because it's so easy to define different types like Point3, RGB, HSV or Normal3, one can create customized code for these types via multiple dispatch. This is great for visualizing data, as you can offer default visualizations based on the type.
Without FixedSizeArrays, this would end up in a lot of types which would all need to define the same functions over and over again.

#### FixedArray abstract types

The package provides several abstract types:

  * `FixedArray{T,NDim,SIZE}` is the abstract base type for all fixed
    arrays.  `T` and `NDim` mirror the eltype and number of dimension type
    parameters in `AbstractArray`.  In addition there's a `SIZE` Tuple which
    defines the extent of each fixed dimension as an integer.

There's some convenient type aliases:

  * `FixedVector{N,T}` is a convenient type alias for a one dimensional fixed
    vector of length `N` and eltype `T`.
  * `FixedMatrix{N,M,T}` is a convenient type alias for a two dimensional fixed
    matrix of size `(N,M)` and eltype `T`.

Finally there's an abstract type `FixedVectorNoTuple{N, T}` for use when you'd
like to name the fields of a `FixedVector` explicitly rather than accessing them
via an index.


#### FixedArray concrete types

The package currently provides three concrete FixedArray types

  * `Vec{N,T}` is a length `N` vector of eltype `T`.
  * `Mat{N,M,T}` is an `N×M` matrix of eltype `T`

These two types are intended to behave the same as `Base.Vector` and
`Base.Matrix`, but with fixed size.  That is, the interface is a convenient
union of elementwise array-like functionality and vector space / linear algebra
operations.  Hopefully we'll have more general higher dimensional fixed size
containers in the future (note that the total number of elements of a higher
dimensional container quickly grows beyond the size where having a fixed stack
allocated container really makes sense).

  * `Point{N,T}` is a position type which is structurally identical to `Vec{N,T}`.

Semantically `Point{N,T}` should be used to represent position in an
`N`-dimensional Cartesian space.  The distinction between this and `Vec` is
particularly relevant when overloading functions which deal with geometric data.
For instance, a geometric transformation applies differently depending on
whether you're transforming a *position* (`Point`) versus a *direction* (`Vec`).


#### User-supplied functions for FixedArray subtypes

Most array functionality comes for free when inheriting from one of the abstract
types `FixedArray`, `FixedVector`, `FixedMatrix`, or `FixedVectorNoTuple`.
However, the user may want to overload a few things.  At the moment,
`similar_type` is the main function you may want to customize.  The signature is

```julia
similar_type{FSA<:FixedArray, T, NDim}(::Type{FSA}, ::Type{T}, sz::NTuple{NDim,Int})
```

This is quite similar to `Base.similar` but the first argument is a type rather
than a value.  Given a custom FixedArray type, eltype and size, this function
should return a similar output type which will be used to store the results of
elementwise operations, general `map()` invocations, etc.

By default, `similar_type` introspects `FSA` to determine whether it can be
reparameterized by both `eltype(FSA) == T` and `size(FSA) == sz`.  If not, the
canonical concrete FixedArray type (a `Vec` or `Mat`) are returned by calling
the fallback `similar_type(FixedArray, T, sz)`.  Sometimes this may not make
sense for your custom FixedArray subtype.

For example, suppose you define the type `RGB{T}` as above, and you'd prefer
relational operators to return a `Vec{3,Bool}` as a mask rather than an
`RGB{Bool}`.  In this case you could write something like:

```julia
function FixedSizeArrays.similar_type{FSA<:RGB,T}(::Type{FSA}, ::Type{T}, n::Tuple{Int})
    n == (3,) && T != Bool ? RGB{T} : similar_type(FixedArray, T, n)
end
```

We then have `RGB(1,2,3) .< RGB(2,2,2) === Vec{3,Bool}(true,false,false)`.

Note that `similar_type` isn't type stable in julia-0.4.  For the internal
use in `FixedSizeArrays` (type deduction inside `@generated` functions) this
isn't a problem, but you may want to annotate your custom overlads with
`Base.@pure` if you're using julia-0.5 and you want to use `similar_type` in a
normal function.


#### Roadmap
* improve coverage
* incorperate https://github.com/StephenVavasis/Modifyfield.jl
* improve API and consistency

#### TODO's

- [ ] Core Array
	- [x] basic array interface
	- [ ] Inherit from DenseArray (a lot of warnings is caused by this)
	- [x] use tuples as a basis
- [ ] Indexing:
	- [x] multidimensional access
	- [x] colon access for matrices
	- [x] multidimensional colon access
	- [ ] setindex!
	- [ ] setindex!/getindex for arrays of FSA (e.g. easy acces to single fields) 
	- [x] access slices e.g. Matrix{RGBA} -> Matrix{Red} (sort of)
- [ ] Constructor
	- [x] generic constructor for arbitrary Nvectors
	- [x] fast constructor for arbitrary types
	- [x] parsing constructor e.g Vec{3, Float32}(["23.", "23.", "0.23"])
	- [x] different constructors for ease of use (zero, eye, from other FSAs, etc...) (could be more)
	- [ ] clean up constructor code (very messy since its hard to write constructors for abstract types)
- [x] Functions
	- [x] all kinds of unary/binary operators
	- [x] matrix multiplication 
	- [x] matrix functions (inv, transpose, etc...) (could be more)




#### Acknowledgements
[ImmutableArrays](https://github.com/twadleigh/ImmutableArrays.jl) by [twadleigh](https://github.com/twadleigh) was the package that got me going and gave the initial inspirations.
There has been quite a few discussions on [JuliaLang/julia#7568](https://github.com/JuliaLang/julia/pull/7568) shaping the implementation.
Also, [aaalexandrov](https://github.com/aaalexandrov) supplied some code and inspirations.
Big thanks to all the other [contributors](https://github.com/SimonDanisch/FixedSizeArrays.jl/graphs/contributors) !

