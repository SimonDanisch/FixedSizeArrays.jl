# FixedSizeArrays

[![Join the chat at https://gitter.im/SimonDanisch/FixedSizeArrays.jl](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/SimonDanisch/FixedSizeArrays.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl.svg?branch=master)](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl)

#### This package is 0.4 only


#### Usage and advantages:
FixedSizeArrays is giving any composite type array like behavior by inheriting from FixedSizeArrays.
So you can do something like this:
```Julia
immutable Vec3{T} <: FixedSizeVector{T, 3}
r::T
g::T
b::T
end
Vec3(0) # constructor with 0 argument already defined. Actually, this only works if you created them
#  with @gen_fixed_size_vector, due to ambiguity issues :(
rand(Vec3{Int})+sin(Vec3(0,2,2)) # a lot of array functions are already defined
immutable Mat3x3{T} FixedMatrix{T, 3,3}
v_1::T
...
v_9::T
end
eye(Mat3x3{Float32}) * rand(Vec3{Float32}) # will also "just work"
```
This is expendable to a lot of other areas.
You can define color types (see ColorTypes.jl) the same way, and arbitrary other point types like normals, vertices etc.
As they all inherit from FixedSizeArray, it's very easy to handle them in the same way. So I'm using this for my GPU array types, which can take any fixedsizearray, if its a color, a point or what not, because I can be sure that all the important functions are defined.
For some more advantages, you can take a look at [MeshIO](https://github.com/JuliaIO/MeshIO.jl)
Also, because it's so easy to define different types like Point3, RGB, HSV or Normal3, one can create customized code for these types via multiple dispatch. This is great for visualizing data, as you can offer default visualizations based on the type.
Without FixedSizeArrays, this would end up in a lot of types which would all need to define the same functions over and over again.



#### TODO's

- [ ] Core Array
	- [x] basic array interface
	- [ ] Inherit from DenseArray (a lot of warnings is caused by this)
	- [ ] use tuples as a basis
- [ ] Indexing:
	- [x] multidimensional access
	- [x] colon access for matrices
	- [ ] multidimensional colon access
	- [ ] setindex!
	- [ ] setindex!/getindex for arrays of FSA (e.g. easy acces to single fields) 
	- [ ] access via dimension type (Red -> redchannel)
- [ ] Constructor
	- [x] generic constructor for arbitrary Nvectors
	- [x] fast constructor for arbitrary types
	- [x] different constructors for ease of use (zero, eye, from other FSAs, etc...) (could be more)
- [ ] Functions
	- [x] all kinds of unary/binary operators
	- [x] matrix multiplication 
	- [x] matrix functions (inv, transpose, etc...) (could be more)




#### Acknowledgements
[ImmutableArrays](https://github.com/twadleigh/ImmutableArrays.jl) by [twadleigh](https://github.com/twadleigh) was the package that got me going and gave the initial inspirations.
There has been quite a few discussions on [JuliaLang/julia#7568](https://github.com/JuliaLang/julia/pull/7568) shaping the implementation.
Also, [aaalexandrov](https://github.com/aaalexandrov) supplied some code and inspirations.
