# FixedSizeArrays

[![Join the chat at https://gitter.im/SimonDanisch/FixedSizeArrays.jl](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/SimonDanisch/FixedSizeArrays.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl.svg?branch=master)](https://travis-ci.org/SimonDanisch/FixedSizeArrays.jl)

- [ ] Core Array
	- [x] basic array interface
	- [ ] Inherit from DenseArray (a lot of warnings is caused by this)
- [ ] Indexing:
	- [x] multidimensional access
	- [x] colon access for matrices
	- [ ] multidimensional colon access
	- [ ] setindex!
	- [ ] setindex!/getindex for arrays of FSA (e.g. easy acces to single fields) 
	- [ ] access via dimension type (Red -> redchannel)
- [ ] Constructor
	- [x] generic constructor for arbitrary Nvectors
	- [ ] fast constructor for arbitrary types (slow inside other staged functions, slow due to varargs!?)
	- [ ] different constructors for ease of use (zero, eye, from other FSAs, etc...)
- [ ] Functions
	- [x] all kinds of unary/binary operators
	- [x] matrix multiplication (speed issue with generic constructors. Issue with staged function inside staged function?!)
	- [ ] matrix functions (inv, transpose, etc...)
- [ ] FSA Wrapper
	- [x] Abstract Wrapper type (for types that wrap other FSAs)
	- [x] Indexing
	- [ ] Map/Reduce (-> so no other functions yet)
	
usage:
```Julia

immutable RGB{T} <: FixedSizeVector{T, 3}
r::T
g::T
b::T
end

immutable FSMatrix{T <: FixedSizeMatrix} <: FixedSizeWrapper{T} 
data::T
end
```
#### Acknowledgements
[ImmutableArrays](https://github.com/twadleigh/ImmutableArrays.jl) by [twadleigh](https://github.com/twadleigh) was the package that got me going and gave the initial inspirations.
There has been quite a few discussions on [JuliaLang/julia#7568](https://github.com/JuliaLang/julia/pull/7568) shaping the implementation.
Also, [aaalexandrov](https://github.com/aaalexandrov) supplied some code and inspirations.
