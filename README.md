# FixedSizeArrays

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
