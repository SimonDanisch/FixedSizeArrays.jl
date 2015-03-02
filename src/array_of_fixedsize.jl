meths = methodswith(DenseArray, true)
unary = filter(meths) do meth
	length(meth.sig) == 1
end
sinm = filter(unary) do x
string(x.env.code.name) == "sin"
end
println(sinm)