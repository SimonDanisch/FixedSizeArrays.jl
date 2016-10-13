using SnoopCompile
SnoopCompile.@snoop "fsa_compiles.csv" begin
    include(Pkg.dir("FixedSizeArrays/test/runtests.jl"))
end

using FixedSizeArrays

str = open("fsa_compiles.csv") do io
    str = readstring(io)
    x = replace(str, r"[0-9]+\t\"<toplevel thunk> [\s\S]+?\)::Any\]\"\n", "")
end
open("fsa_compiles.csv", "w") do io
    seekstart(io)
    print(io, str)
end
data = SnoopCompile.read("fsa_compiles.csv")
blacklist = ["MIME"]
pc = SnoopCompile.format_userimg(data[end:-1:1,2], blacklist=blacklist)
SnoopCompile.write(Pkg.dir("FixedSizeArrays", "src", "glv_userimg.jl"), pc)
