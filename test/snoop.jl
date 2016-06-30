using SnoopCompile

SnoopCompile.@snoop "fsa_compiles.csv" begin
    include(Pkg.dir("FixedSizeArrays", "test", "runtests.jl"))
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

# IMPORTANT: we must have the module(s) defined for the parcelation
# step, otherwise we will get no precompiles for the Images module
using FixedSizeArrays

data = SnoopCompile.read("fsa_compiles.csv")

blacklist = ["MIME"]

pc = SnoopCompile.format_userimg(data[end:-1:1,2], subst=subst, blacklist=blacklist)
SnoopCompile.write("userimg_fsa.jl", pc)