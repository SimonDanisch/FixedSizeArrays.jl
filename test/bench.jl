using BenchmarkTools, FixedSizeArrays
immutable Normal{N, T} <: FixedVector{N, T}
    values::NTuple{N, T}
end

function test()
    result = []
    sizehint!(result, 2*2*2*3*3)
    for N=(1,10)
        for VT=[Point, Vec], VT2=[Normal, Vec], ET=[Float32, Int, UInt], ET2=[Float64, UInt, Float32]
            rand_range  = ET(1):ET(10)
            rand_range2 = ET2(1):ET2(10)
            rn = rand(rand_range, N)
            v0 = VT(rn)
            push!(result, v0)
            
            # parse constructor:
            VT{N, ET}(map(string, rn))
            # multi constructor
            v1 = VT{N, ET}(rn...)
            length(typeof(v1))
            eltype(typeof(v1))
            push!(result, v1)

            for i=1:N
                v1[i]
            end
            # from other FSA without parameters
            v2 = VT2(v1)
            push!(result, v2)

            for i=1:N
                v2[i]
            end
            # from other FSA with parameters
            for i=1:N
                v3 = VT2{i, ET2}(v1)
                for i=1:i
                    v3[i]
                end
            end
            push!(result, v3)

            # from single
            r  = rand(rand_range)
            r2 = rand(rand_range2)
            v1 = VT{N, ET}(r)
            v2 = VT{N, ET2}(r)
            v3 = VT{N, ET}(r2)
            v4 = VT{N, ET2}(r2)
            push!(result, v4)
            push!(result, v1)


            for i=1:N
                v1[i]
                v2[i]
                v3[i]
                v4[i]
            end
            x = VT{N, ET}[VT{N, ET}(1) for i=1:10]
            x1 = VT2{N, ET}[VT{N, ET}(1) for i=1:10]
            x2 = map(VT2, x)
            x3 = map(VT, x2)
            push!(result, x)
            push!(result, x1)
            push!(result, x2)
            push!(result, x3)
            push!(result, VT{N}(ntuple(x->ET(1), N)...))
            push!(result,  VT{N}(ET(1)))

            # Construction with only N, issue #56
           
            
        end
    end
    result
end