#!/usr/bin/env julia

# Script to compute the vertical component of the eigenvectors of a Tensor
#
#           input tensor           | output |
# dot.jl T11 T12 T13 T22 T23 T33   direc
#          input tensor               |output |
# dot.jl T11 T12 T13 T22 T23 v1 v2 v3   direc
#
#          input vectors        | output |
# dot.jl   v1  v2  v3 u1 u2 u3    direc



using LinearAlgebra, GlobalFileHelper, FluidTensors, FluidFields

function each_dot(S,v,o1)
    Threads.@threads for i in eachindex(S)
        @inbounds begin
            o1[i] = S[i]⋅v[i]
        end
    end
end


function run(args)
    V = VectorField(args[1:3]...)
    V2 = VectorField(args[4:6]...)
    each_dot(V.rr,V2.rr,V.rr.x)
    write(args[7],V.rr.x)
    return 0
end

run(ARGS)