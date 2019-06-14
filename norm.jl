#!/usr/bin/env julia

# Script to compute the vertical component of the eigenvectors of a Tensor
#
#           input tensor           | output |
# norm.jl T11 T12 T13 T22 T23 T33   direc
#          input tensor       |output |
# norm.jl T11 T12 T13 T22 T23   direc
#
#          input vector  | output |
# norm.jl   v1  v2  v3     direc



using LinearAlgebra, GlobalFileHelper, FluidTensors, FluidFields

function each_norm(S,o1)
    Threads.@threads for i in eachindex(S)
        @inbounds begin
            o1[i] = norm(S[i])
        end
    end
end


function run(args)
    if length(args) == 7
        T = SymTenField(args[1:6]...)
        each_norm(T.rr,T.rr.xx)
        write(args[7],T.rr.xx)
    elseif length(args) == 6
        T = SymTrTenField(args[1:5]...)
        each_norm(T.rr,T.rr.xx)
        write(args[6],T.rr.xx)
    elseif length(args) == 4
        V = VectorField(args[1:3]...)
        each_norm(V.rr,V.rr.x)
        write(args[4],V.rr.x)
    else
        error("Wrong usage")
    end
    return 0
end

run(ARGS)