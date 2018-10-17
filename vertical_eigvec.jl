#!/usr/bin/env julia

# Script to compute the vertical component of the eigenvectors of a Tensor
#
#                     input tensor comp          output
# vertical_eigvec.jl T11 T12 T13 T22 T23 T33 v1ez v2ez v3ez

using LinearAlgebra, ReadGlobal, FluidTensors

function veg(S,o1,o2,o3)
    Threads.@threads for i in eachindex(S)
        @inbounds begin
            _,e = eigvec(S[i])
            o1[i] = abs(e[1].z) # most compressing (negative)
            o2[i] = abs(e[2].z) # intermediary
            o3[i] = abs(e[3].z) # most extensive (positive)
        end
    end
end

function run()
    inp = (SymTenArray=>(ARGS[1:6]...,),)
    out = (Vector{Float64}=>(ARGS[7],),Vector{Float64}=>(ARGS[8],),Vector{Float64}=>(ARGS[9],))
    doinchunks(veg,input=inp,output=out)
end

run()