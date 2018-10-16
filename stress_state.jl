#!/usr/bin/env julia

using LinearAlgebra, ReadGlobal, FluidTensors

function stress(S,o)
    Threads.@threads for i in LinearIndices(S)
        @inbounds o[i] = stress_state(S[i])
    end
end

function run()
    inp = (SymTenArray=>(ARGS[1:6]...,),)
    out = (Vector{Float64}=>(ARGS[7],),)
    doinchunks(stress,input=inp,output=out)
end

run()