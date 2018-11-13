#!/usr/bin/env julia

using LinearAlgebra, ReadGlobal, FluidTensors

function stress(S,o)
    Threads.@threads for i in eachindex(S)
        @inbounds o[i] = stress_state(S[i])
    end
end

function run()
    if length(ARGS) == 7
        inp = (SymTenArray=>(ARGS[1:6]...,),)
        out = (Vector{Float64}=>(ARGS[7],),)
    elseif length(ARGS) == 6
        inp = (SymTrTenArray=>(ARGS[1:5]...,),)
        out = (Vector{Float64}=>(ARGS[6],),)
    end
    doinchunks(stress,input=inp,output=out)
end

run()