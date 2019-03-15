#!/usr/bin/env julia

# Script to compute the energy redistribution tensor
#
# Usage: (Note that the 33 elements should not be inputed)
#         |------------ input -----------||- output -| 
# makeP.jl S11 S12 S13 S22 S23 W12 W13 W23 Pname
using FluidTensors, GlobalFileHelper, LinearAlgebra

function calc(S,W,P)
    Threads.@threads for i in eachindex(P)
        @inbounds begin
            P[i] = Lie(S[i],W[i])
        end
    end
end

function run()
    nx,ny,nz,lx,ly,lz = getdimsize()
    dtp = checkinput.(ARGS[1:8],nx,ny,nz)
    dtypes = getindex.(dtp,1)
    paddeds = getindex.(dtp,2)
    @assert all(dtypes[1] .=== dtypes)
    T = dtypes[1]

    inp = (SymTrTenArray{T} => (ARGS[1:5]...,), AntiSymTenArray{T} => (ARGS[6:8]...,))
    out = (SymTenArray{T} => (ARGS[9] .* ("11","12","13","22","23","33")),)
    doinchunks(calc,input=inp,output=out)
end

run()
