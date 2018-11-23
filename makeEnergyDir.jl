#!/usr/bin/env julia

# Script to compute the energy redistribution tensor
#
# Usage: (Note that the 33 elements should not be inputed)
#                  |----------------- input --------------------------||---------- output -------------| 
# makeEnergyDir.jl S11 S12 S13 S22 S23 W12 W13 W23 T11 T12 T13 T22 T23 dir pr
using FluidTensors, GlobalFileHelper, LinearAlgebra

function calc(S,W,T,dir,pr)
    Threads.@threads for i in eachindex(pr)
        @inbounds begin
            R = SymTen(T[i]⋅(S[i]+W[i]) + (S[i]-W[i])⋅T[i]) 
            pr[i] = tr(R)/2
            dir[i] = traceless(R)
        end
    end
end

function run()
    nx,ny,nz,lx,ly,lz = getdimsize()
    dtp = checkinput.(ARGS[1:13],nx,ny,nz)
    dtypes = getindex.(dtp,1)
    paddeds = getindex.(dtp,2)
    @assert all(dtypes[1] .=== dtypes)
    T = dtypes[1]

    inp = (SymTrTenArray{T} => (ARGS[1:5]...,), AntiSymTenArray{T} => (ARGS[6:8]...,), SymTrTenArray{T} => (ARGS[9:13]...,))
    out = (SymTenArray{T} => (ARGS[14] .* ("11","12","13","22","23","33")), Vector{T} =>(ARGS[15],))
    doinchunks(calc,input=inp,output=out)
end

#S = SymTrTenArray(ARGS[1:5]...)
#W = AntiSymTenArray(ARGS[6:8]...)
#T = SymTrTenArray(ARGS[9:13]...)

#@. T = traceless(SymTen(T⋅(S+W) + (S-W)⋅T))

#write((ARGS[14:18]...,),T)
#write(string(SubString(ARGS[18],1,lastindex(ARGS[18])-2),"33"),(-).(T.xx .+ T.yy))

run()
