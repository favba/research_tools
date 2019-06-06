#!/usr/bin/env julia

# Script to compute the energy redistribution tensor
#
# Usage: (Note that the 33 elements should not be inputed)
#                  |----------------- input --------------------------||- output -| 
# makeEnergyDir.jl S11 S12 S13 S22 S23 W12 W13 W23 T11 T12 T13 T22 T23   dir pr outposfix
using FluidFields, GlobalFileHelper, LinearAlgebra, FluidTensors

function calc(S,W,T,dir,pr)
    Threads.@threads for i in eachindex(pr)
        @inbounds begin
            R = symmetric(T[i]⋅(S[i]+W[i]))
            pr[i] = tr(R)
            dir[i] = traceless(R)
        end
    end
end

function run()

    S = SymTrTenField(ARGS[1:5]...)
    W = AntiSymTenField(ARGS[6:8]...)
    Tr = SymTrTenField(ARGS[9:13]...)

    calc(S.rr,W.rr,Tr.rr,Tr.rr,W.rr.xy)

    write.(string.(ARGS[14],("11","12","13","22","23"),ARGS[16]),getfield.(Ref(Tr.rr),(:xx,:xy,:xz,:yy,:yz)))
    write(string(ARGS[15],ARGS[16]),W.rr.xy)

    return 0
end

run()
