#!/usr/bin/env julia

# Script to compute the energy redistribution tensor
#
# Usage: (Note that the 33 elements should not be inputed)
#                  |----------------- input --------------------------||---------- output ---------| 
# makeEnergyDir.jl S11 S12 S13 S22 S23 W12 W13 W23 T11 T12 T13 T22 T23 dir11 dir12 dir13 dir22 dir23
using FluidTensors, GlobalFileHelper, LinearAlgebra

S = SymTrTenArray(ARGS[1:5]...)
W = AntiSymTenArray(ARGS[6:8]...)
T = SymTrTenArray(ARGS[9:13]...)

@. T = traceless(SymTen(T⋅(S+W) + (S-W)⋅T))

write((ARGS[14:18]...,),T)
write(string(SubString(ARGS[18],1,lastindex(ARGS[18])-2),"33"),(-).(T.xx .+ T.yy))
