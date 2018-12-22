#!/usr/bin/env julia

# Script to compute the nonlinearmodel decomposition
#
#                     input         
# vertical_eigvec.jl T11 T12 T13 T22 T23 S11 S12 S13 S22 S23 W12 W13 W23

using LinearAlgebra, ReadGlobal, FluidTensors

function nonlinearmodel(T,S,W,alpha,ratio,model,alpha2,ratio2,model2)
    Threads.@threads for i in eachindex(S)
        @inbounds @fastmath begin
            ∇u = S[i]+W[i]
            nmT = traceless(symmetric(dot(∇u,∇u')))
            lalpha,Tm = prop_decomp(T[i],nmT)
            ratio[i] = (nmT:nmT)/(T[i]:T[i])
            alpha[i] = lalpha
            model[i] = nmT

            nmT2 = traceless(symmetric(dot(∇u',∇u)))
            lalpha2,Tm2 = prop_decomp(T[i],nmT2)
            ratio2[i] = (nmT2:nmT2)/(T[i]:T[i])
            alpha2[i] = lalpha2
            model2[i] = nmT2
       end
    end
end

function run()
    isdir("Model_nonlinear") || mkdir("Model_nonlinear")
    isdir("Model_nonlinear2") || mkdir("Model_nonlinear2")
    inp = (SymTrTenArray=>(ARGS[1:5]...,),SymTrTenArray=>(ARGS[6:10]...,),AntiSymTenArray=>(ARGS[11:13]...,))
    out = (Vector{Float64}=>("Model_nonlinear/alpha",),Vector{Float64}=>("Model_nonlinear/ratio",),SymTenArray{Float64}=>("Model_nonlinear/" .* ("Tm11","Tm12","Tm13","Tm22","Tm23","Tm33")), 
    Vector{Float64}=>("Model_nonlinear2/alpha",),Vector{Float64}=>("Model_nonlinear2/ratio",),SymTenArray{Float64}=>("Model_nonlinear2/" .* ("Tm11","Tm12","Tm13","Tm22","Tm23","Tm33")), 
    )
    doinchunks(nonlinearmodel,input=inp,output=out)
end

run()