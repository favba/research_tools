#!/usr/bin/env julia

# Script to compute the local Richardson number of fields
#
#             
# localRichardson.jl S13 S23 W13 W23 rho
using FluidFields, GlobalFileHelper

function dz!(out::AbstractArray{<:Complex,3},f::AbstractArray{<:Complex,3}) 
    isrealspace(f) && fourier!(f)
    isrealspace(out) && setfourier!(out)
    kz = f.kz
    Threads.@threads for k in axes(f,3)
        KZ = kz[k]
        for j in axes(f,2)
            @inbounds @simd for i in axes(f,1)
                out[i,j,k] = f[i,j,k]*im*KZ
            end
        end
    end
end

function calc!(out,s13,s23,w13,w23,dρdz,grho0,mdρdz)
    isrealspace(s13) || real!(s13)
    isrealspace(s23) || real!(s23)
    isrealspace(w13) || real!(w13)
    isrealspace(w23) || real!(w23)
    isrealspace(dρdz) || real!(dρdz)
    isrealspace(out) || setreal!(out)

    s13r = s13.r
    s23r = s23.r
    w23r = w23.r
    w13r = w13.r
    dρdzr = dρdz.r
    outr = out.r

    Threads.@threads for k in axes(s13r,3)
        for j in axes(s13r,2)
            @inbounds @simd for i in axes(s13r,1)
                outr[i,j,k] = max(-1.0, min(1.0,(-grho0*(dρdzr[i,j,k] + mdρdz))/((s13r[i,j,k] + w13r[i,j,k])^2 + (s23r[i,j,k] + w23r[i,j,k])^2 + eps())))
            end
        end
    end

end

function main(ARGS::Vector{String})
    s13 = ScalarField(ARGS[1])
    s23 = ScalarField(ARGS[2])
    w13 = ScalarField(ARGS[3])
    w23 = ScalarField(ARGS[4])
    ρ = ScalarField(ARGS[5])
    d = readglobal()
    g = parse(Float64,d[:zAcceleration])
    ρ0 = parse(Float64,d[:referenceDensity])
    dρdz = parse(Float64,d[:densityGradient])
    dz!(ρ,ρ)
    calc!(ρ,s13,s23,w13,w23,ρ,g/ρ0,dρdz)
    out = length(ARGS) == 6 ? ARGS[6] : "localRichardson"
    write(out,ρ)
end

main(ARGS)
