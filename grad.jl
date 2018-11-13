#!/usr/bin/env julia

# Script to compute the gradient of fields
#
#             
# grad.jl input-scalar-file output-name
# grad.jl input-vector-file1 input-vector-file3 input-vector-file3 output_symmetric-name output-antisymmetric-name

using FluidFields, GlobalFileHelper, FluidTensors, LinearAlgebra

function grad_scalar!(vec::VectorField,scalar::ScalarField)
    isrealspace(scalar) && fourier!(scalar)
    setfourier!(vec)
    Threads.@threads for l in axes(scalar,3)
        @inbounds for j in axes(scalar,2)
            for i in axes(scalar,1)
                ∇ = im*scalar.k[i,j,l]
                vec[i,j,l] = ∇*scalar[i,j,l]
            end
        end
    end
    real!(vec)
    return vec
end

function grad_scalar(scalar::ScalarField)
    vec = VectorField(similar(scalar),similar(scalar),similar(scalar))
    return grad_scalar!(vec,scalar)
end

function grad_vec!(sym::SymTrTenField,asym::VectorField,vec::VectorField)
    isrealspace(vec) && fourier!(vec)
    setfourier!(sym)
    setfourier!(asym)
    Threads.@threads for l in axes(vec,3)
        @inbounds for j in axes(vec,2)
            for i in axes(vec,1)
                ∇ = im*vec.k[i,j,l]
                v = vec[i,j,l]
                sym[i,j,l] = symouter(∇,v)
                asym[i,j,l] = 0.5*∇×v
            end
        end
    end
    real!(sym)
    real!(asym)
    return sym, asym
end

function grad_vec(vec::VectorField)
    asym = similar(vec)
    sym = SymTrTenField(similar(vec.c.x),similar(vec.c.x),similar(vec.c.x),similar(vec.c.x),similar(vec.c.x))
    return grad_vec!(sym,asym,vec)
end

function grad(s::AbstractString,v::AbstractString)
    scalar = ScalarField(s)
    vec = grad_scalar(scalar)
    write(string(v,"1"),vec.c.x.r)
    write(string(v,"2"),vec.c.y.r)
    write(string(v,"3"),vec.c.z.r)
    return nothing
end

function grad(v1::AbstractString,v2::AbstractString,v3::AbstractString,S::AbstractString,W::AbstractString)
    vec = VectorField{Float64}(v1,v2,v3)
    sym,asym = grad_vec(vec)
    write(string(S,"11"),sym.c.xx.r)
    write(string(S,"12"),sym.c.xy.r)
    write(string(S,"13"),sym.c.xz.r)
    write(string(S,"22"),sym.c.yy.r)
    write(string(S,"23"),sym.c.yz.r)
    
    write(string(W,"12"),asym.c.z.r)
    write(string(W,"13"),-asym.c.y.r)
    write(string(W,"23"),asym.c.x.r)

    return nothing
end

grad(ARGS...)
