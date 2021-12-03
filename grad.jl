#!/usr/bin/env julia

# Script to compute the gradient of fields
#
#             
# grad.jl input-scalar-file output-name [output-posfix]
# grad.jl input-vector-file1 input-vector-file3 input-vector-file3 output_symmetric-name output-antisymmetric-name [output-posfix]

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

function grad_vec!(sym::SymTenField,asym::VectorField,vec::VectorField)
    isrealspace(vec) && fourier!(vec)
    setfourier!(sym)
    setfourier!(asym)
    Threads.@threads for l in axes(vec,3)
        @inbounds for j in axes(vec,2)
            for i in axes(vec,1)
                ∇ = im*vec.k[i,j,l]
                v = vec[i,j,l]
                sym[i,j,l] = symouter(∇,v)
                asym[i,j,l] = (∇×v)/2
            end
        end
    end
    real!(sym)
    real!(asym)
    return sym, asym
end

function grad_vec(vec::VectorField)
    asym = similar(vec)
    sym = SymTenField(similar(vec.c.x),similar(vec.c.x),similar(vec.c.x),similar(vec.c.x),similar(vec.c.x),similar(vec.c.x))
    return grad_vec!(sym,asym,vec)
end

function grad(s::AbstractString,v::AbstractString,p::AbstractString)
    scalar = ScalarField(s)
    vec = grad_scalar(scalar)
    write(string(v,"1",p),vec.c.x.r)
    write(string(v,"2",p),vec.c.y.r)
    write(string(v,"3",p),vec.c.z.r)
    return nothing
end

function grad(v1::AbstractString,v2::AbstractString,v3::AbstractString,S::AbstractString,W::AbstractString,p::AbstractString)
    vec = VectorField(v1,v2,v3)
    sym,asym = grad_vec(vec)
    write(string(S,"11",p),sym.c.xx.r)
    write(string(S,"12",p),sym.c.xy.r)
    write(string(S,"13",p),sym.c.xz.r)
    write(string(S,"22",p),sym.c.yy.r)
    write(string(S,"23",p),sym.c.yz.r)
    write(string(S,"33",p),sym.c.zz.r)
    
    write(string(W,"12",p),asym.c.z.r)
    write(string(W,"13",p),-asym.c.y.r)
    write(string(W,"23",p),asym.c.x.r)

    return nothing
end

#grad(ARGS...)
function main(args)
    if length(args) == 2
        grad(args[1],args[2],"")
    elseif length(args) == 3
        grad(args[1],args[2],args[3])
    elseif length(args) == 5
        grad(args[1],args[2],args[3],args[4],args[5],"")
    elseif length(args) == 6
        grad(args[1],args[2],args[3],args[4],args[5],args[6])
    else
        error("Wrong number of arguments")
    end
    return 0
end

main(ARGS)
