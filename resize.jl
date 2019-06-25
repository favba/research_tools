#!/usr/bin/env julia
using FluidFields, GlobalFileHelper

function resize(input,u1_big) where {T,L}
    nx,ny,nz = size(input)
    nxb,nyb,nzb = size(u1_big)
  
    nxs = min(nx,nxb)
    nys = min(ny,nyb)
    nzs = min(nz,nzb)
  
    fourier!(input)
    @views begin
        u1_big[1:nxs,1:(div(nys,2)),1:(div(nzs,2))] .= input[1:nxs,1:(div(nys,2)),1:(div(nzs,2))]
        u1_big[1:nxs,(nyb-div(nys,2)):end,1:(div(nzs,2))] .= input[1:nxs,(ny-div(nys,2)):end,1:(div(nzs,2))]
        u1_big[1:nxs,(nyb-div(nys,2)):end,(nzb-div(nzs,2)):end] .= input[1:nxs,(ny-div(nys,2)):end,(nz-div(nzs,2)):end]
        u1_big[1:nxs,1:(div(nys,2)),(nzb-div(nzs,2)):end] .= input[1:nxs,1:div(nys,2),(nz-div(nzs,2)):end]

        u1_big[end,:,:] .= false
    end
    real!(u1_big)
    return u1_big
end
  
function resize(input::String,nnx::Int,nny::Int,nnz::Int)
    nx,ny,nz,lx,ly,lz = getdimsize()
    us = ScalarField(input)
    out = ScalarField{eltype(us.rr)}((nnx,nny,nnz),(lx,ly,lz))
    setfourier!(out)
    resize(us,out)
    write(input*"_resized",out)
end


resize(ARGS[1],parse(Int,ARGS[2]),parse(Int,ARGS[3]),parse(Int,ARGS[4]))