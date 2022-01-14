#!/usr/bin/env julia
using ReadGlobal

function sgs!(r::Array{T,N},u::Array{T,N},v::Array{T,N}) where {T<:Number,N}
  Threads.@threads for i::Int in 1:length(u)
    @inbounds r[i] = r[i] - u[i]*v[i]
  end
  return r
end

function sum!(r::Array{T,N},u::Array{T,N},v::Array{T,N}) where {T<:Number,N}
  Threads.@threads for i::Int in 1:length(u)
    @inbounds r[i] = r[i]+u[i]+v[i]
  end
  return r
end

function traceless!(r::Array{T,N},u::Array{T,N}) where {T<:Number,N}
  third::T = 1/3
  Threads.@threads for i::Int in 1:length(u)
    @inbounds r[i] = r[i] - third*u[i]
  end
  return r
end

function main()

place = split(pwd(),"/")
N = place[end][2:end]
Fil = place[end-1][1]

nx, ny, nz, lx, ly, lz = getdimsize()

aux1 = readfield("$(Fil)u1u1_N$N",nx,ny,nz)
aux2 = readfield("$(Fil)u1_N$N",nx,ny,nz)
write("$(Fil)Tr11_N$N",sgs!(aux1,aux2,aux2))

read!("$(Fil)u1u2_N$N",aux1)
aux3 = readfield("$(Fil)u2_N$N",nx,ny,nz)
write("$(Fil)T12_N$N",sgs!(aux1,aux2,aux3))

read!("$(Fil)u1u3_N$N",aux1)
read!("$(Fil)u3_N$N",aux3)
write("$(Fil)T13_N$N",sgs!(aux1,aux2,aux3))

read!("$(Fil)u2u2_N$N",aux1)
read!("$(Fil)u2_N$N",aux2)
write("$(Fil)Tr22_N$N",sgs!(aux1,aux2,aux2))

read!("$(Fil)u2u3_N$N",aux1)
write("$(Fil)T23_N$N",sgs!(aux1,aux2,aux3))

read!("$(Fil)u3u3_N$N",aux1)
write("$(Fil)Tr33_N$N",sgs!(aux1,aux3,aux3))

read!("$(Fil)Tr11_N$N",aux2)
read!("$(Fil)Tr22_N$N",aux3)

write("2k",sum!(aux1,aux2,aux3))

write("$(Fil)T11_N$N",traceless!(aux2,aux1))
write("$(Fil)T22_N$N",traceless!(aux3,aux1))

read!("$(Fil)Tr33_N$N",aux3)
write("$(Fil)T33_N$N",traceless!(aux3,aux1))

return 0
end

main()
