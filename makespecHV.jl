#!/usr/bin/env julia
using ReadGlobal, Statistics, DelimitedFiles, FluidFields

function calculate_specHV(h,v,out,KX,KY,KZ,NX,YRANGE,ZRANGE)
    DKZ = KZ[2]
    MAXDKH = max(KY[2],KX[2])
    fill!(h,0.0)
    fill!(v,0.0)
    @inbounds for l in ZRANGE
        KZ2 = KZ[l]^2
        nk = round(Int, sqrt(KZ2)/DKZ) + 1
        for j in YRANGE
            KY2 = KY[j]^2
            nj = round(Int, sqrt(KY2)/MAXDKH) + 1
            ee = out[1,j,l]
            h[1] += 0.5*ee
            h[nj] += 0.5*ee
            v[nk] += ee
            for i in 2:NX
                k = sqrt(KX[i]*KX[i])

                nx = round(Int, k/MAXDKH) + 1

                ee = 2*out[i,j,l]
                h[nx] += 0.5*ee*MAXDKH
                h[nj] += 0.5*ee*MAXDKH
                v[nk] += ee*DKZ
 
            end
        end
    end
end

function makeHV(file::String,kc::Float64=-1.0)

    NRX,NY,NZ,LX,LY,LZ = getdimsize()

    NX = div(NRX,2)+1

    KX = FluidFields.SRKvec(NRX,LX)#(kxp...,)
    KY = FluidFields.SKvec(NY,LY)#(kyp...,)
    KZ = FluidFields.SKvec(NZ,LZ)#(kzp...,)
    KRZ = FluidFields.SRKvec(NZ,LZ)
    KH = (KX[2] >= KY[2]) ? KX : FluidFields.SRKvec(NY,LY)

    f3D = zeros((NX,NY,NZ))
    read!(file,f3D)

    specH = zeros(Float64,length(KH))
    specV = zeros(Float64,length(KRZ))
    YRANGE = eachindex(KY)
    ZRANGE = eachindex(KZ)

    if kc > 0.0
        IX = findall(x->abs(x)>abs(kc),KX)
        IY = findall(x->abs(x)>abs(kc),KY)
        IZ = findall(x->abs(x)>abs(kc),KZ)

        f3D[IX,:,:] .= 0.0
        f3D[:,IY,:] .= 0.0
        f3D[:,:,IZ] .= 0.0
    end

    calculate_specHV(specH,specV,f3D,KX,KY,KZ,NX,YRANGE,ZRANGE)
    
    sf = split(file,".")
    n = sf[1]
    i = sf[end]
    writedlm("$(n).specH_kc$(kc).$i",zip(KH,specH))
    writedlm("$(n).specV_kc$(kc).$i",zip(KRZ,specV))

end

makeHV(args::Vector{String}) = length(args) == 2 ? makeHV(args[1],parse(Float64,args[2])) : makeHV(args[1])

makeHV(ARGS)