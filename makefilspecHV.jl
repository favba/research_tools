#!/usr/bin/env julia
using ReadGlobal, Statistics, DelimitedFiles, FluidFields, FieldFilters

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

function makeHVf(file::String,filN::String,fil)

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

    @inbounds for k in axes(f3D,3)
        for j in axes(f3D,2)
            for i in axes(f3D,1)
                f3D[i,j,k] *= fil(KX[i]^2 + KY[j]^2 + KZ[k]^2)^2
            end
        end
    end

    kc = inv(sqrt(fil.Δ²))*π

    IX = findall(x->abs(x)>abs(kc),KX)
    IY = findall(x->abs(x)>abs(kc),KY)
    IZ = findall(x->abs(x)>abs(kc),KZ)

    f3D[IX,:,:] .= 0.0
    f3D[:,IY,:] .= 0.0
    f3D[:,:,IZ] .= 0.0

    calculate_specHV(specH,specV,f3D,KX,KY,KZ,NX,YRANGE,ZRANGE)
    
    sf = split(file,".")
    n = sf[1]
    i = sf[end]
    writedlm("$(n).specH_fil_$(filN).$i",zip(KH,specH))
    writedlm("$(n).specV_fil_$(filN).$i",zip(KRZ,specV))

end

function makeHV(args::Vector{String})
    N = parse(Float64,ARGS[2])
    fil = ARGS[3]

    nx,ny,nz,lx,ly,lz = getdimsize()
    boxdim = N*2*lz*π/nz

    if fil == "G"
        filterkernel = GaussianFilter(boxdim^2)
    elseif fil =="C"
        filterkernel = SpectralCutoffFilter(boxdim^2)
    elseif fil =="B"
        filterkernel = BoxFilter(boxdim^2)
    elseif fil =="E"
        filterkernel = EyinkFilter(boxdim^2)
    elseif fil =="S"
        filterkernel = SpectralBarrier(boxdim^2)
   end

   makeHVf(ARGS[1],fil,filterkernel)
 
end

makeHV(ARGS)