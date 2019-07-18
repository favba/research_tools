#!/usr/bin/env julia
using ReadGlobal, Statistics, Glob, DelimitedFiles

function main(R::StepRange)

    nR = length(R)

    filesH = getindex.(split.(glob("*.spec3D.0"),"."),1)

    nx,ny,nz,_,_,_ = getdimsize()

    lh = (div(nx,2)+1)*ny*nz

    auxh = zeros(lh)
    aux2 = zeros(lh)

    for f in filesH
        fill!(auxh,0.0)
        fh = string(f,".spec3D.")
        for r in R
            auxh .+= read!(string(fh,r),aux2)
        end
        auxh ./= nR
        write(string(fh,"avg","_",R.start,"-",R.stop),auxh)
    end

end

main(args::Vector{String}) = main((parse(Int,args[1]):parse(Int,args[2]):parse(Int,args[3])))

main(ARGS)