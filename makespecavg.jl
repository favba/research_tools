#!/usr/bin/env julia
using ReadGlobal, Statistics, Glob, DelimitedFiles

function main(n1,n2,n3)

    #R = 5000:50:20000
    R = n1:n2:n3
    nR = length(R)

    filesH = getindex.(split.(glob("*.specH.0"),"."),1)

    kh = vec(readdlm(filesH[1]*".specH.0")[:,1])
    lh = length(kh)

    auxh = zeros(lh)

    for f in filesH
        fill!(auxh,0.0)
        fh = string(f,".specH.")
        for r in R
            auxh .+= readdlm(string(fh,r),Float64,dims=(lh,2))[:,2]
        end
        auxh ./= nR
        writedlm(string(fh,"avg_",n1,"-",n3),zip(kh,auxh))
    end

    filesV = getindex.(split.(glob("*.specV.0"),"."),1)
    kz = vec(readdlm(filesV[1]*".specV.0")[:,1])
    lz = length(kz)
    auxv = zeros(lz)

    for f in filesV
        fill!(auxv,0.0)
        fv = string(f,".specV.")
        for r in R
            auxv .+= readdlm(string(fv,r),Float64,dims=(lz,2))[:,2]
        end
        auxv ./= nR
        writedlm(string(fv,"avg_",n1,"-",n3),zip(kz,auxv))
    end

end

main(parse(Int,ARGS[1]),parse(Int,ARGS[2]),parse(Int,ARGS[3]))