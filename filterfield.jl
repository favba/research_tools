#!/usr/bin/env julia

using FluidFields, GlobalFileHelper, FieldFilters

function main(N::Int,fil::String,inputfiles::NTuple{NT,AbstractString}) where {NT}

    nx,ny,nz,lx,ly,lz = getdimsize()
    boxdim = N*2*lz*π/nz

    if fil == "G"
        pathoutput = "./Filtered_Fields/Gaussian/N$N"
        filterkernel = GaussianFilter(boxdim^2)
    elseif fil =="C"
        pathoutput = "./Filtered_Fields/CutOff/N$N"
        filterkernel = SpectralCutoffFilter(boxdim^2)
    elseif fil =="B"
        pathoutput = "./Filtered_Fields/Box/N$N"
        filterkernel = BoxFilter(boxdim^2)
    elseif fil =="E"
        pathoutput = "./Filtered_Fields/Eyink/N$N"
        filterkernel = EyinkFilter(boxdim^2)
    end
    
    field = ScalarField(inputfiles[1])

    @info("Filtering File $(inputfiles[1])")
    fourier!(field)
    filterfield!(field, filterkernel)
    real!(field)
    @info("Done")


    isdir(pathoutput) || mkpath(pathoutput)
    outfilename = "$(pathoutput)/$(inputfiles[1])"
    @info("Saving filtered field in $outfilename")
    write(outfilename,field.r)
 
    N > 1 && for inputfile in inputfiles[2:end]

        readfield!(inputfile,field.rr)

        @info("Filtering File $inputfile")
        fourier!(field)
        filterfield!(field, filterkernel)
        real!(field)
        @info("Done")


        outfilename = "$(pathoutput)/$inputfile"
        @info("Saving filtered field in $outfilename")
        write(outfilename,field.r)
    end
end

main(parse(Int,ARGS[1]),ARGS[2],(ARGS[3:end]...,))