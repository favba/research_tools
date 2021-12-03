#!/usr/bin/env julia

using Statistics, ReadGlobal, Distributions

function stats(file::String)
    if testinput(file)
        field = readfield(file)
        vf = vec(field)
        med = mean(vf)
        mod = mode(round.(vf,digits=4))
        mdevi = mean(x->abs(x-med),vf)
        devi = std(vf,mean=med)
        fmin, imin = findmin(field)
        fmax, imax = findmax(field)

        println("Mean: ",med)
        println("Mode (heuristics): ",mod)
        println("Mean dev: ",mdevi)
        println("Std: ",devi)
        println("Tail size: ",devi/mdevi - 1)
        println("Max: ",fmax," at ",Tuple(imax))
        println("Min: ",fmin," at ",Tuple(imin))
    else
        nb = Base.Filesystem.filesize(file)
        nn = nb÷8
        v = Vector{Float64}(undef,nn)
        read!(file,v)
        med = mean(v)
        mod = mode(round.(v,digits=3))
        mdevi = mean(x->abs(x-med),v)
        devi = std(v,mean=med)
        println("Mean: ",med)
        println("Mode (heuristics): ",mod)
        println("Mean dev: ",mdevi)
        println("Std: ",devi)
        println("Tail size: ",devi/mdevi - 1)
    end
end

stats(ARGS[1])
