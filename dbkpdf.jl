#!/usr/bin/env julia

using ReadGlobal, DelimitedFiles

function splitrange(lr,nt)
    a = UnitRange{Int}[]
    sizehint!(a,nt)
    n = lr÷nt
    r = lr%nt
    stop = 0
    init = 1
    for i=1:r
        stop=init+n
        push!(a,init:stop)
        init = stop+1
    end
    for i=1:(nt-r)
        stop=init+n-1
        push!(a,init:stop)
        init = stop+1
    end
    return a
end

function dbkhist(field)
    vf = vec(field)
    sort!(vf)
    x = zeros(100)
    dx = zeros(99)
    y = zeros(100)
    errx =  zeros(100)
    r = splitrange(length(vf),100)
    
    @inbounds for i in 1:100
        elr = r[i]
        ms = zero(eltype(vf))

        for j in elr
            ms += vf[j]
        end

        le = length(elr)
        ms /= le
        ds = vf[elr[end]] - vf[elr[1]]

        x[i] = ms
        y[i] =1/ds

        md = zero(eltype(vf))

        for j in elr
            md += abs(ms-vf[j])
        end

        md /= le
        errx[i] = md
    end

    dx .= @views x[2:end] - x[1:end-1]

    tv = zero(eltype(vf))

    for i in 1:99
        tv += y[i]*dx[i]
    end

    y ./= tv

    return x, y, errx
end

function main(file::String)
    if testinput(file)
        field = readfield(file)
    else
        nb = Base.Filesystem.filesize(file)
        nn = nb÷8
        field = Vector{Float64}(undef,nn)
        read!(file,field)
    end
    x,y,xr = dbkhist(field)
    writedlm(file*".dbkpdf.txt",zip(x,y,xr))
end

main(ARGS[1])