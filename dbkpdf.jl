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

function dbkhist(field,nb::Integer)
    vf = vec(field)
    sort!(vf)
    x = zeros(nb)
    dx = zeros((nb-1))
    y = zeros(nb)
    errx =  zeros(nb)
    r = splitrange(length(vf),nb)
    
    @inbounds for i in 1:nb
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

    for i in 1:(nb-1)
        tv += y[i]*dx[i]
    end

    y ./= tv

    return x, y, errx
end

function main(file::String,nbins::String="100")
    if testinput(file)
        field = readfield(file)
    else
        nb = Base.Filesystem.filesize(file)
        nn = nb÷8
        field = Vector{Float64}(undef,nn)
        read!(file,field)
    end
    x,y,xr = dbkhist(field,parse(Int,nbins))
    writedlm(file*".dbkpdf.txt",zip(x,y,xr))
end

main(ARGS...)