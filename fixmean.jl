#!/usr/bin/env julia

using FluidFields, GlobalFileHelper, Statistics

function fixmean!(fab,fa,fb,ma,mb)
    @inbounds for i in eachindex(fab)
        fab[i] = fab[i] - mb*fa[i] - ma*fb[i] + ma*mb
    end
end

function main()
    fab = ScalarField(ARGS[1])
    fa = ScalarField(ARGS[2])
    fb = ScalarField(ARGS[3])
    ma = mean(fa.r)
    mb = mean(fb.r)
    fixmean!(fab.rr,fa.rr,fb.rr,ma,mb)
    write(ARGS[4],fab)
end

main()