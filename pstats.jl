#!/usr/bin/env julia

using ReadGlobal, UnicodePlots

function main(ARGS::Vector{String})
    d = readcsv("Stats.txt")
    plt = lineplot(d.time,getproperty(d,Symbol(ARGS[1])),xlabel="time",name=ARGS[1])
    for i=2:length(ARGS)
        lineplot!(plt,d.time,getproperty(d,Symbol(ARGS[i])),name=ARGS[i])
    end
    display(plt)
end

main(ARGS)
