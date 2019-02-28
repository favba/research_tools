#!/usr/bin/env julia

using ReadGlobal, UnicodePlots

function main()
    d = readcsv("Stats.txt")
    show(lineplot(d.time,d[Symbol(ARGS[1])],xlabel="time",ylabel=ARGS[1]))
end

main()
