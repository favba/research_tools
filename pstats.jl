#!/usr/bin/env julia

using ReadGlobal,GR

d = readcsv("Stats.txt")
plot(d.time,d[Symbol(ARGS[1])])
readline(stdin)