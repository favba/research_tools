#!/usr/bin/env julia

using ReadGlobal, Gaston

d = readcsv("Stats.txt")
plot(d.time,d[Symbol(ARGS[1])],xlabel="time",ylabel=ARGS[1])
Gaston.llplot()
readline(stdin)