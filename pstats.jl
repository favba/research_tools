#!/usr/bin/env julia

using ReadGlobal,PyPlot

d = readcsv("Stats.txt")
xlabel("time")
ylabel(ARGS[1])
plot(d.time,d[Symbol(ARGS[1])])
PyPlot.show()
#readline(stdin)