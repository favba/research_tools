#!/usr/bin/env julia
using ReadGlobal, Statistics

function main()

    stats = readcsv("Stats.txt")
    statsheaders = getfield(stats,:names)
    statsdata = getfield(stats,:data)
    meanstats = zeros(length(statsdata))

    @inbounds for l in eachindex(statsdata)
        meanstats[l] = mean(view(statsdata[l],202:length(statsdata[l])))
    end

    open("StatsAvg.txt",write=true) do io
        print(io,join(statsheaders,","),"\n",join(meanstats,","))
    end

    scales = readcsv("Scales.txt")
    scalesheaders = getfield(scales,:names)
    scalesdata = getfield(scales,:data)
    meanscales = zeros(length(scalesdata))
    @inbounds for l in eachindex(scalesdata)
        meanscales[l] = mean(view(scalesdata[l],202:length(scalesdata[l])))
    end

    open("ScalesAvg.txt",write=true) do io
        print(io,join(scalesheaders,","),"\n",join(meanscales,","))
    end

    numbers = readcsv("Numbers.txt")
    numbersheaders = getfield(numbers,:names)
    numbersdata = getfield(numbers,:data)
    meannumbers = zeros(length(numbersdata))
    @inbounds for l in eachindex(numbersdata)
        meannumbers[l] = mean(view(numbersdata[l],202:length(numbersdata[l])))
    end

    open("NumbersAvg.txt",write=true) do io
        print(io,join(numbersheaders,","),"\n",join(meannumbers,","))
    end

end

main()