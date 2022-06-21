export readstartupFile

#-------------------------------------------------------------------------------
function readstartupFile(startupfile::String)

    #
    if isfile(startupfile)
        fid = open(startupfile, "r")
    else
        error("$(startupfile) does not exist, please try again.")
    end

    #
    surveyType = []
    datafile = []
    mclimits = initMCPrior()
    while !eof(fid)
        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while cline[1] == '#' || isempty(cline)
            cline = strip(readline(fid))
        end
        # cline = lowercase(cline)

        # data filename
        if occursin("datafile:", cline)
            tmp = split(cline)
            datafile = string(tmp[end])

        elseif occursin("surveytype:", cline)
            tmp = split(cline)
            surveyType = string(tmp[end])

        # markov chain parameter
        elseif occursin("burninsamples:", cline)
            tmp = split(cline)
            mclimits.burninsamples = parse(Int, tmp[end])

        elseif occursin("totalsamples:", cline)
            tmp = split(cline)
            mclimits.totalsamples = parse(Int, tmp[end])

        # prior parameter
        elseif occursin("numberoflayer:", cline)
            tmp = split(cline)
            mclimits.nlayermin = parse(Int, tmp[end-1])
            mclimits.nlayermax = parse(Int, tmp[end])

        elseif occursin("zcoordinate(m):", cline)
            tmp = split(cline)
            zmin = parse(Float64, tmp[end-2])
            zmax = parse(Float64, tmp[end-1])
            zstd = parse(Float64, tmp[end])
            #
            mclimits.zmin = log10(zmin)
            mclimits.zmax = log10(zmax)
            mclimits.zstd = zstd * (mclimits.zmax - mclimits.zmin)

        elseif occursin("minthickness(m):",cline)
            tmp = split(cline)
            mclimits.hmin = parse(Float64, tmp[end])

        # proposal parameter
        elseif occursin("resistivity:", cline)
            tmp = split(cline)
            rhomin = parse(Float64, tmp[end-2])
            rhomax = parse(Float64, tmp[end-1])
            rhostd = parse(Float64, tmp[end])
            #
            mclimits.rhomin = log10(rhomin)
            mclimits.rhomax = log10(rhomax)
            mclimits.rhostd = rhostd * (mclimits.rhomax - mclimits.rhomin)
            mclimits.mrhostd = rhostd * (mclimits.rhomax - mclimits.rhomin)

        # parameters for post analysis
        elseif occursin("numberofbins:", cline)
            tmp = split(cline)
            mclimits.nzBins = parse(Float64, tmp[end-1])
            mclimits.npBins = parse(Float64, tmp[end])

        elseif occursin("credinterval:", cline)
            tmp = split(cline)
            mclimits.credInterval = parse(Float64, tmp[end])

        else
            @warn("$(cline) is not supported!")

        end

    end # while
    close(fid)

    # read observed data
    if surveyType == "tem"
        println("TEM data is loading ...")
        emData = readTEMData(datafile)
    elseif surveyType == "mt"
        println("MT data is loading ...")
        if occursin(".edi", datafile)
            emData = readMTEDIFile(datafile, 3, 0.05)
        else
            emData = readMTData(datafile)
        end
    elseif surveyType == "csem"
        println("CSEM data is loading ...")
        emData = readCSEMData(datafile)
    end

    return emData, mclimits

end
