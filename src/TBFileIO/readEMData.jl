export readTEMData, readMTData, readCSEMData
export readMTEDIFile, readEDISection
#-------------------------------------------------------------------------------
"""
    `readTEMData(datafile)`

"""
function readTEMData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again.")
    end

    loopxLen  = []; loopyLen = []
    timePoint = []; obsData  = []; dataErr = []
    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while cline[1] == '#' || isempty(cline)
            cline = strip(readline(fid))
        end

        #
        if occursin("LoopWidth(m)", cline)
            tmp = split(cline)
            loopxLen = parse(Float64, tmp[end-1])
            loopyLen = parse(Float64, tmp[end-1])

        elseif occursin("Data Block:", cline)
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            timePoint = zeros(Float64, nData)
            obsData = zeros(Float64, nData)
            dataErr = zeros(nData)
            while nd < nData
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                timePoint[nd] = parse(Float64, cline[1])
                obsData[nd]   = parse(Float64, cline[2])
                dataErr[nd]   = parse(Float64, cline[3])

            end

        end

    end
    close(fid)

    temData = TEMData(loopxLen,loopyLen,timePoint,obsData,dataErr)
    return temData

end


#-------------------------------------------------------------------------------
"""
    readMTData(datafile)

reads MT data file for MT problem.
"""
function readMTData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again.")
    end

    rxLoc = []; freqID = []
    dtID  = []; obsData = []
    dataErr   = []
    dataType  = []
    freqArray = []
    nf  = 0
    nDt = 0
    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while cline[1] == '#' || isempty(cline)
            cline = strip(readline(fid))
        end

        # data format
        if occursin("Format", cline)
            tmp = split(cline)
            format = tmp[2]

        # receiver location
        elseif occursin("Receiver Location", cline)
            tmp = split(cline)
            nr  = parse(Int, tmp[end])
            nd  = 0
            rxLoc = zeros(nr, 3)
            while nd < nr
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:3
                    rxLoc[nd, j] = parse(Float64, cline[j])
                end
            end

        # frequencies
        elseif occursin("Frequencies", cline)
            tmp = split(cline)
            nf  = parse(Int, tmp[end])
            nd  = 0
            freqArray = zeros(nf)
            while nd < nf
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqArray[nd] = parse(Float64, cline)
            end

        # data type
        elseif occursin("DataType", cline)
            tmp = split(cline)
            nDt = parse(Int, tmp[end])
            dataType = Array{String}(undef, nDt)
            nd = 0
            while nd < nDt
                cline = strip(readline(fid))
                nd = nd + 1
                dataType[nd] = cline
            end

        # data block
        elseif occursin("Data Block", cline)
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            rxID  = zeros(Int, nData)
            dtID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            obsData = zeros(Float64, nData)
            dataErr = zeros(nData)

            while nd < nData
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                dtID[nd]   = parse(Int, cline[2])
                obsData[nd] = parse(Float64, cline[3])
                dataErr[nd] = parse(Float64, cline[4])

            end

        end

    end
    close(fid)

    # define dataID
    dataID = zeros(Bool,nDt,nf)
    for k = 1:length(obsData)
        fidx = freqID[k]
        didx = dtID[k]
        dataID[didx,fidx] = true
    end
    dataID = vec(dataID)
    mtData = MTData(freqArray,dataType,dataID,obsData,dataErr)

    return mtData

end


#-------------------------------------------------------------------------------
"""
    readCSEMData(datafile)

read CSEM data information from file.

"""
function readCSEMData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again.")
    end

    txLoc = [];    rxLoc = [];
    txID  = [];    rxID  = [];   freqID = []
    dtID  = [];    obsData = []; dataErr = []
    dataType  = []
    freqArray = []
    phaseCon  = []
    dpLen    = 0.0
    seawater = zeros(0)
    ns = 0; nr = 0; nf = 0; nDt = 0
    while !eof(fid)

        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with #
        while cline[1] == '#' || isempty(cline)
            cline = strip(readline(fid))
        end

        # data format
        if occursin("Format", cline)
            tmp = split(cline)
            format = tmp[2]

        # source type
        elseif occursin("SeaLayer", cline)
            tmp = split(cline)
            seawater = zeros(2)
            seawater[1] = parse(Float64, tmp[end-1])
            seawater[2] = parse(Float64, tmp[end])

        # dipole length, optional
        elseif occursin("Dipole Length", cline)
            tmp = split(cline)
            dpLen = parse(Float64, tmp[end])

        # phase convention
        elseif occursin("Phase Convention", cline)
            tmp = split(cline)
            phaseCon = string(tmp[end])

        # source location
        elseif occursin("Source Location", cline)
            tmp = split(cline)
            ns  = parse(Int, tmp[end])
            nd  = 0
            txLoc = zeros(ns, 5)

            while nd < ns
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:5
                    txLoc[nd,j] = parse(Float64, cline[j])
                end
            end

        # receiver location
        elseif occursin("Receiver Location", cline)
            tmp = split(cline)
            nr  = parse(Int, tmp[end])
            nd  = 0
            rxLoc = zeros(nr, 3)
            while nd < nr
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:3
                    rxLoc[nd,j] = parse(Float64, cline[j])
                end
            end

        # frequencies
        elseif occursin("Frequencies", cline)
            tmp = split(cline)
            nf  = parse(Int, tmp[end])
            nd  = 0
            freqArray = zeros(nf)
            while nd < nf
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqArray[nd] = parse(Float64, cline)
            end

        # data type
        elseif occursin("DataType", cline)
            tmp = split(cline)
            nDt = parse(Int, tmp[end])
            dataType = Array{String}(undef, nDt)
            nd = 0
            while nd < nDt
                cline = strip(readline(fid))
                nd = nd + 1
                dataType[nd] = cline
            end

        # data block
        elseif occursin("Data Block", cline)
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            txID  = zeros(Int, nData)
            rxID  = zeros(Int, nData)
            dtID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            obsData = zeros(Float64, nData)
            dataErr = zeros(nData)

            while nd < nData
                cline = strip(readline(fid))
                while cline[1] == '#' || isempty(cline)
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                txID[nd]   = parse(Int, cline[2])
                rxID[nd]   = parse(Int, cline[3])
                dtID[nd]   = parse(Int, cline[4])
                obsData[nd] = parse(Float64, cline[5])
                dataErr[nd] = parse(Float64, cline[6])
            end

        end

    end
    close(fid)

    # define dataID
    dataID = zeros(Bool,nDt,nr,ns,nf)
    for k = 1:length(obsData)
        tidx = txID[k]
        ridx = rxID[k]
        fidx = freqID[k]
        didx = dtID[k]
        dataID[didx,ridx,tidx,fidx] = true
    end
    dataID   = vec(dataID)
    csemData = CSEMData(seawater, txLoc, rxLoc, dpLen, freqArray, dataType,
                        dataID, obsData, dataErr)

    return csemData

end


#-------------------------------------------------------------------------------
"""
    `readMTEDIFile(edifile, dataType)`

 read MT impedance from edi file

"""
function readMTEDIFile(edifile::String, dataType::Int=1, errLvl=0.05)

    #
    if isfile(edifile)
        fid = open(edifile, "r")
    else
        error("$(edifile) does not exist, please try again.")
    end

    #
    freqs = readEDISection(fid, ">FREQ")
    if dataType == 1
        zxyr  = readEDISection(fid, ">ZXYR")
        zxyi  = readEDISection(fid, ">ZXYI")
        zxyE  = readEDISection(fid, ">ZXY.VAR")
        obsData = hcat(zxyr, zxyi)
        obsData = vec(copy(transpose(obsData)))

    elseif dataType == 2
        zyxr  = readEDISection(fid, ">ZYXR")
        zyxi  = readEDISection(fid, ">ZYXI")
        zyxE  = readEDISection(fid, ">ZYX.VAR")
        obsData = hcat(zyxr, zyxi)
        obsData = vec(copy(transpose(obsData)))

    elseif dataType == 3
        zxxr  = readEDISection(fid, ">ZXXR")
        zxxi  = readEDISection(fid, ">ZXXI")
        zxxE  = readEDISection(fid, ">ZXX.VAR")
        zxx   = zxxr + 1im * zxxi

        zxyr  = readEDISection(fid, ">ZXYR")
        zxyi  = readEDISection(fid, ">ZXYI")
        zxyE  = readEDISection(fid, ">ZXY.VAR")
        zxy   = zxyr + 1im * zxyi

        zyxr  = readEDISection(fid, ">ZYXR")
        zyxi  = readEDISection(fid, ">ZYXI")
        zyxE  = readEDISection(fid, ">ZYX.VAR")
        zyx   = zyxr + 1im * zyxi

        zyyr  = readEDISection(fid, ">ZYYR")
        zyyi  = readEDISection(fid, ">ZYYI")
        zyyE  = readEDISection(fid, ">ZYY.VAR")
        zyy   = zyyr + 1im * zyyi
        zdet  = sqrt.(zxx .* zyy - zxy .* zyx)
        zdetr = real.(zdet)
        zdeti = imag.(zdet)
        obsData = hcat(zdetr, zdeti)
        obsData = vec(copy(transpose(obsData)))

    end

    # convert impedance data from mv/km/nT to Ω
    factor = 796
    obsData /= factor
    dataErr  = abs.(obsData) * errLvl
    datastr = ["realZ", "imagZ"]
    nf  = length(freqs)
    nDt = length(datastr)
    dataID = ones(Bool,nDt*nf)
    mtData = MTData(freqs,datastr,dataID,obsData,dataErr)

end


#-------------------------------------------------------------------------------
function readEDISection(fid::IOStream, token::String)

    data = []
    while !eof(fid)
        cline = strip(readline(fid))
        if occursin(token, cline)
            tmp  = split(cline, "//")
            nd   = parse(Int, tmp[end])
            data = zeros(nd)
            k = 0
            while k < nd
                cline = strip(readline(fid))
                cline = split(cline)
                num   = length(cline)
                for i = 1:num
                    k = k + 1
                    data[k] = parse(Float64, cline[i])
                end
            end
            break
        end
    end

    return data

end
