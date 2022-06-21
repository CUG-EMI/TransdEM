#------------------------------------------------------------------------------
# synthetic data for CSEM case
#
#------------------------------------------------------------------------------
# push!(LOAD_PATH, "/yourpath/TransdEM")
#
using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBFileIO
using TransdEM.TBFwdSolver

#------------------------------------------------------------------------------
# 1D synthetic model from Constable et al., 2006.
#
rho    = [1.0, 100.0, 1.0]
hLayer = [1000.0, 100.0]
sigAir   = 1e-12
airLayer = -1000000.
sigSea   = 3.2
seaLayer = 1000.
#
sigma    = vcat(sigAir,sigSea,1 ./rho)
seaDepth = cumsum(vcat(0., hLayer))
depth1D  = vcat(airLayer, 0., seaDepth .+ seaLayer)

# simulated data
dpLen  = 0.
txLoc  = [0. 0. 950. 0. 0.]
nTx    = size(txLoc, 1)
nRx    = 20
rxLoc  = zeros(nRx, 3)
rxLoc[:, 1]  = collect(500:500:10000)
rxLoc[:, 3] .= 1000.0
freqs  = [0.25, 0.5, 0.75]
nFreq  = length(freqs)
dataType = ["ampEx", "phsEx"]
nDt    = length(dataType)

#
tvec = zeros(0)
dataID = ones(Bool,nDt,nRx,nTx,nFreq)
dataID = vec(dataID)
seawater = [seaLayer,sigSea]
emData = CSEMData(seawater,txLoc,rxLoc,0.0,freqs,dataType,dataID,tvec,tvec)

@time predData = compCSEM1dResp(emData, sigma, depth1D)

# add noise
errlvl  = 0.05
nData   = length(predData)
dataErr = zeros(nData)
ampErr  = abs.(predData[1:2:end]) * errlvl
#
errFloor =  1e-15
idx = findall(ampErr .< errFloor) 
ampErr[idx] .= errFloor
dataErr[1:2:end]  = ampErr
dataErr[2:2:end] .= errlvl * 180/pi
obsData = predData + dataErr .* randn(nData)

# output
datafile = "csem1data.dat"
writeCSEMData(datafile, emData, obsData, dataErr)
