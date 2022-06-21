#------------------------------------------------------------------------------
# synthetic data for TEM case
#
#------------------------------------------------------------------------------
# push!(LOAD_PATH, "/yourpath/TransdEM")
#
using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBFileIO
using TransdEM.TBFwdSolver

#------------------------------------------------------------------------------
# 1D synthetic model 
#
rho    = [20., 1e4, 10., 100.0]
hLayer = [10, 10, 10, 1e3]
depth1D = vcat(0., cumsum(hLayer))

# simulated data
loopxLen = 20.0
loopyLen = 20.0

strtime = -6
endtime = -3
ntime   = (endtime - strtime)*10 + 1
timePoint = 10 .^ LinRange(strtime, endtime, ntime)

#
@time predData = compTEM1dResponse(loopxLen,loopyLen,timePoint,rho,depth1D)
errlvl   = 0.05
errFloor = 1e-9
dataErr  = predData * errlvl
obsData  = predData + dataErr .* randn(ntime) .+ errFloor
dataErr  = dataErr .+ errFloor

# output
temData = TEMData(loopxLen,loopyLen,timePoint,obsData,dataErr)
datafile = "tem1data.dat"
writeTEMData(datafile,temData,obsData,dataErr)
