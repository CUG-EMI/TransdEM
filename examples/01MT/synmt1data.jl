#------------------------------------------------------------------------------
# synthetic data for MT case
#
#------------------------------------------------------------------------------
using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBFileIO
using TransdEM.TBFwdSolver

#------------------------------------------------------------------------------
# 1D synthetic model adapted from Guo et al., 2011.
#
sigma  = [0.004, 0.04, 0.01, 0.1, 0.04]
hLayer = [600, 1400, 4e3, 4e3, 1e4]
depth1D = vcat(0., cumsum(hLayer))

# simulated data
period = 10 .^ LinRange(log10(1/320), log10(80), 40)
freqs  = 1 ./ period
zimp   = compMT1DImpedance(freqs, sigma, depth1D)

# add 5% noise
nfreq = length(freqs)
predData = zeros(2*nfreq)
predData[1:2:end] = real.(zimp)
predData[2:2:end] = imag.(zimp)
errlvl  = 0.05
dataErr = abs.(predData) * errlvl
obsData = predData + dataErr .* randn(2*nfreq)

# output
dataType = ["realZxy", "imagZxy"]
dataID   = ones(Bool,2,nfreq)
dataID   = vec(dataID)
mtData   = MTData(freqs,dataType,dataID,obsData,dataErr)
datafile = "mt1data.dat"
writeMTData(datafile,mtData,obsData,dataErr)
