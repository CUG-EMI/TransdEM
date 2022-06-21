#------------------------------------------------------------------------------
# script to perform MCMC sampling
#
#------------------------------------------------------------------------------
# push!(LOAD_PATH, "/yourpath/code")
#
using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBFileIO
using TransdEM.TBChain
using TransdEM.TBFwdSolver
using TransdEM.TBAnalysis

#------------------------------------------------------------------------------
printstyled("read datafile and startupfile ...\n", color=:blue)
startup = "startupfile"
(emData,mclimits) = readstartupFile(startup)

#
printstyled("perform MCMC sampling ...\n", color=:blue)
@time (mcArray,mcstatus) = runMCMC(mclimits, emData)

#
printstyled("Output Markov chain ...\n", color=:blue)
@time outputMCMCsamples(mcArray)

printstyled("Extract statistical quantities ...\n", color=:blue)
filestr = "mt"
@time sampleStatistics(mcArray, mclimits, filestr)

#
printstyled("extract model parameters from the ensemble ...\n", color=:blue)
nLayer = 5
extractModelParam(mcArray, mclimits, nLayer)

println("===============================")
#
