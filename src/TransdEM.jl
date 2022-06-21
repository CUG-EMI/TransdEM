#-------------------------------------------------------------------------------
#
# A framework for transdimensional Bayesian inversion of electromagnetic data
# in stratified media
#
# (c) Ronghua Peng and Bo Han, China University of Geosciences, Wuhan, 2020-2022.
#
#-------------------------------------------------------------------------------
VERSION == v"1.0.0"

module TransdEM

include("TBUtility/TBUtility.jl")
include("TBStruct/TBStruct.jl")
include("TBFileIO/TBFileIO.jl")
include("TBFwdSolver/TBFwdSolver.jl")
include("TBChain/TBChain.jl")
include("TBAnalysis/TBAnalysis.jl")

end
