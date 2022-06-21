#-------------------------------------------------------------------------------
#
# module `TBStruct` defines data structure and routines for transdimensional
# MCMC sampling.
#
#
#-------------------------------------------------------------------------------
module TBStruct

using  TransdEM.TBUtility
export MCParameter, MCPrior
export MCStatus, MChainStep
export MCDataFit
export MChainArray
export initMCPrior, initMCParameter
export initMCDataFit, initMChainArray
export initModelParameter

#-------------------------------------------------------------------------------
"""
    Struct `MCParameter` defines parameter structure for MCMC chains.

"""
mutable struct MCParameter{Ti<:Int, Tv<:Float64}

    #
    nLayer::Ti               # current layer number
    layeridx::Ti             # selected layer index to perturb
    zNode::Vector{Tv}        # interface depth
    rho::Vector{Tv}          # layer resistivity
    inBound::Bool            # indicator for inside or outside prior bounds

end # MCParameter


#-------------------------------------------------------------------------------
"""
    Struct `MCPrior` defines prior bounds of variables.

"""
mutable struct MCPrior{Ti<:Int, Tv<:Float64}

    #
    burninsamples::Ti
    totalsamples::Ti

    # min/max number of layers
    nlayermin::Ti
    nlayermax::Ti

    # lower/upper bounds for interface depth
    zmin::Tv
    zmax::Tv

    # minimum layer thickness
    hmin::Tv

    # lower/upper bounds for layer resistivity
    rhomin::Tv
    rhomax::Tv

    # standard derivation associated with variables
    zstd::Tv
    rhostd::Tv
    mrhostd::Tv

    # parameters for post analysis
    nzBins::Ti
    npBins::Ti
    credInterval::Tv

end # RJPrior


#-------------------------------------------------------------------------------
"""
    Struct `MCStatus` stores proposal statistics of MCMC chain.

"""
mutable struct MCStatus{T<:Int}

    #
    accepted::Bool
    acceptstats::Vector{T}   # acceptance statistics
    rejectstats::Vector{T}   # rejection statistics

end # RJStatus


#-------------------------------------------------------------------------------
"""
    Struct `MChainStep` records the movement of MCMC chain at every step.

"""
mutable struct MChainStep{T<:Bool}

    isBirth::T
    isDeath::T
    isMove::T
    isPerturb::T

end


#-------------------------------------------------------------------------------
"""
    Struct `MCDataFit` stores data from current and proposed models of MCMC
chain at every step.

"""
mutable struct MCDataFit{T<:Real}

    # current model
    currPredData::Vector{T}
    currLikelihood::T
    currMisfit::T

    # proposed model
    propPredData::Vector{T}
    propLikelihood::T
    propMisfit::T

    # data fitting
    dataResiduals::Vector{T}

end


#-------------------------------------------------------------------------------
"""
    Struct `MChainArray` stores status of MCMC chain at every step.

"""
mutable struct MChainArray{Ti<:Int, Tv<:Real}

    # current sample number
    nsample::Ti

    # starting model
    startModel::MCParameter
    chainstep::Array{Ti, 2}
    chainvalue::Array{Tv, 2}
    chainindice::Vector{Ti}
    chainnlayer::Vector{Ti}
    chainMisfit::Array{Tv}
    residuals::Array{Tv}

end


#-------------------------------------------------------------------------------
"""
    `initMCParameter()` initializes parameter struct for MCMC chain.

"""
function initMCParameter()

    #
    itmp = zero(0)
    pvec = zeros(0)
    mcParam = MCParameter(itmp, itmp, pvec, pvec, false)
    return mcParam

end # initMCParameter


#-------------------------------------------------------------------------------
"""
    `initMCDataFit()` initializes datafit struct for MCMC chain.

"""
function initMCDataFit()

    #
    it = 0.0
    iv = zeros(0)
    mcDatafit = MCDataFit(iv, it, it, iv, it, it, iv)

    return mcDatafit

end # initMCDataFit

#------------------------------------------------------------------------------
"""
    `initMCPrior()`

"""
function initMCPrior()

    #
    dtmp = 0.05
    mclimits = MCPrior(100, 1000, 2, 10, 0., 5.0, 0., 1.0, 5.0, dtmp, dtmp,dtmp,
                       400, 400, 0.9)

    return mclimits

end # initMCPrior


#-------------------------------------------------------------------------------
"""
    `initMChainArray(tblimits::MCPrior)` initializes struct `MChainArray`.

"""
function initMChainArray(mclimits::MCPrior)

    #
    currsample = 0
    startModel = initMCParameter()
    #
    nsample     = mclimits.totalsamples
    chainstep   = zeros(Int, nsample, 2)
    chainvalue  = zeros(Float64, nsample, 2)
    chainindice = zeros(Int, nsample)
    chainncell  = zeros(Int, nsample)
    chainMisfit = zeros(Float64, nsample)
    residuals   = zeros(Float64, 0)

    mcArray = MChainArray(currsample, startModel, chainstep, chainvalue,
                  chainindice, chainncell, chainMisfit, residuals)

    return mcArray

end # initMChainArray


#------------------------------------------------------------------------------
"""
    `initModelParameter(mclimits)`

initialize model parameter.

"""
function initModelParameter(mclimits::MCPrior)

    #
    mcParam = initMCParameter()

    # set number of layer
    currLayerNumber = 3#unirandInteger(mclimits.nlayermin, mclimits.nlayermax)
    mcParam.nLayer  = currLayerNumber

    # set the locations of interface depth randomly
    mcParam.zNode = zeros(Float64, mclimits.nlayermax)
    mcParam.rho   = zeros(Float64, mclimits.nlayermax)
    for i = 1:currLayerNumber
        mcParam.zNode[i] = unirandDouble(mclimits.zmin, mclimits.zmax)

        # set the layer resistivity randomly
        mcParam.rho[i] = unirandDouble(mclimits.rhomin, mclimits.rhomax)

    end

    return mcParam

end


end # TBStruct
