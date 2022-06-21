#-------------------------------------------------------------------------------
#
# module `TBChain` defines routines to perform transdimensional MCMC sampling.
#
# (c) Ronghua Peng and Bo Han, China University of Geosciences, Wuhan, 2020-2022.
#
#-------------------------------------------------------------------------------
module TBChain

using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBFileIO
using TransdEM.TBFwdSolver

#
using LinearAlgebra, Printf, Random
using Distributed, DistributedArrays

export mcBirth!, mcDeath!, mcMoveLocation!, mcPerturbProperty!
export checkPropertyBounds
export duplicateMCParameter, updateMCParameter!
export selectChainStep!

#-------------------------------------------------------------------------------
"""
    `duplicateMCParameter(mcParam)`

"""
function duplicateMCParameter(mcParam::MCParameter)

    #
    nLayer   = mcParam.nLayer
    layeridx = mcParam.layeridx
    zNode    = copy(mcParam.zNode)
    rho      = copy(mcParam.rho)
    inBound  = mcParam.inBound
    mcParamNew = MCParameter(nLayer, layeridx, zNode, rho, inBound)

    return mcParamNew

end


#-------------------------------------------------------------------------------
"""
    `updateMCParameter(mcParam, mcParamNew)`

"""
function updateMCParameter!(mcParam::MCParameter, mcParamNew::MCParameter)

    #
    mcParamNew.nLayer   = mcParam.nLayer
    mcParamNew.layeridx = mcParam.layeridx
    mcParamNew.zNode    = copy(mcParam.zNode)
    mcParamNew.rho      = copy(mcParam.rho)
    mcParamNew.inBound  = mcParam.inBound

    return mcParamNew

end


#-------------------------------------------------------------------------------
"""
    `mcBirth!(mcParam, mcParamProposed, mclimits)`

add a new layer to current model.

"""
function mcBirth!(mcParam::MCParameter,
                  mcParamProposed::MCParameter,
                  mclimits::MCPrior)

    # check if current layer number is within the prior bound
    currLayerNumber = mcParam.nLayer
    if currLayerNumber == mclimits.nlayermax
        mcParamProposed.inBound = false
        return
    end

    #
    propLayerNumber  = currLayerNumber + 1
    mcParamProposed.nLayer = propLayerNumber

    # get random location for the new layer interface
    zLoc = 0.0
    iter = 0
    while iter < 1000
        iter += 1
        zLoc = unirandDouble(mclimits.zmin, mclimits.zmax)
        mcParamProposed.zNode[propLayerNumber] = zLoc

        # check if meet the mimimum layer thickness constraint
        if mclimits.hmin > 0.
            hmin = mclimits.hmin
            zNode = copy(mcParamProposed.zNode[1:propLayerNumber])
            sort!(zNode)
            zNode = 10 .^ zNode
            if minimum(diff(zNode)) > hmin
                break
            end
        else
            break
        end
    end

    # get resistivity for the new layer by first finding the location nearest to
    # the new layer in the previous model, then perturb its resistivity
    layeridx = findLocation1D(zLoc, mcParam.zNode[1:currLayerNumber])
    mcParam.layeridx = layeridx
    mcParamProposed.layeridx = layeridx

    # assign resistivity of the new layer
    drho = mclimits.rhostd * randn()
    mcParamProposed.rho[propLayerNumber] =  mcParam.rho[layeridx] + drho

    # check if the new resistivity value in prior bounds
    mcParamProposed.inBound = checkPropertyBounds(mcParamProposed, mclimits)

    if iter == 1000
        mcParamProposed.inBound = false
    end

    return mcParamProposed

end


#-------------------------------------------------------------------------------
"""
    `mcDeath!(mcParam, mcParamProposed, mclimits)`

delete an exsiting layer from the current model.

"""
function mcDeath!(mcParam::MCParameter,
                  mcParamProposed::MCParameter,
                  mclimits::MCPrior)

    #
    currLayerNumber = mcParam.nLayer
    if currLayerNumber == mclimits.nlayermin
        mcParamProposed.inBound = false
        return
    end

    #
    propLayerNumber = currLayerNumber - 1
    mcParamProposed.nLayer = propLayerNumber

    # randomly choose an old layer to delete
    layeridx = unirandInteger(1, currLayerNumber)
    mcParam.layeridx = layeridx
    mcParamProposed.layeridx = layeridx

    # replace the deleted layer with the last layer
    zLoc = mcParam.zNode[currLayerNumber]
    mcParamProposed.zNode[layeridx] = zLoc
    mcParamProposed.rho[layeridx]   = mcParam.rho[currLayerNumber]

    # delete the chosen layer
    mcParamProposed.zNode[currLayerNumber:mclimits.nlayermax] .= 0.
    mcParamProposed.rho[currLayerNumber:mclimits.nlayermax]   .= 0.

    return mcParamProposed
end


#-------------------------------------------------------------------------------
"""
    `mcMoveLocation!(mcParam, mcParamProposed)`

move the location of one layer interface in the current model.

"""
function mcMoveLocation!(mcParam::MCParameter,
                         mcParamProposed::MCParameter,
                         mclimits::MCPrior)

    # randomly choose a layer to perturb
    layeridx = unirandInteger(1, mcParam.nLayer)
    mcParam.layeridx = layeridx
    mcParamProposed.layeridx = layeridx

    # perturb the location of the chosen layer interface
    mcParamProposed.zNode[layeridx] = mcParam.zNode[layeridx] + mclimits.zstd * randn()

    # check if the new location in prior bounds
    mcParamProposed.inBound = checkPropertyBounds(mcParamProposed, mclimits)

    return mcParamProposed
end


#-------------------------------------------------------------------------------
"""
    `mcPerturbProperty!(mcParam, mcParamProposed)`

perturb the resistivity of one layer in the current model

"""
function mcPerturbProperty!(mcParam::MCParameter,
                           mcParamProposed::MCParameter,
                           mclimits::MCPrior)

    # randomly choose a layer to perturb
    layeridx = unirandInteger(1, mcParam.nLayer)
    mcParam.layeridx = layeridx
    mcParamProposed.layeridx = layeridx

    # perturb the resistivity of the chosen layer
    mcParamProposed.rho[layeridx] = mcParam.rho[layeridx] + mclimits.mrhostd * randn()

    # check if the new location in prior bounds
    mcParamProposed.inBound = checkPropertyBounds(mcParamProposed, mclimits)

    return mcParamProposed

end


#-------------------------------------------------------------------------------
"""
    `checkPropertyBounds(mcParam, mclimits)`

check if the value of corresponding variables within specified bounds.

"""
function checkPropertyBounds(mcParam::MCParameter, mclimits::MCPrior)

    #
    inBound = true
    nLayer  = mcParam.nLayer

    # check location of layer interface
    if !all(x->(mclimits.zmin<=x<=mclimits.zmax), mcParam.zNode[1:nLayer])
        inBound = false
    end

    # check layer resistivity
    if !all(x->(mclimits.rhomin<=x<=mclimits.rhomax), mcParam.rho[1:nLayer])
        inBound = false
    end

    return inBound

end


#------------------------------------------------------------------------------
"""
    `selectChainStep(iterNo)`

randomly select a type of model perturbation for model proposal.

"""
function selectChainStep(iterNo::Int)

    #
    num = unirandInteger(1, 16)
    if num < 4
        step = 1
    elseif 3 < num < 7
        step = 2
    elseif 6 < num < 12
        step = 3
    else
        step = 4
    end

    return step

end #


function selectChainStep!(chainstep::MChainStep, iterNo::Int=1)

    #
    chainstep.isBirth = false
    chainstep.isDeath = false
    chainstep.isMove  = false
    chainstep.isPerturb = false

    #
    step = selectChainStep(iterNo)

    if step == 1
        chainstep.isBirth = true
    elseif step == 2
        chainstep.isDeath = true
    elseif step == 3
        chainstep.isMove  = true
    elseif step == 4
        chainstep.isPerturb = true
    else
        error("step $(step) is not supported")
    end

    return chainstep

end


#-------------------------------------------------------------------------------
include("acceptanceRatio.jl")
include("runMCMC.jl")
include("parallelTempering.jl")

end # module RJChains
