export runMCMC
export MCMCsampling!
export resetChainModel!
export updateChainModel!
export updateChainArray!
export continueMCMCsampling
export readsamplefile

#------------------------------------------------------------------------------
function runMCMC(mclimits::MCPrior, emData::EMData)

    # initialize model parameter
    mcParamCurrent = initModelParameter(mclimits)
    mcDatafit = initMCDataFit()
    mcstep    = MChainStep(false, false, false, false)
    mcstatus  = MCStatus(false,zeros(Int,4),zeros(Int,4))

    # compute predicted data
    currPredData = getPredData(emData, mcParamCurrent)
    residuals = emData.obsData - currPredData

    # calculate likelihood and data misfit
    (currLikelihood,currMisfit) = compDataMisfit(emData, currPredData)

    mcDatafit.currPredData   = currPredData
    mcDatafit.currLikelihood = currLikelihood
    mcDatafit.currMisfit = currMisfit

    # initialize proposed model
    mcParamProposed = duplicateMCParameter(mcParamCurrent)

    # start MCMC sampling
    mcArray = MCMCsampling!(mcParamCurrent, mcParamProposed, mcstep,
                            mclimits, mcstatus, mcDatafit, emData)

    return mcArray, mcstatus

end


#------------------------------------------------------------------------------
"""
    `MCMCsampling(mcParamCurrent,mcParamProposed,mcstep,mclimits,mcstatus,
    mcDatafit,emData)`

produces an ensemble of samples with rjMCMC algorithm.

"""
function MCMCsampling!(mcParamCurrent::T,
                      mcParamProposed::T,
                      mcstep::MChainStep,
                      mclimits::MCPrior,
                      mcstatus::MCStatus,
                      mcDatafit::MCDataFit,
                      emData::EMData) where {T<:MCParameter}

    # mcmc information
    mcArray = initMChainArray(mclimits)

    # record starting model
    mcArray.startModel = duplicateMCParameter(mcParamCurrent)

    iterNo = 0
    cycleNumber  = 1000
    totalsamples = mclimits.totalsamples
    while (iterNo < totalsamples)

        iterNo = iterNo + 1
        # reset proposed models
        resetChainModel!(mcParamCurrent, mcParamProposed, mcDatafit)

        # randomly select a type of model perturbation
        mcstep = selectChainStep!(mcstep, iterNo)

        # which step do we have
        if mcstep.isBirth # birth step
            mcBirth!(mcParamCurrent, mcParamProposed, mclimits)

        elseif mcstep.isDeath # death step
            mcDeath!(mcParamCurrent, mcParamProposed, mclimits)

        elseif mcstep.isMove # move step
            mcMoveLocation!(mcParamCurrent, mcParamProposed, mclimits)

        elseif mcstep.isPerturb # perturb step
            mcPerturbProperty!(mcParamCurrent, mcParamProposed, mclimits)

        end

        # check if proposed model is within the prior bounds
        if mcParamProposed.inBound

            if mod(iterNo, cycleNumber) == 0
                println("iterNo=$(iterNo), dataMisfit=$(mcDatafit.currMisfit)")
            end

            # get the predicted data for the proposed model
            mcDatafit.propPredData = getPredData(emData, mcParamProposed)

            # get model likelihood and data misfit
            (mcDatafit.propLikelihood,mcDatafit.propMisfit) = compDataMisfit(emData,
            mcDatafit.propPredData)

            # calculate acceptance ratio of rjMcMC chain
            propAlpha = compAcceptanceRatio(mcParamCurrent, mcParamProposed,
            mcDatafit, mclimits, mcstep)

            # check if Metropolis-Hasting acceptance criterion is meeted
            checkMHCriterion!(propAlpha, mcstatus, mcstep)

        end #

        # update markov chain
        if mcParamProposed.inBound & mcstatus.accepted
            updateChainModel!(mcParamCurrent, mcParamProposed, mcDatafit)

            # data residuals
            mcDatafit.dataResiduals = emData.obsData - mcDatafit.propPredData
        end

        # record in chainarray
        updateChainArray!(mcArray, mcParamCurrent, mcstep, mcstatus, mcDatafit)

    end

    return mcArray

end


#------------------------------------------------------------------------------
"""
    `resetChainModel(mcParam, mcParamProposed, mcDatafit)`

reset proposed chain model.

"""
function resetChainModel!(mcParam::T,
                          mcParamProposed::T,
                          mcDatafit::MCDataFit) where{T<:MCParameter}

    # reset proposed model
    updateMCParameter!(mcParam, mcParamProposed)

    # update predicted data for proposed model
    mcDatafit.propPredData    = copy(mcDatafit.currPredData)
    mcDatafit.propLikelihood  = mcDatafit.currLikelihood
    mcDatafit.propMisfit  = mcDatafit.currMisfit

    return mcParamProposed, mcDatafit

end


#-------------------------------------------------------------------------------
"""
    `updateChainModel!(mcParamCurrent, rjdatafit)`

"""
function updateChainModel!(mcParamCurrent::MCParameter,
                          mcParamProposed::MCParameter,
                          rjdatafit::MCDataFit)

    # model update
    mcParamCurrent = updateMCParameter!(mcParamProposed, mcParamCurrent)

    #
    rjdatafit.currPredData   = copy(rjdatafit.propPredData)
    rjdatafit.currLikelihood = copy(rjdatafit.propLikelihood)
    rjdatafit.currMisfit     = copy(rjdatafit.propMisfit)

    return mcParamCurrent, rjdatafit

end


#-------------------------------------------------------------------------------
"""
    `updateChainArray(mcArray, iterNo, mcParamCurrent, mcstep, mcstatus, mcDatafit)`

"""
function updateChainArray!(mcArray::MChainArray,
                          mcParamCurrent::MCParameter,
                          mcstep::MChainStep,
                          mcstatus::MCStatus,
                          mcDatafit::MCDataFit)

    #
    mcArray.nsample += 1
    iterNo   = mcArray.nsample
    nLayer   = mcParamCurrent.nLayer
    layeridx = mcParamCurrent.layeridx
    mcArray.chainnlayer[iterNo]  = nLayer
    mcArray.chainMisfit[iterNo]  = mcDatafit.currMisfit

    if mcstep.isBirth
        mcArray.chainstep[iterNo, 1] = 1
        if mcstatus.accepted
            mcArray.chainstep[iterNo, 2] = 1
            mcArray.chainindice[iterNo]  = 0
            mcArray.chainvalue[iterNo, 1] = mcParamCurrent.zNode[nLayer]
            mcArray.chainvalue[iterNo, 2] = mcParamCurrent.rho[nLayer]
        end

    elseif mcstep.isDeath
        mcArray.chainstep[iterNo, 1] = 2
        if mcstatus.accepted
            mcArray.chainstep[iterNo, 2] = 1
            mcArray.chainindice[iterNo]  = layeridx
            mcArray.chainvalue[iterNo, 1] = 0.
            mcArray.chainvalue[iterNo, 2] = 0.
        end

    elseif mcstep.isMove
        mcArray.chainstep[iterNo, 1] = 3
        if mcstatus.accepted
            mcArray.chainstep[iterNo, 2] = 1
            mcArray.chainindice[iterNo] = layeridx
            mcArray.chainvalue[iterNo, 1] = mcParamCurrent.zNode[layeridx]
            mcArray.chainvalue[iterNo, 2] = 0.
        end

    elseif mcstep.isPerturb
        mcArray.chainstep[iterNo, 1] = 4
        if mcstatus.accepted
            mcArray.chainstep[iterNo, 2] = 1
            mcArray.chainindice[iterNo] = layeridx
            mcArray.chainvalue[iterNo, 1] = 0.
            mcArray.chainvalue[iterNo, 2] = mcParamCurrent.rho[layeridx]
        end

    end

    # data residuals
    # if mcstatus.accepted
    #     mcArray.residuals[:, iterNo] = mcDatafit.dataResiduals
    # end

    # reset status
    mcstatus.accepted = false

end


#------------------------------------------------------------------------------
"""
    `continueMCMCsampling(samplefile, mclimits, emData)`

continue a MCMC sampling from previous Markov chain.

"""
function continueMCMCsampling(samplefile::String, mclimits::MCPrior, emData::EMData)

    # initialize model parameter
    mcDatafit = initMCDataFit()
    mcstep    = MChainStep(false, false, false, false)
    mcstatus  = MCStatus(false,zeros(Int,4),zeros(Int,4))

    # get the last sample of previous Markov chain
    (mcParamCurrent, mcArrayPrev) = readsamplefile(samplefile)

    # compute predicted data
    currPredData = getPredData(emData, mcParamCurrent)
    residuals = emData.obsData - currPredData

    # calculate likelihood and data misfit
    (currLikelihood,currMisfit) = compDataMisfit(emData, currPredData)

    mcDatafit.currPredData   = currPredData
    mcDatafit.currLikelihood = currLikelihood
    mcDatafit.currMisfit = currMisfit

    # initialize proposed model
    mcParamProposed = duplicateMCParameter(mcParamCurrent)

    # start MCMC sampling
    mcArray = MCMCsampling!(mcParamCurrent, mcParamProposed, mcstep,
                            mclimits, mcstatus, mcDatafit, emData)

    return mcArray, mcstatus


end


#------------------------------------------------------------------------------
"""
    `readsamplefile(samplefile, mclimits)`

get the last sample of the Markov chain, which will be used as the current sample
for a new Markov chain.

"""
function readsamplefile(samplefile::String, mclimits::MCPrior)

    #
    if isfile(samplefile)
        fid = open(samplefile, "r")
    else
        error("$(samplefile) does not exist, please try again.")
    end

    # mcmc information
    mcArray = initMChainArray(mclimits)
    zNode   = []
    rho     = []
    nsample = 0
    while !eof(fid)
        cline = strip(readline(fid))

        # ignore all comments: empty line, or line preceded with !
        while cline[1] == '!' || isempty(cline)
            cline = strip(readline(fid))
        end
        cline = lowercase(cline)

        if occursin("zcoordinate", cline)
            tmp = split(cline)
            nz  = parse(Int, tmp[end])
            zNode = zeros(Float64, nz)
            cline = strip(readline(fid))
            cline = split(cline)
            for k = 1:nz
                zNode[k] = parse(Float64, cline[k])
            end

        elseif occursin("resistivity", cline)
            tmp = split(cline)
            nz  = parse(Int, tmp[end])
            rho = zeros(Float64, nz)
            cline = strip(readline(fid))
            cline = split(cline)
            for k = 1:nz
                rho[k] = parse(Float64, cline[k])
            end

        elseif occursin("samplelist", cline)
            tmp = split(cline)
            nsample = parse(Int, tmp[end])
            for k = 1:nsample
                cline = strip(readline(fid))
                cline = split(cline)
                mcArray.chainstep[k, :] = parse.(Float64, cline[2:3])
                mcArray.chainvalue[k,:] = parse.(Float64, cline[4:5])
                mcArray.chainindice[k]  = parse(Float64, cline[6])
                mcArray.chainncell[k]   = parse(Float64, cline[7])
                mcArray.chainMisfit[k]  = parse(Float64, cline[8])
            end

        end #

    end # while

    close(fid)

    # update model parameters
    for k = 1:nsample
        updatesampleModel!(zNode, rho, k, mcArray)
    end

    mcParam = initModelParameter(mclimits)
    mcParam.zNode = copy(zNode)
    mcParam.rho   = copy(rho)

    return mcParam, mcArray

end
