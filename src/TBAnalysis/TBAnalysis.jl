#-------------------------------------------------------------------------------
#
# module `TBAnalysis` defines routines to perform posterior analysis of model
# parameters.
#
#
#-------------------------------------------------------------------------------
module TBAnalysis

using TransdEM.TBUtility
using TransdEM.TBStruct
using TransdEM.TBChain

using LinearAlgebra, Printf, Statistics
using Distributed, DistributedArrays

export sampleStatistics
export multipleStatistics
export updatesampleModel!
export outputMCMCsamples
export outputEnsembleModel
export getCredibleIntervalModel
export mksampleArray
export extractModelParam

export outputPTDataMisfit
export outputPTchainHeader
export outputPTchains
export samplePTStatistics
export clearPTrefence

#-------------------------------------------------------------------------------
"""
    `outputMCMCsamples(mcArray)`

"""
function outputMCMCsamples(mcArray::MChainArray, ichain::Int=1)

    nsample = length(mcArray.chainMisfit)
    chainnlayer = mcArray.chainnlayer

    # chain values
    valuefile = "chainsamples_id$(ichain).dat"
    valID = open(valuefile, "w")
    chainstep = mcArray.chainstep
    chainindice = mcArray.chainindice
    chainvalue  = mcArray.chainvalue

    # sart model
    @printf(valID,"! starting model \n")
    startModel  = mcArray.startModel
    nlayer   = startModel.nLayer
    zNode   = startModel.zNode
    rho = startModel.rho
    @printf(valID, "zcoordinate: %5d \n", nlayer)
    for k = 1:nlayer
        @printf(valID, "%8g ", zNode[k])
    end
    @printf(valID, "\n")
    @printf(valID, "resistivity: %5d \n", nlayer)
    for k = 1:nlayer
        @printf(valID, "%8g ", rho[k])
    end
    @printf(valID, "\n")

    @printf(valID,"! #sample  step  accept layeridx zcoord(lg10)  rho(lg10)  nLayer  dataMisfit\n")
    @printf(valID, "samplelist: %8d\n", nsample)
    for i = 1:nsample
        @printf(valID, "%8d %6d %6d %6d ", i, chainstep[i,1], chainstep[i,2], chainindice[i])
        @printf(valID, "%8g %8g ", chainvalue[i, 1], chainvalue[i, 2])
        @printf(valID, "%8d %12g \n", chainnlayer[i], mcArray.chainMisfit[i])
    end
    close(valID)

end


#-------------------------------------------------------------------------------
"""
    `sampleStatistics(mcArray, mclimits)`

"""
function sampleStatistics(mcArray::MChainArray, mclimits::MCPrior,
                          filestr::String="mc")

    #
    burnin    = mclimits.burninsamples
    nsamples  = mclimits.totalsamples
    nlayermin = mclimits.nlayermin
    nlayermax = mclimits.nlayermax
    zmin     = mclimits.zmin
    zmax     = mclimits.zmax
    rhomin   = mclimits.rhomin
    rhomax   = mclimits.rhomax

    # start model
    startModel  = mcArray.startModel
    zNode   = copy(startModel.zNode)
    rho     = copy(startModel.rho)

    # intervals
    nzBins = mclimits.nzBins
    npBins = mclimits.npBins
    (zLocBins, zspace) = mksampleArray(zmin, zmax, nzBins)
    (rhoBins,rhospace) = mksampleArray(rhomin, rhomax, npBins)

    for k = 1:burnin
        updatesampleModel!(zNode, rho, k, mcArray)
    end

    # calculate posteriors
    rhoModel = zeros(nzBins)
    nlayerPosterior = zeros(Int, nlayermax)
    depthPosterior  = zeros(Int, nzBins)
    valuePosterior  = zeros(Int, nzBins, npBins)

    # statistical moment
    sm01 = zeros(nzBins)
    sm02 = zeros(nzBins)

    # collecte mcmc samples after the burnin period
    samplestep = 100
    # number of samples in ensemble
    nEnsemble  = 0
    for i in burnin+1:nsamples
        updatesampleModel!(zNode, rho, i, mcArray)

        #
        if mod(i, samplestep) > 0
            continue
        end

        nEnsemble += 1
        # posterior layer number
        nlayer  = mcArray.chainnlayer[i]
        nlayerPosterior[nlayer] += 1

        # get depth posterior
        for j = 1:nlayer
            zLoc = zNode[j]
            cidx = findLocation1D(zLoc, zLocBins)
            depthPosterior[cidx] += 1
        end

        # get values at sampling location
        for j = 1:nzBins
            zLoc = zLocBins[j]
            cidx = findNearest1D(zLoc, zNode[1:nlayer])
            rhoModel[j] = rho[cidx]
            #
            idx = findLocation1D(rho[cidx], rhoBins)
            valuePosterior[j, idx] += 1
        end

        # statistical moments
        sm01 += rhoModel
        sm02 += rhoModel .^ 2

    end

    # evaluate statistical variables
    sm01 = sm01 / nEnsemble
    sm02 = sm02 / nEnsemble

    # arithmetic mean
    valueMean = sm01

    # standard deviation
    std = sm02 - sm01 .^ 2
    std[std .< 0] .= eps(1.0)
    valueStd  = sqrt.(std)

    # mode = maximum a posterior
    valueMode = zeros(nzBins)
    for i = 1:nzBins
        (val,idx) = findmax(valuePosterior[i, :])
        valueMode[i] = rhoBins[idx]
    end

    # median
    valueMedian = zeros(nzBins)
    for i = 1:nzBins
        count = 0
        for j = 1:npBins
            count += valuePosterior[i, j]
            if count > nEnsemble / 2
                valueMedian[i] = rhoBins[j]
                break
            end
        end
    end

    # output results
    outputPosterior(zLocBins, rhoBins, nlayerPosterior, depthPosterior,
    valuePosterior, filestr)

    # credible min and max model
    ci = mclimits.credInterval
    (valueMin,valueMax) = getCredibleIntervalModel(valuePosterior, rhoBins,
                                                   nEnsemble, ci)
    #
    filename = "posteriorModel-"*filestr*".dat"
    results  = hcat(zLocBins,valueMean,valueMedian,valueMode,valueMin,valueMax,valueStd)
    outputEnsembleModel(filename, results)

end # sampleStatistics


#-------------------------------------------------------------------------------
"""
    `updatesampleModel!(zNode, rho, k, mcArray)`

"""
function updatesampleModel!(zNode::Vector{T}, rho::Vector{T}, k::Int,
                           mcArray::MChainArray) where{T<:Float64}

    mcstep = mcArray.chainstep[k, 1]
    accept = mcArray.chainstep[k, 2]
    nlayer = mcArray.chainnlayer[k]
    cidx   = mcArray.chainindice[k]
    nlayermax = length(zNode)

    # birth step
    if mcstep==1 && accept>0
        zNode[nlayer] = mcArray.chainvalue[k, 1]
        rho[nlayer]   = mcArray.chainvalue[k, 2]
        zNode[nlayer+1:nlayermax] .= 0.
        rho[nlayer+1:nlayermax]   .= 0.

    # death step
    elseif mcstep==2 && accept>0
        zNode[cidx] = zNode[nlayer+1]
        rho[cidx]   = rho[nlayer+1]
        zNode[nlayer+1:nlayermax] .= 0.
        rho[nlayer+1:nlayermax]   .= 0.

    # move step
    elseif mcstep==3 && accept>0
        zNode[cidx] = mcArray.chainvalue[k, 1]

    # perturb step
    elseif mcstep==4 && accept>0
        rho[cidx] = mcArray.chainvalue[k, 2]

    end

end # updatesampleModel


#-------------------------------------------------------------------------------
"""
    `outputEnsembleModel(filename, coords, results)`

"""
function outputEnsembleModel(filename::String, coords::Vector, results::Vector)

    fileID = open(filename, "w")

    for k = 1:length(coords)
        @printf(fileID, "%12g %12g \n", coords[k], results[k])
    end

    close(fileID)

end


#-------------------------------------------------------------------------------
"""
    `outputEnsembleModel(filename, results, comments)`

"""
function outputEnsembleModel(filename::String, results::Array{Float64}, comments::String="")

    fileID = open(filename, "w")
    if !isempty(comments)
        @printf(fileID, "%s\n", comments)
    end

    (nBins, nm) = size(results)
    for i = 1:nBins
        for j = 1:nm
            @printf(fileID, "%10g ", results[i,j])
        end
        @printf(fileID, "\n")
    end

    close(fileID)

end


#-------------------------------------------------------------------------------
"""
    `outputPosterior(coordBins, rhoBins, nlayerPosterior, depthPosterior,
    valuePosterior)`

"""
function outputPosterior(coordBins, rhoBins, nlayerPosterior,
    depthPosterior, valuePosterior, filestr)

    filename ="depthBins-"*filestr*".dat"
    fid = open(filename, "w")
    for k = 1:length(coordBins)
        @printf(fid, "%12g\n", coordBins[k])
    end
    close(fid)

    filename ="rhoBins-"*filestr*".dat"
    fid = open(filename, "w")
    for k = 1:length(rhoBins)
        @printf(fid, "%12g\n", rhoBins[k])
    end
    close(fid)

    filename ="nlayerHistogram-"*filestr*".dat"
    fid = open(filename, "w")
    for k = 1:length(nlayerPosterior)
        @printf(fid, "%4d %8d\n", k, nlayerPosterior[k])
    end
    close(fid)

    filename = "depthHistogram-"*filestr*".dat"
    fid = open(filename, "w")
    for k = 1:length(depthPosterior)
        @printf(fid, "%8d \n", depthPosterior[k])
    end
    close(fid)

    filename = "depthrhoHistogram-"*filestr*".dat"
    fid = open(filename, "w")
    for i = 1:length(rhoBins)
        for j = 1:length(coordBins)
            @printf(fid, "%8d ", valuePosterior[j, i])
        end
        @printf(fid, "\n")
    end
    close(fid)

end


#-------------------------------------------------------------------------------
"""
    `getCredibleIntervalModel(valuePosterior, rhoBins, nEnsemble, credInterval)`

get the minimum and maximum credible model.

"""
function getCredibleIntervalModel(valuePosterior::Array{Ti},
                                  rhoBins::Vector{Tv},
                                  nEnsemble::Ti,
                                  credInterval::Tv) where {Ti<:Int,Tv<:Float64}

    #
    (nzBins,npBins) = size(valuePosterior)
    valueMin = zeros(nzBins)
    valueMax = zeros(nzBins)
    credsample = nEnsemble * (1 - credInterval) / 2

    # minimum credible model
    for i = 1:nzBins
        count = 0
        for j = 1:npBins
            count += valuePosterior[i, j]
            if count > credsample
                valueMin[i] = rhoBins[j]
                break
            end
        end
    end

    # maximum credible model
    for i = 1:nzBins
        count = 0
        for j = npBins:-1:1
            count += valuePosterior[i, j]
            if count > credsample
                valueMax[i] = rhoBins[j]
                break
            end
        end
    end

    return valueMin, valueMax

end


#-------------------------------------------------------------------------------
"""
    `mksampleArray(vmin, vmax, nsample)`

"""
function mksampleArray(vmin::T, vmax::T, nsample::Int) where{T<:Number}

    dx = (vmax-vmin) / nsample
    x1 = vmin + 0.5 * dx
    x2 = vmax - 0.5 * dx
    val = LinRange(x1, x2, nsample)

    return collect(val), dx

end


#-------------------------------------------------------------------------------
"""
    `multipleStatistics(mcArray, mclimits, filestr)`

"""
function multipleStatistics(results::Array, mclimits::MCPrior,
                            filestr::String="mc", rscale::Float64=1.3)

    #
    burnin    = mclimits.burninsamples
    nsamples  = mclimits.totalsamples
    nlayermin = mclimits.nlayermin
    nlayermax = mclimits.nlayermax
    zmin     = mclimits.zmin
    zmax     = mclimits.zmax
    rhomin   = mclimits.rhomin
    rhomax   = mclimits.rhomax

    # multiple chains
    nchain = length(results)

    # calculate posteriors
    nzBins = mclimits.nzBins
    npBins = mclimits.npBins
    (zLocBins, zspace) = mksampleArray(zmin, zmax, nzBins)
    (rhoBins,rhospace) = mksampleArray(rhomin, rhomax, npBins)

    # calculate posteriors
    rhoModel = zeros(nzBins)
    nlayerPosterior = zeros(Int, nlayermax)
    depthPosterior  = zeros(Int, nzBins)
    valuePosterior  = zeros(Int, nzBins, npBins)

    # statistical moment
    sm01 = zeros(nzBins)
    sm02 = zeros(nzBins)

    # collecte mcmc samples after the burnin period
    samplestep = 100

    # number of samples in ensemble
    nEnsemble  = 0

    # compute average data misfit
    aveMisfit = zeros(nchain)
    for ic = 1:nchain
        aveMisfit[ic] = mean(results[ic].chainMisfit[burnin:nsamples])
    end
    aveMisfit = median(aveMisfit)

    for ic = 1:nchain

        mcArray = results[ic]
        # exclude chains that are not convergent
        mdatafit = mean(mcArray.chainMisfit[burnin:nsamples])
        if mdatafit > rscale * aveMisfit
            println("chain No.$(ic) has been excluded!")
            continue
        end
        #
        println("Output sampling results for chain No.$(ic) ...")
        outputMCMCsamples(mcArray, ic)

        # start model
        startModel  = mcArray.startModel
        zNode   = copy(startModel.zNode)
        rho     = copy(startModel.rho)

        # update model parameters
        for k = 1:burnin
            updatesampleModel!(zNode, rho, k, mcArray)
        end

        for i in burnin+1:nsamples
            updatesampleModel!(zNode, rho, i, mcArray)
            #
            if mod(i, samplestep) > 0; continue; end
            nEnsemble += 1

            # posterior cell number
            nlayer  = mcArray.chainnlayer[i]
            nlayerPosterior[nlayer] += 1

            # get depth posterior
            for j = 1:nlayer
                zLoc = zNode[j]
                cidx = findLocation1D(zLoc, zLocBins)
                depthPosterior[cidx] += 1
            end

            # get values at sampling location
            for j = 1:nzBins
                zLoc = zLocBins[j]
                cidx = findNearest1D(zLoc, zNode[1:nlayer])
                rhoModel[j] = rho[cidx]
                idx = findLocation1D(rho[cidx], rhoBins)
                valuePosterior[j, idx] += 1

            end

            # statistical moments
            sm01 += rhoModel
            sm02 += rhoModel .^ 2

        end

    end # nchain

    # evaluate statistical variables
    sm01 = sm01 / nEnsemble
    sm02 = sm02 / nEnsemble

    # arithmetic mean
    valueMean = sm01

    # standard deviation
    std = sm02 - sm01 .^ 2
    std[std .< 0] .= eps(1.0)
    valueStd  = sqrt.(std)

    # mode = maximum a posterior
    valueMode = zeros(nzBins)
    for i = 1:nzBins
        (val,idx) = findmax(valuePosterior[i,:])
        valueMode[i] = rhoBins[idx]
    end

    # meadian
    valueMedian = zeros(nzBins)
    for i = 1:nzBins
        count = 0
        for j = 1:npBins
            count += valuePosterior[i, j]
            if count > nEnsemble / 2
                valueMedian[i] = rhoBins[j]
                break
            end
        end
    end

    # output results
    outputPosterior(zLocBins, rhoBins, nlayerPosterior, depthPosterior,
    valuePosterior, filestr)

    # credible min and max model
    ci = mclimits.credInterval
    (valueMin,valueMax) = getCredibleIntervalModel(valuePosterior, rhoBins,
                                                   nEnsemble, ci)

    #
    filename = "posteriorModel-"*filestr*".dat"
    results  = hcat(zLocBins,valueMean,valueMedian,valueMode,valueMin,valueMax,valueStd)
    outputEnsembleModel(filename, results)

end # multipleStatistics


#-------------------------------------------------------------------------------
"""
    `extractModelParam(mcArray, mclimits, nLayer)`

"""
function extractModelParam(mcArray::MChainArray, mclimits::MCPrior, nLayer::Int)

    #
    burnin    = mclimits.burninsamples
    nsamples  = mclimits.totalsamples

    # start model
    startModel  = mcArray.startModel
    zNode   = copy(startModel.zNode)
    rho     = copy(startModel.rho)

    for k = 1:burnin
        updatesampleModel!(zNode, rho, k, mcArray)
    end

    # collecte mcmc samples after the burnin period
    samplestep = 100

    #
    nidx = findall(mcArray.chainnlayer[burnin:samplestep:nsamples] .== nLayer)
    num  = length(nidx)
    nrow = 2*nLayer - 1
    km   = 1
    modparam = zeros(nrow, num)
    for i in burnin+1:nsamples

        updatesampleModel!(zNode, rho, i, mcArray)
        if mod(i, samplestep) > 0; continue; end

        # posterior layer number
        pnlayer = mcArray.chainnlayer[i]
        if pnlayer == nLayer
            idx = sortperm(zNode[1:nLayer])
            dep1D = zNode[idx]
            modparam[1:nLayer, km]     = copy(rho[idx])
            modparam[nLayer+1:nrow,km] = copy(dep1D[1:nLayer-1])
            km = km + 1
        end

    end

    return modparam

end


#-------------------------------------------------------------------------------
"""
    `extractModelParam(mcArray, mclimits, nLayer)`

"""
function extractModelParam(results::Array, mclimits::MCPrior, nLayer::Int)

    # number of chains
    nchain = length(results)

    filename = "modelParam_$(nLayer)layer.dat"
    fileID   = open(filename, "w")
    for ic = 1:nchain
        modparam = extractModelParam(results[ic], mclimits, nLayer)
        (nr,nc)  = size(modparam)
        for j = 1:nc
            for k = 1:nr
                @printf(fileID, "%8g ", modparam[k,j])
            end
            @printf(fileID, "\n")
        end
    end
    close(fileID)

end


#-------------------------------------------------------------------------------
"""
    `outputPTDataMisfit(mcArrayRef, tempData, mclimits, filestr)`

output data misfit of the samples in each Markov chain for convergence analysis.

"""
function outputPTDataMisfit(mcArrayRef::DArray, tempData::TBTemperature,
                            mclimits::MCPrior, filestr::String="mc")

    #
    nsamples  = mclimits.totalsamples
    nchain = length(mcArrayRef)

    #
    dataMisfit = zeros(nsamples, nchain)
    swapchain  = tempData.swapchain

    for ic = 1:nchain
        dataMisfit[:, ic] = mcArrayRef[ic].chainMisfit
    end

    #
    for k = 1:nsamples
        for ic = 1:nchain
            if swapchain[3,ic,k] > 0
                idx01 = [swapchain[1,ic,k],swapchain[2,ic,k]]
                idx02 = [swapchain[2,ic,k],swapchain[1,ic,k]]
                chainidx[idx01] = chainidx[idx02]
            end
        end
        ind = indexin(sortchain, chainidx)
        dataMisfit[k, :] = dataMisfit[k, ind]
    end

    # output
    filename = "chainMisfit-"*filestr*".dat"
    fileID   = open(filename, "w")
    for k = 1:nsamples
        for ic = 1:nchain
            @printf(fileID, "%8g ", dataMisfit[k,ic])
        end
        @printf(fileID, "\n")
    end

    close(fileID)

end


#-------------------------------------------------------------------------------
"""
    `outputPTchainHeader(mcArray, chainidx)`

"""
function outputPTchainHeader(mcArray::Array{MChainArray}, chainidx::Vector, filestr::String="mc")

    #
    nchain = length(chainidx)
    valID  = Array{IOStream}(undef, nchain)
    for ic = 1:nchain

        id = chainidx[ic]
        valuefile = filestr*"_chainsamples_id$(ic).dat"
        valID[ic] = open(valuefile, "w")
        chainstep = mcArray[id].chainstep
        chainindice = mcArray[id].chainindice
        chainvalue  = mcArray[id].chainvalue
        nsamples    = length(mcArray[id].chainMisfit)
        # sart model
        @printf(valID[ic],"! starting model \n")
        startModel  = mcArray[id].startModel
        nlayer   = startModel.nLayer
        zNode   = startModel.zNode
        rho = startModel.rho
        @printf(valID[ic], "zcoordinate: %5d \n", nlayer)
        for k = 1:nlayer
            @printf(valID[ic], "%8g ", zNode[k])
        end
        @printf(valID[ic], "\n")
        @printf(valID[ic], "resistivity: %5d \n", nlayer)
        for k = 1:nlayer
            @printf(valID[ic], "%8g ", rho[k])
        end
        @printf(valID[ic], "\n")

        #
        @printf(valID[ic],"! #sample  step  accept layeridx zcoord(lg10)  rho(lg10)  nLayer  dataMisfit\n")
        @printf(valID[ic], "samplelist: %8d\n", nsamples)
        k = 1
        @printf(valID[ic], "%8d %6d %6d %6d ", k, mcArray[id].chainstep[k, 1],
        mcArray[ic].chainstep[k, 2], mcArray[id].chainindice[k])
        @printf(valID[ic], "%8g %8g ", mcArray[id].chainvalue[k,1],
        mcArray[id].chainvalue[k, 2])
        @printf(valID[ic], "%8d %12g \n", mcArray[id].chainnlayer[k],
        mcArray[id].chainMisfit[k])

    end

    return valID

end


#-------------------------------------------------------------------------------
"""
    `outputPTchains(mcArrayRef, tempData)`

write out Markov chains in parallel tempering into file.

"""
function outputPTchains(mcArrayRef::DArray, tempData::TBTemperature, filestr::String="mc")

    #
    nchain     = length(mcArrayRef)
    swapchain  = tempData.swapchain
    chainidx   = collect(1:nchain)
    sortchain  = collect(1:nchain)
    mcArray    = Array{MChainArray}(undef, nchain)
    valID      = Array{IOStream}(undef, nchain)

    #
    for ic = 1:nchain
        mcArray[ic] = mcArrayRef[ic]
    end
    nsamples = length(mcArray[1].chainMisfit)
    for k = 1:nsamples

        for ic = 1:nchain
            if swapchain[3,ic,k] > 0
                idx01 = [ swapchain[1,ic,k], swapchain[2,ic,k] ]
                idx02 = [ swapchain[2,ic,k], swapchain[1,ic,k] ]
                chainidx[idx01] = chainidx[idx02]
            end
        end
        ind = indexin(sortchain, chainidx)

        # output samples
        if k < 2
            valID = outputPTchainHeader(mcArray, ind, filestr)
        else
            for ic = 1:nchain
                id = ind[ic]
                @printf(valID[ic], "%8d %6d %6d %6d ", k, mcArray[id].chainstep[k, 1],
                mcArray[ic].chainstep[k, 2], mcArray[id].chainindice[k])
                @printf(valID[ic], "%8g %8g ", mcArray[id].chainvalue[k,1],
                mcArray[id].chainvalue[k, 2])
                @printf(valID[ic], "%8d %12g \n", mcArray[id].chainnlayer[k],
                mcArray[id].chainMisfit[k])
            end
        end
    end

    #
    for ic = 1:nchain
        close(valID[ic])
    end

end


#-------------------------------------------------------------------------------
"""
    `samplePTStatistics(mcArrayRef, tempData, mclimits, filestr)`

calculate statistical models from parallel tempering MCMC results.

"""
function samplePTStatistics(mcArrayRef::DArray, tempData::TBTemperature,
                            mclimits::MCPrior, filestr::String="mc")

    #
    burnin    = mclimits.burninsamples
    nsamples  = mclimits.totalsamples
    nlayermin = mclimits.nlayermin
    nlayermax = mclimits.nlayermax
    zmin     = mclimits.zmin
    zmax     = mclimits.zmax
    rhomin   = mclimits.rhomin
    rhomax   = mclimits.rhomax

    # intervals
    nzBins = mclimits.nzBins
    npBins = mclimits.npBins
    (zLocBins, zspace) = mksampleArray(zmin, zmax, nzBins)
    (rhoBins,rhospace) = mksampleArray(rhomin, rhomax, npBins)

    # multiple chains
    nchain = length(mcArrayRef)

    # calculate posteriors
    rhoModel = zeros(nzBins)
    nlayerPosterior = zeros(Int, nlayermax)
    depthPosterior  = zeros(Int, nzBins)
    valuePosterior  = zeros(Int, nzBins, npBins)

    # statistical moment
    sm01 = zeros(nzBins)
    sm02 = zeros(nzBins)

    # collecte mcmc samples after the burnin period
    samplestep = 100

    # number of samples in ensemble
    nEnsemble  = 0
    tempLadder = copy(tempData.tempLadder)
    chainidx   = collect(1:nchain)
    sortchain  = collect(1:nchain)

    #
    mcArray    = Array{MChainArray}(undef, nchain)
    dataMisfit = zeros(nsamples, nchain)
    swapchain  = tempData.swapchain
    zNodeArray = Array{Any}(undef, nchain)
    rhoArray   = Array{Any}(undef, nchain)
    for ic = 1:nchain
        mcArray[ic] = mcArrayRef[ic]
        dataMisfit[:, ic] = mcArray[ic].chainMisfit
        zNodeArray[ic] = copy(mcArray[ic].startModel.zNode)
        rhoArray[ic]   = copy(mcArray[ic].startModel.rho)
    end

    #
    for k = 1:nsamples

        for ic = 1:nchain
            # update model parameters
            updatesampleModel!(zNodeArray[ic], rhoArray[ic], k, mcArray[ic])
            if swapchain[3,ic,k] > 0
                idx01 = [ swapchain[1,ic,k], swapchain[2,ic,k] ]
                idx02 = [ swapchain[2,ic,k], swapchain[1,ic,k] ]
                chainidx[idx01] = chainidx[idx02]

                # update chain temperature
                tempLadder[idx01] = tempLadder[idx02]
            end
        end
        ind = indexin(sortchain, chainidx)
        dataMisfit[k, :] = dataMisfit[k, ind]

        # collect samples after burnin period
        if k <= burnin
            continue
        else
            if mod(k, samplestep) > 0; continue; end

            for ic = 1:nchain
                if tempLadder[ic] > 1.0; continue; end
                nEnsemble += 1

                # posterior cell number
                nlayer  = mcArray[ic].chainnlayer[k]
                nlayerPosterior[nlayer] += 1

                # get depth posterior
                for j = 1:nlayer
                    zLoc = zNodeArray[ic][j]
                    cidx = findLocation1D(zLoc, zLocBins)
                    depthPosterior[cidx] += 1
                end

                # get values at sampling location
                for j = 1:nzBins
                    zLoc = zLocBins[j]
                    cidx = findNearest1D(zLoc, zNodeArray[ic][1:nlayer])
                    rhoModel[j] = rhoArray[ic][cidx]
                    idx = findLocation1D(rhoArray[ic][cidx], rhoBins)
                    valuePosterior[j, idx] += 1

                end

                # statistical moments
                sm01 += rhoModel
                sm02 += rhoModel .^ 2

            end # ic

        end # if

    end # k

    # evaluate statistical variables
    sm01 = sm01 / nEnsemble
    sm02 = sm02 / nEnsemble

    # arithmetic mean
    valueMean = sm01

    # standard deviation
    std = sm02 - sm01 .^ 2
    std[std .< 0] .= eps(1.0)
    valueStd  = sqrt.(std)

    # mode = maximum a posterior
    valueMode = zeros(nzBins)
    for i = 1:nzBins
        (val,idx) = findmax(valuePosterior[i,:])
        valueMode[i] = rhoBins[idx]
    end

    # meadian
    valueMedian = zeros(nzBins)
    for i = 1:nzBins
        count = 0
        for j = 1:npBins
            count += valuePosterior[i, j]
            if count > nEnsemble / 2
                valueMedian[i] = rhoBins[j]
                break
            end
        end
    end

    # output results
    outputPosterior(zLocBins, rhoBins, nlayerPosterior, depthPosterior,
    valuePosterior, filestr)

    # credible min and max model
    ci = mclimits.credInterval
    (valueMin,valueMax) = getCredibleIntervalModel(valuePosterior, rhoBins,
                                                   nEnsemble, ci)

    # output results
    filename = "posteriorModel-"*filestr*".dat"
    results  = hcat(zLocBins,valueMean,valueMedian,valueMode,valueMin,valueMax,valueStd)
    comments = "#z(lg10) mean(lg10) median(lg10) mode(lg10) credmin(lg10) credmax(lg10) std"
    outputEnsembleModel(filename, results, comments)

    # data misfit
    filename = "chainMisfit-"*filestr*".dat"
    outputEnsembleModel(filename, dataMisfit)

    return nEnsemble

end

#------------------------------------------------------------------------------
"""
    `clearPTrefence()`

release all darrays created during the inversion explicitly.

"""
function clearPTrefence(mcRef::DArray, stRef::DArray, dfRef::DArray, paramRef::DArray)

    #
    close(mcRef)
    close(stRef)
    close(dfRef)
    close(paramRef)
    mcRef = []
    stRef = []
    dfRef = []
    paramRef = []
    d_closeall()

end

#-------------------------------------------------------------------------------

end # TBAnalysis
