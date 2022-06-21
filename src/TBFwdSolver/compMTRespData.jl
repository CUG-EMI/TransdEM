export compMT1dResp, compMT1DImpedance
#-------------------------------------------------------------------------------
function compMT1dResp(emData::MTData, sig1D::Array, depth1D::Array)

    # frequencies
    freqs = emData.freqArray
    zimp  = compMT1DImpedance(freqs, sig1D, depth1D)

    nfreq = length(freqs)
    predData = zeros(2*nfreq)

    predData[1:2:end] = real.(zimp)
    predData[2:2:end] = imag.(zimp)

    dataID   = emData.dataID
    predData = predData[dataID]

    return predData

end


#-------------------------------------------------------------------------------
"""
    compMT1DImpedance(freqArray, sigma, depth1D)

calculates the impedance at the earth's surface for MT 1D layered model.
Note: length(sigma) = length(depth1D)-1

#Arguments
 - Input:
    * `freqArray` =::Array, frequencies
    * `sigma`     =::Array, conductivity model
    * `depth1D`     =::Array, location of top of each layer

 - Output:
    * `z1d`  =::Array, impedance at surface

"""
function compMT1DImpedance(freqArray::Array{T}, sigma::Array{T},
                           depth1D::Array{T}) where {T<:Float64}

    if length(sigma) != length(depth1D)-1
        error("layer's conductivity is not the same size with its depth.")
    end

    # physical constant
    mu0   = 4 * pi * 1e-7
    eps0  = 8.85 * 1e-12
    nFreq = length(freqArray)
    zLen  = diff(depth1D)
    nLayer = length(depth1D)

    # impedances at surface
    z1d = zeros(ComplexF64, nFreq)
    for i = 1:nFreq

        omega = 2 * pi * freqArray[i]
        k = sqrt(mu0 * eps0 * omega^2 - mu0 * sigma[end] * omega *1im)

        # calculate the impedance at the bottom layer
        ztmp = omega * mu0 / k

        for j = 2:nLayer
            ind = nLayer - j + 1
            k = sqrt(mu0 * eps0 * omega^2 - mu0 * sigma[ind] * omega * 1im)
            zp = omega * mu0 / k
            ztmp = zp * (ztmp + zp * tanh(k*zLen[ind]*1im)) / (zp + ztmp * tanh(k*zLen[ind]*1im))
        end
        z1d[i] = ztmp

    end

    return z1d

end
