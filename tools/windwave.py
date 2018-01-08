import numpy as np
import math

def main():

    ustokes, vstokes = stokes_drift(freq, spec, xcmp, ycmp, z)

def stokes_drift(freq, spec, xcmp, ycmp, z):
    """Calculate Stokes drift from spectrum

    :freq: frequency
    :spec: spectrum
    :xcmp: x-component
    :ycmp: y-component
    :z: vertical levels
    :returns: ustokes, vstokes

    """
    gravity = 9.81
    nfreq = np.size(freq)
    nlev = np.size(z)

    # calculate some common factors
    factor = np.zeros(nfreq, 1)
    factor2 = np.zeros(nfreq, 1)
    for i in range(nfreq):
        factor2[i] = 8.*np.pi*freq[i]**2./gravity
        factor[i] = 2.*np.pi*freq[i]*factor2[i]

    # initialize Stokes drift
    ustokes = np.zeros(nlev, 1)
    vstokes = np.zeros(nlev, 1)
    # integration
    for k in range(nlev):
        for i in range(nfreq):
            tmp = factor[i]*spec[i]*np.exp(factor2[i]*z[k])
            ustokes[k] += tmp*xcmp[i]
            vstokes[k] += tmp*ycmp[i]
        # add contribution from a f^-5 tail
        tmp = np.pi*freqs[nfreq]*factor[nfreq]*spec[nfreq]
             *(np.exp(factor2[nfreq]*z[k])-np.sqrt(np.pi*factor2[nfreq]*np.abs(z[k]))
             *(1.-math.erf(np.sqrt(factor2[nfreq]*np.abs(z[k])))))
        ustokes[k] += tmp*xcmp[i]
        vstokes[k] += tmp*ycmp[i]
    # return Stokes drift
    return ustokes, vstokes

if __name__ == "__main__":
    main()
