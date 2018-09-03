#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scipy.special as ss
import numpy as np
import argparse
import datetime
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2index
from gotmtool import nctime_to_datetime, nctime_indices

# constants
gravity = 9.81 # m s**-2
alpha_phil = 8.3e-3

def main():
    # process the input arguments
    parser = argparse.ArgumentParser(description="""
            Test functions in $(prog)s.""")
    parser.add_argument('-i', '--input', action='store', dest='fname_in',
            metavar='NCFILENAME', help='Input netCDF filename')
    parser.add_argument('-o', '--output', action='store', dest='fname_out',
            metavar='DATFILENAME', help='Output filename')
    parser.add_argument('-ds', '--date_start', action='store', dest='date_start',
            metavar='STARTDATE', help='Starting date of input data, in the format of YYYYMMDD')
    parser.add_argument('-de', '--date_end', action='store', dest='date_end',
            metavar='ENDDATE', help='Ending date of input data, in the format of YYYYMMDD')
    parser.add_argument('-t', '--test', action='store_true', dest='l_test',
            help='Run a few tests.')
    parser.add_argument('--version', action='version', version='%(prog)s: 1.0')
    # parsing arguments and save to args
    args=parser.parse_args()

    # test
    if args.l_test:
        print('Running some tests...')
        _test()
        exit(0)
    else:
        if any(opt is None for opt in [args.fname_in, args.fname_out, args.date_start, args.date_end]):
            # print help and exit
            parser.print_help()
            sys.exit(2)

    # set input
    fname_in = args.fname_in
    fname_out = args.fname_out
    date_start = args.date_start
    date_end = args.date_end

    # read data
    infile = Dataset(fname_in, 'r')
    # wave time
    nctime = infile.variables['waveTime']
    dttime = nctime_to_datetime(nctime) # time in datetime format
    # get starting and ending indices
    tidx_start, tidx_end = nctime_indices(nctime, date_start, date_end)
    tdat = [dttime[i].strftime('%Y-%m-%d %H:%M:%S')
            for i in range(tidx_start, tidx_end+1)] # truncated to seconds

    # band center frequency
    freq = infile.variables['waveFrequency'][:]
    # frequency bandwidth
    dfreq = infile.variables['waveBandwidth'][:]
    # band energy density
    spec = infile.variables['waveEnergyDensity'][:]
    # band mean direction that wave is coming from, in degree clockwise from the true North
    mdir = infile.variables['waveMeanDirection'][:]
    theta = 90.0-mdir # angle in degree counterclockwise from East
    d2r = np.pi/180.0
    xcmp = np.cos(theta*d2r)
    ycmp = np.sin(theta*d2r)
    # print(freq)
    # print(dfreq)

    # depth
    dz = 0.1
    nz = 100
    z = np.arange(0,-dz*nz,-dz)-dz/2
    # print(z)
    # print(z.shape)

    i = tidx_start
    # plot spectrum
    # plot_spec(freq, spec[i,:])

    # # Stokes drift
    # spec_h2 = spec*dfreq
    # ustokes, vstokes = stokes_drift_spec(freq, spec_h2[i,:], xcmp[i,:], ycmp[i,:], z)
    # ustokes_nt, vstokes_nt = stokes_drift_spec(freq, spec_h2[i,:], xcmp[i,:],
    #         ycmp[i,:], z, tail=False)
    # ustokes_int, vstokes_int = stokes_drift_spec(freq, spec_h2[i,:], xcmp[i,:],
    #         ycmp[i,:], z, cellavg=False)

    # # plot Stokes drift
    # plot_stokes_drift(z, ustokes, ustokes_nt)

###############################################################################
#                               Visualizations                                #
###############################################################################

def plot_spec(freq, spec, figname=None):
    """Plot spectrum

    :freq: frequency
    :spec: spectrum
    :figname: (optional) if present, save figure in a file
    :returns: None

    """
    plt.plot(freq, spec)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Energy density (m$^2$ Hz$^{-1}$)')
    xmin, xmax = plt.xlim()
    plt.xlim(0, xmax)
    ymin, ymax = plt.ylim()
    plt.ylim(0, ymax)
    if figname:
        plt.savefig(figname)
    else:
        plt.show()

def plot_stokes_drift(z, ustokes, vstokes, figname=None):
    """Plot profile data

    :z: depth
    :ustokes: Stokes drift, x-component
    :vstokes: Stokes drift, y-component
    :figname: (optional) if present, save figure in a file
    :returns: None

    """
    plt.plot(ustokes, z, '-k')
    plt.plot(vstokes, z, '--k')
    plt.xlabel('Stokes drift (m s$^{-1}$)')
    plt.ylabel('z (m)')
    if figname:
        plt.savefig(figname)
    else:
        plt.show()

###############################################################################
#                             Stokes drift, etc.                              #
###############################################################################

def stokes_drift_spec(freq, spec, xcmp, ycmp, z, tail=True, cellavg=True):
    """Calculate Stokes drift from spectrum. (Numerical integration)

    :freq: frequency
    :spec: spectrum*dfreq
    :xcmp: x-component
    :ycmp: y-component
    :z: vertical levels
    :tail: (optional) add tail contributions if true (default)
    :cellavg: (optional) cell averaged if true (default), interpolated otherwise
    :returns: ustokes, vstokes

    """
    nfreq = np.size(freq)
    nlev = np.size(z)

    # calculate some common factors
    const = 8.*np.pi**2./gravity
    factor2 = const*freq**2.
    factor = 2.*np.pi*freq*factor2

    # cutoff frequency
    freqc = 1.5*freq[-1]-0.5*freq[-2]
    dfreqc = freq[-1]-freq[-2]

    # initialize Stokes drift
    ustokes = np.zeros(nlev)
    vstokes = np.zeros(nlev)
    # integration
    if cellavg:
        # get cell
        zi = np.zeros(nlev+1)
        zi[1:nlev] = 0.5*(z[0:-1]+z[1:])
        zi[-1] = 2.*z[-1]-zi[-2]
        dz = zi[0:nlev]-zi[1:]
        for k in range(nlev):
            kdz = factor2*dz[k]/2.
            kflt = np.sinh(kdz)/kdz
            tmp = kflt*factor*spec*np.exp(factor2*z[k])
            ustokes[k] = np.sum(tmp*xcmp)
            vstokes[k] = np.sum(tmp*ycmp)
            if tail:
                # add contribution from a f^-5 tail
                aplus = max(1.e-9, -const*freqc**2.*zi[k])
                aminus = -const*freqc**2.*zi[k+1]
                iplus  = 2.*aplus/3.*(np.sqrt(np.pi*aplus)*ss.erfc(np.sqrt(aplus))\
                       -(1.-0.5/aplus)*np.exp(-aplus))
                iminus = 2.*aminus/3.*(np.sqrt(np.pi*aminus)*ss.erfc(np.sqrt(aminus))\
                       -(1.-0.5/aminus)*np.exp(-aminus))
                tmp = 2.*np.pi*freqc**2./dz[k]*spec[-1]/dfreqc*(iplus-iminus)
                ustokes[k] += tmp*xcmp[-1]
                vstokes[k] += tmp*ycmp[-1]
    else:
        for k in range(nlev):
            tmp = factor*spec*np.exp(factor2*z[k])
            ustokes[k] = np.sum(tmp*xcmp)
            vstokes[k] = np.sum(tmp*ycmp)
        if tail:
            # add contribution from a f^-5 tail
            factor2c = const*freqc**2
            factorc = 2.*np.pi*freqc*factor2c
            tmp = freqc*factorc*spec[-1]/dfreqc\
                *(np.exp(factor2c*z)-np.sqrt(np.pi*factor2c*np.abs(z))\
                *ss.erfc(np.sqrt(factor2c*np.abs(z))))
            ustokes += tmp*xcmp[-1]
            vstokes += tmp*ycmp[-1]
    # return Stokes drift
    return ustokes, vstokes

def stokes_drift_usp(freq, ussp, vssp, z, zi):
    """Calculate the Stokes drift profile from partitioned surface Stokes
    drift.

    :freq: frequency
    :ussp: partitioned surface Stokes drift, x-component
    :vssp: partitioned surface Stokes drift, y-component
    :z: vertical levels
    :returns: ustokes, vstokes

    """
    nfreq = np.size(freq)
    nlev = np.size(z)

    # calculate some common factors
    const = 8.*np.pi**2./gravity
    factor = const*freq**2.

    # initialize Stokes drift
    ustokes = np.zeros(nlev)
    vstokes = np.zeros(nlev)

    # get cell
    zi = np.zeros(nlev+1)
    zi[1:nlev] = 0.5*(z[0:-1]+z[1:])
    zi[-1] = 2.*z[-1]-zi[-2]
    dz = zi[0:nlev]-zi[1:]
    for k in range(nlev):
        kdz = factor*dz[k]/2.
        kflt = np.sinh(kdz)/kdz
        tmp = kflt*np.exp(factor*z[k])
        ustokes[k] = np.sum(tmp*ussp)
        vstokes[k] = np.sum(tmp*vssp)
    # return Stokes drift
    return ustokes, vstokes


###############################################################################
#                             Empirical spectra                               #
###############################################################################

def spectrum_Phillips(f, fp):
    """Returns the Phillips spectrum (Phillips, 1958) at given frequency.

    :f: frequency
    :fp: peak frequency
    :returns: spectrum at given frequency

    """
    nf = np.size(f)
    spec = np.zeros(nf)
    pi2 = 2.*np.pi
    spec[f>=fp] = alpha_phil*gravity**2./pi2**4./f[f>=fp]**5.
    return spec

def stokes_drift_Phillips(z, fp):
    """Returns the Stokes drift calculated from the Phillips spectrum
    (Analytical).

    :z: depth
    :fp: peak frequency
    :returns: Stokes drift

    """
    omegap = 2.*np.pi*fp
    us0 = 2.*alpha_phil*gravity/omegap
    kp = omegap**2./gravity
    T1 = np.exp(2.*kp*z)
    T2 = np.sqrt(2.*np.pi*kp*np.abs(z))*ss.erfc(np.sqrt(2.*kp*np.abs(z)))
    us = us0*(T1-T2)
    return us

###############################################################################
#                                  Utilities                                  #
###############################################################################

def get_delta(x):
    """Get the intervals of the input array.

    :x: a one-dimensional array
    :returns: an array of the intervals

    """
    nx = np.size(x)
    dx = np.zeros(nx)
    dx[1:-1]= 0.5*(x[2:]-x[0:-2])
    dx[0] = dx[1]
    dx[-1] = dx[-2]
    return dx

###############################################################################
#                                    Tests                                    #
###############################################################################

def _test():
    """A wrapper for all the tests.
    :returns: None

    """
    # Test 1: Test the integration along frequency in the function
    # stokes_drift_spec() with Phillips spectrum
    # high resolution versus low resolution in frequency
    _test_stokes_drift_freq(dz=0.1, dfreq=0.005, label='LRes')
    _test_stokes_drift_freq(dz=0.1, dfreq=0.0005, label='HRes')

def _test_stokes_drift_freq(dz=0.1, dfreq=0.005, label='LRes', outdir='./'):
    """Routine to test the integration along frequency in the function
    stokes_drift_spec() using Phillips spectrum.

    :dz: depth interval
    :dfreq: frequency interval
    :label: label
    :outdir: (optional) if specified, save figures to the directory
    :returns: None

    """
    # depth
    h = 10 # (m)
    nz = h/dz
    z = np.arange(0,-dz*nz,-dz)-dz/2
    # frequency
    f = np.arange(0.005, 0.5+dfreq, dfreq)
    df = get_delta(f)
    fp = 0.2
    nf = np.size(f)
    # directions
    xcmp_phil = np.ones(nf)
    ycmp_phil = np.zeros(nf)
    # get spectrum
    spec_phil = spectrum_Phillips(f, fp)
    # plot spectrum
    plt.figure()
    plt.plot(f, spec_phil, '-k')
    plt.step(f, spec_phil, where='mid', color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Energy density (m$^2$ Hz$^{-1}$)')
    xmin, xmax = plt.xlim()
    plt.xlim(0, xmax)
    ymin, ymax = plt.ylim()
    plt.ylim(0, ymax)
    suffix = label
    plt.savefig(outdir+'test_spectrum_'+suffix+'.pdf')
    # get Stokes drift
    us_phil = stokes_drift_Phillips(z, fp)
    us_cal, vs_cal = stokes_drift_spec(f, spec_phil*df, xcmp_phil, ycmp_phil, z, cellavg = False)
    # plot Stokes drift
    plt.figure()
    plt.plot(us_phil, z, '-k', label='Analytical')
    plt.plot(us_cal, z, '--r', label='Numerical')
    plt.legend(loc='lower right')
    plt.xlabel('Stokes drift (m s$^{-1}$)')
    plt.ylabel('z (m)')
    plt.savefig(outdir+'test_stokes_drift_'+suffix+'.pdf')

if __name__ == "__main__":
    main()
