#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import argparse
import datetime
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2index
from gotmtool import nctime_to_datetime, nctime_indices

def main():
    # process the input arguments
    parser = argparse.ArgumentParser(description="""
            Test functions in $(prog)s.""")
    parser.add_argument('-i', '--input', action='store', dest='fname_in',
            metavar='NCFILENAME', required=True, help='Input netCDF filename')
    parser.add_argument('-o', '--output', action='store', dest='fname_out',
            metavar='DATFILENAME', required=True, help='Output filename')
    parser.add_argument('-ds', '--date_start', action='store', dest='date_start',
            metavar='STARTDATE',
            help='Starting date of input data, in the format of YYYYMMDD')
    parser.add_argument('-de', '--date_end', action='store', dest='date_end',
            metavar='ENDDATE',
            help='Ending date of input data, in the format of YYYYMMDD')
    parser.add_argument('--version', action='version', version='%(prog)s: 1.0')
    # parsing arguments and save to args
    args=parser.parse_args()

    #  TODO: Check if input data exist <17-12-17, Qing Li> #
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
    tdat = [dttime[i].isoformat(' ', 'seconds')
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
    # wave spectrum
    spec_h2 = spec*dfreq

    i = tidx_start
    # plot spectrum
    plot_spec(freq, spec[i,:])

    # depth
    z = np.arange(0,-20,-0.2)

    # Stokes drift
    ustokes, vstokes = stokes_drift(freq, spec_h2[i,:], xcmp[i,:], ycmp[i,:], z)
    ustokes_nt, vstokes_nt = stokes_drift(freq, spec_h2[i,:], xcmp[i,:],
            ycmp[i,:], z, tail=False)

    # plot Stokes drift
    plot_stokes_drift(z, ustokes, vstokes, 'stokesdrift.pdf')

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

def stokes_drift(freq, spec, xcmp, ycmp, z, tail=True):
    """Calculate Stokes drift from spectrum

    :freq: frequency
    :spec: spectrum
    :xcmp: x-component
    :ycmp: y-component
    :z: vertical levels
    :tail: (optional) not add contributions from spectral tail if set to False
    :returns: ustokes, vstokes

    """
    gravity = 9.81
    nfreq = np.size(freq)
    nlev = np.size(z)

    # calculate some common factors
    factor = np.zeros([nfreq, 1])
    factor2 = np.zeros([nfreq, 1])
    for i in range(nfreq):
        factor2[i] = 8.*np.pi*freq[i]**2./gravity
        factor[i] = 2.*np.pi*freq[i]*factor2[i]

    # initialize Stokes drift
    ustokes = np.zeros([nlev, 1])
    vstokes = np.zeros([nlev, 1])
    # integration
    for k in range(nlev):
        for i in range(nfreq):
            tmp = factor[i]*spec[i]*np.exp(factor2[i]*z[k])
            ustokes[k] += tmp*xcmp[i]
            vstokes[k] += tmp*ycmp[i]
    if tail:
        # add contribution from a f^-5 tail
        for k in range(nlev):
            tmp = np.pi*freq[-1]*factor[-1]*spec[-1]\
                *(np.exp(factor2[-1]*z[k])-np.sqrt(np.pi*factor2[-1]*np.abs(z[k]))\
                *(1.-math.erf(np.sqrt(factor2[-1]*np.abs(z[k])))))
            ustokes[k] += tmp*xcmp[i]
            vstokes[k] += tmp*ycmp[i]
    # return Stokes drift
    return ustokes, vstokes

if __name__ == "__main__":
    main()
