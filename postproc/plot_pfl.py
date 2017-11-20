#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset, num2date
import argparse

def main():
    # process input arguments
    parser = argparse.ArgumentParser(description="""
        Plot time series of profile from GOTM together with the observation.""")
    parser.add_argument('-f', '--file', action='store', dest='fname',
            metavar='FILENAME', help='Input GOTM data')
    parser.add_argument('-v', '--variable', action='store', dest='vname',
            metavar='VARNAME', help='Variable name')
    parser.add_argument('-obs', '--obs', action='store_false', dest='lobs',
            help='Also plot data from the observation.')
    # parsing arguments and save to args
    args = parser.parse_args()

    # read data
    fin = Dataset(args.fname, 'r')

    nctime = fin.variables['time']
    z = fin.variables['z'][0,:,0,0]
    mdl = fin.variables[args.vname][:,:,0,0]
    if args.lobs:
        obs = fin.variables[args.vname+'_obs'][:,:,0,0]

    # nctime -> datetime
    dttime = num2date(nctime[:], units=nctime.units, calendar=nctime.calendar)

    # subplot, share x axis
    if args.lobs:
        fig, axarr = plt.subplots(2, sharex=True)
        im1 = axarr[0].contourf(dttime, z, np.transpose(obs))
        im2 = axarr[1].contourf(dttime, z, np.transpose(mdl))
        plt.ylabel('Depth (m)')
        plt.gcf().autofmt_xdate()
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im2, cax=cbar_ax)
    else:
        im = plt.contourf(dttime, z, np.transpose(mdl))
        plt.ylabel('Depth (m)')
        plt.colorbar(im)

    # save figure
    figname = args.vname+'.pdf'
    plt.savefig(figname)

if __name__ == "__main__":
    main()
