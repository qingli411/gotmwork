#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset, num2date
import argparse

def main():

    # process input arguments
    parser = argparse.ArgumentParser(description="""
        Plot time series from GOTM output.""")
    parser.add_argument('-f', '--file', action='store', dest='fname',
            metavar='FILENAME', help='Input GOTM data')
    parser.add_argument('-v', '--variable', action='store', dest='vname',
            metavar='VARNAME', nargs='+', help='Variable name')
    parser.add_argument('-o', '--output', action='store', dest='fname_out',
            metavar='FIGNAME', help='Output figure name')
    # parsing arguments and save to args
    args = parser.parse_args()

    # read data
    fin = Dataset(args.fname, 'r')

    # read time
    nctime = fin.variables['time']
    # nctime -> datetime
    dttime = num2date(nctime[:], units=nctime.units, calendar=nctime.calendar)

    # subplot, share x axis
    nv = len(args.vname)
    f, axarr = plt.subplots(nv, sharex=True)
    for i in range(nv):
        dat = fin.variables[args.vname[i]][:,0,0]
        axarr[i].plot(dttime, dat, '-k', linewidth=1.5)
    plt.gcf().autofmt_xdate()

    # save figure
    figname = args.fname_out
    plt.savefig(figname)

if __name__ == "__main__":
    main()
