import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
from netCDF4 import Dataset
from gotmplot import plot_map_scatter

def main():
    # forcing method
    forcmethod = 'COREII'
    # turbulent method
    turbmethod = 'KPP-CVMix'
    # starting and ending date (yyyymmdd), empty if find all
    date_start = '20080601'
    date_end = '20081231'

    # data directory
    datadir = os.environ['HOME']+'/scratch'
    outdir = os.environ['HOME']+'/work/gotmfigures/'+forcmethod
    # netcdf file name
    filename = 'gotm_out.nc'

    # format of case name to match
    # casefmt = r'^\w+_LAT-{0,1}\d{1,2}_LON\d{1,3}_\w+-{0,1}\w*_\d{8}-\d{8}$'
    if date_start.strip() and date_end.strip():
        datefmt = r'^\d{8}$'
        if re.search(datefmt, date_start) and re.search(datefmt, date_end):
            casefmt = r'^'+forcmethod+r'_LAT-{0,1}\d{1,2}_LON\d{1,3}_'+turbmethod+r'_'+date_start+'-'+date_end+r'$'
        else:
            print('The starting date and ending date should be in yyyymmdd format, got {} and {}. Stop.'.format(date_start, date_end))
            sys.exit(1)
    else:
        casefmt = r'^'+forcmethod+r'_LAT-{0,1}\d{1,2}_LON\d{1,3}_'+turbmethod+r'_\d{8}-\d{8}$'

    # print(casefmt)

    # search for directories that match the case name pattern
    lslist = os.listdir(datadir)

    r = re.compile(casefmt)
    caselist = list(filter(r.search, lslist))
    ncase = len(caselist)
    print('Total number of cases: {}'.format(ncase))

    # initialize array
    lat = np.zeros(ncase)
    lon = np.zeros(ncase)
    dat = np.zeros(ncase)

    # loop over all cases
    for i in np.arange(ncase):
        ncfile = datadir+'/'+caselist[i]+'/'+filename
        if os.path.isfile(ncfile):
            infile = Dataset(ncfile, 'r')
            lat[i] = infile.variables['lat'][:]
            lon[i] = infile.variables['lon'][:]
            tmp = infile.variables['sst'][:]
            # processing the data
            dat[i] = np.mean(tmp)
        else:
            lat[i] = None
            lon[i] = None
            dat[i] = None

    # update longitude in the range (0, 360)
    lon[lon<0] = lon[lon<0]+360

    # plot figures
    plt.figure()
    outfig = outdir+'/sst.png'
    plot_map_scatter(lon, lat, dat, outfig, vmax=30, vmin=0)

if __name__ == "__main__":
    main()
