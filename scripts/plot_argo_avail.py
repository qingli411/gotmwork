import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
from gotmtool import plot_map_scatter

def main():
    datadir = os.environ['HOME']+'/work/gotmfigures/argo_availability/'
    filename = 'profile_map.dat'
    infile = datadir+filename

    # count line number
    nlines = sum(1 for line in open(infile, 'r'))

    # initialize array
    npfl = np.zeros(nlines)
    lat = np.zeros(nlines)
    lon = np.zeros(nlines)
    rlat = np.zeros(nlines)
    rlon = np.zeros(nlines)
    ndt = np.zeros(nlines)
    temp = np.zeros(nlines)

    # date format
    dtformat = '%Y-%m-%d %H:%M:%S'
    dt_str_ref = '2008-06-01 00:00:00'
    i = 0
    with open(infile, 'r') as f:
        for line in f:
            [date, time, lon_str, lat_str, npfl_str, rlon_str, rlat_str, temp_str] = line.split()
            lon[i] = float(lon_str)
            lat[i] = float(lat_str)
            rlon[i] = float(rlon_str)
            rlat[i] = float(rlat_str)
            npfl[i] = float(npfl_str)
            temp[i] = float(temp_str)
            dt_str_tmp = date+' '+time
            dt_tmp = datetime.datetime.strptime(dt_str_tmp, dtformat).replace(year=2008)
            dt_ref = datetime.datetime.strptime(dt_str_ref, dtformat)
            ndt[i] = (dt_tmp-dt_ref).days
            i += 1

    # print some message
    print('Total number grid points: {}'.format(nlines))
    print('Average number of available profiles per grid: {:6.2f}'.format(np.mean(npfl)))
    print('Root mean square time difference in days: {:6.2f}'.format(np.sqrt(np.mean(ndt**2))))

    # update longitude in the range (0, 360)
    rlon[rlon<0] = rlon[rlon<0]+360

    # plot figures
    plt.figure()
    outfig = datadir+'/npfl.png'
    plot_map_scatter(rlon, rlat, npfl, vmax=100, vmin=0)
    plt.savefig(outfig, dpi = 300)

    plt.figure()
    outfig = datadir+'/ndt.png'
    plot_map_scatter(rlon, rlat, ndt, vmax=30, vmin=-30, cmap='RdBu')
    plt.savefig(outfig, dpi = 300)

    plt.figure()
    outfig = datadir+'/temp.png'
    plot_map_scatter(rlon, rlat, temp, vmax=30, vmin=-2)
    plt.savefig(outfig, dpi = 300)

if __name__ == "__main__":
    main()
