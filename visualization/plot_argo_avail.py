import numpy as np
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap

def main():
    datadir = './'
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

    # date format
    dtformat = '%Y-%m-%d %H:%M:%S'
    dt_str_ref = '2008-06-01 00:00:00'
    i = 0
    with open(infile, 'r') as f:
        for line in f:
            [date, time, lon_str, lat_str, npfl_str, rlon_str, rlat_str] = line.split()
            lon[i] = float(lon_str)
            lat[i] = float(lat_str)
            rlon[i] = float(rlon_str)
            rlat[i] = float(rlat_str)
            npfl[i] = float(npfl_str)
            dt_str_tmp = date+' '+time
            dt_tmp = datetime.datetime.strptime(dt_str_tmp, dtformat)
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
    plot_map_scatter(rlon, rlat, npfl, 'npfl.png', vmax=40, vmin=0)

    plt.figure()
    plot_map_scatter(rlon, rlat, ndt, 'ndt.png', vmax=30, vmin=-30, cmap='RdBu')

def plot_map_scatter(rlon, rlat, dat, figname, vmax=None, vmin=None, cmap='rainbow'):
    """Plot scatters on a map

    :rlon: (Numpy array) 1D array of longitude
    :rlat: (Numpy array) 1D array of latitude
    :dat: (Numpy array) 1D array of data to plot
    :return: none
    """
    # plot map
    m = Basemap(projection='cyl', llcrnrlat=-70,urcrnrlat=70,llcrnrlon=0,urcrnrlon=360)
    # plot coastlines, draw label meridians and parallels.
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='gray',lake_color='lightgray')
    m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,1,1])
    m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,1,1])
    x, y = m(rlon, rlat)
    m.scatter(x, y, marker='.', s=18, c=dat, cmap=plt.cm.get_cmap(cmap), vmin=vmin, vmax=vmax)
    m.colorbar()
    # set figure size
    f = plt.gcf()
    f.set_size_inches(10, 5)
    # save figure
    plt.savefig(figname, dpi=300)

if __name__ == "__main__":
    main()
