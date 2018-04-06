import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot_map_scatter(rlon, rlat, dat, figname, vmax=None, vmin=None, cmap='rainbow'):
    """Plot scatters on a map

    :rlon: (Numpy array) 1D array of longitude
    :rlat: (Numpy array) 1D array of latitude
    :dat: (Numpy array) 1D array of data to plot
    :return: none
    """
    # plot map
    m = Basemap(projection='cyl', llcrnrlat=-72,urcrnrlat=72,llcrnrlon=0,urcrnrlon=360)
    # plot coastlines, draw label meridians and parallels.
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='gray',lake_color='lightgray')
    m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,1,1])
    m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,1,1])
    x, y = m(rlon, rlat)
    m.scatter(x, y, marker='.', s=32, c=dat, cmap=plt.cm.get_cmap(cmap), vmin=vmin, vmax=vmax)
    m.colorbar()
    # set figure size
    f = plt.gcf()
    f.set_size_inches(8, 4)
    # save figure
    plt.savefig(figname, dpi=300)

