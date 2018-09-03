import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from netCDF4 import Dataset, num2date
tooldir_default = os.environ['GOTMWORK_ROOT']+'/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import *

def main():
    # case name
    # casename = 'OSMOSIS_winter_OSMOSIS'
    casename = 'COREII_LAT-46_LON222_20080601-20080731'
    forcmethod = 'COREII'

    # data directory
    datadir = os.environ['GOTMRUN_ROOT']+'/VR5m_DT600s/KPP-CVMix'
    outdir = os.environ['GOTMFIG_ROOT']+'/'+forcmethod
    # netcdf file name
    filename = 'gotm_out.nc'

    # loop over all cases
    ncfile = datadir+'/'+casename+'/'+filename
    infile = Dataset(ncfile, 'r')
    lat = infile.variables['lat'][:]
    lon = infile.variables['lon'][:]
    # update longitude in the range (0, 360)
    lon[lon<0] = lon[lon<0]+360
    # read time
    t_cal = 'standard'
    nctime = infile.variables['time']
    dttime = num2date(nctime[:], units=nctime.units, calendar=t_cal)
    z = infile.variables['z'][0,:,0,0]
    fld = infile.variables['rho'][:,:,0,0]

    # mld1 = get_mld('maxNsqr')(infile)
    # mld2 = get_mld('deltaT')(infile)
    # mld3 = get_mld('deltaR')(infile)

    mld1 = get_mld('stratification')(infile)
    mld2 = get_mld('temperature')(infile)
    mld3 = get_mld('density')(infile)
    # tidx_start = 1050
    # tidx_end = 1500
    # mld_ts = get_mld('maxNsqr')(infile, tidx_start=tidx_start, tidx_end=tidx_end+1)
    # mld_mean = do_analysis('mldMean')(infile, tidx_start=tidx_start, tidx_end=tidx_end+1)
    mld_ts = get_mld('maxNsqr')(infile)
    mld_mean = do_analysis('mldMean')(infile)
    print(mld_ts.shape)
    # print(tidx_end-tidx_start+1)
    print('Timeseries of MLD = ')
    print(mld_ts)
    print('Mean MLD = {}'.format(mld_mean))
    print(mld3)

    # plot figures

    # contour levels
    c_max = 30
    c_min = 0
    c_int = (c_max-c_min)/20
    levels = np.arange(c_min, c_max+c_int, c_int)

    # contourf by default
    im = plt.contourf(dttime, z, np.transpose(fld), extend='both', cmap='rainbow')
    cb = plt.colorbar(im)
    plt.ylabel('Depth (m)')
    plt.plot(dttime, mld1, '-', color='white')
    plt.plot(dttime, mld2, '--', color='gray')
    plt.plot(dttime, mld3, ':', color='black')

    # auto adjust the x-axis label
    plt.gcf().autofmt_xdate()

    # save figure
    # figname = 'test.png'
    # plt.savefig(figname, dpi = 300)
    plt.show()

if __name__ == "__main__":
    main()
