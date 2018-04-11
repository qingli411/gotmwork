"""
Shared functions.

Qing Li, 20171213
"""

import datetime
import numpy as np
from netCDF4 import num2date, date2index
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

###############################################################################
#                                preprocessing                                #
###############################################################################

#--------------------------------
# time
#--------------------------------

def nctime_indices(nctime, date_start, date_end):
    """Return indices corresponding to the starting and ending date & time.

    :nctime: (netCDF time object) nctime object
    :date_start: (str) starting date YYYYMMDD
    :date_end: (str) ending date YYYYMMDD
    :return: ([int,int]) starting and ending indices

    Returns [0, iend] if date_start is earlier than the first
    date & time in nctime, and [istart, ntime-1] if date_end is
    later than the last date & time in nctime. Returns [0, ntime-1]
    if no date_start or date_end is specified.

    """

    ntime = len(nctime[:])
    dtformat = '%Y%m%d'
    # get time range indices
    if date_start and date_end:
        dt_start = datetime.datetime.strptime(date_start, dtformat)
        dt_end = datetime.datetime.strptime(date_end, dtformat)
        tidx_start = date2index(dt_start, nctime, calendar=None, select='before')
        tidx_end = date2index(dt_end, nctime, calendar=None, select='after')
    else:
        tidx_start = 0
        tidx_end = ntime-1
    # return the indices
    return [tidx_start, tidx_end]

def nctime_to_datetime(nctime, tidx_start=None, tidx_end=None):
    """Convert from nctime object to datetime object.

    :nctime: (netCDF time object) nctime object
    :tidx_start: (int) starting index
    :tidx_end: (int) ending index
    :returns: (datetime object) datetime object

    """
    # check if attributes exist
    t_units = nctime.units
    try:
        t_cal = nctime.calendar
    except AttributeError :
        t_cal = 'standard'
    # return sliced datetime
    return num2date(nctime[tidx_start:tidx_end], units=t_units, calendar=t_cal)

def print_dttime_range(dttime):
    """Print the range of dttime.

    :dttime: (datetime object)
    :returns: none

    """
    dt_format = '%Y%m%d'
    str_start = datetime.datetime.strftime(dttime[0], dt_format)
    str_end   = datetime.datetime.strftime(dttime[-1],   dt_format)
    print('Timeseries from {} to {}...'.format(str_start, str_end))

def ncread_dim_time(infile, date_start, date_end):
    """Read the time dimension within the given range from
    a netCDF file.

    :infile: (netCDF4 Dataset) input netCDF file
    :date_start: (str) starting date, in the format of YYYYMMDD
    :date_end: (str) ending date, in the format of YYYYMMDD
    :returns: dttime, tidx_start, tidx_end
    ::dttime: (datetime obj) time in datetime format
    ::tidx_start: (int) starting index
    ::tidx_end: (int) ending index

    """
    # read time
    varlist = infile.variables.keys()
    if 'time' in varlist:
        nctime = infile.variables['time']
    elif 'TIME' in varlist:
        nctime = infile.variables['TIME']
    else:
        print('Time dimension is required and should have the name \"time\" or \"TIME\"')

    # get starting and ending indices
    tidx_start, tidx_end = nctime_indices(nctime, date_start, date_end)

    # nctime -> datetime
    dttime = nctime_to_datetime(nctime, tidx_start=tidx_start, tidx_end=tidx_end+1)

    # print some message
    print_dttime_range(dttime)

    # return
    return dttime, tidx_start, tidx_end

#--------------------------------
# write to .dat
#--------------------------------

def write_ts(fnout, tdat, vdat, mask=None):
    """Write time series in GOTMv5 format.

    :fnout: (str) filename of output file
    :tdat: (list) array of time
    :vdat: (list) array of variables
    :mask: (float, optional) value in vdat to skip
    :returns: none

    """
    nt = len(tdat)   # size of time
    with open(fnout, 'w') as fout:
        if mask is None:
            # no mask is applied
            for i in range(nt):
                # time
                out_str = '{}'.format(tdat[i])
                # variables
                for var in vdat:
                    out_str += '  {:10.6g}'.format(var[i])
                # newline
                out_str += '\n'
                fout.write(out_str)
        else:
            # skip if the value of any variable matches the mask value
            # or is NaN
            for i in range(nt):
                if not any(var[i] == mask or np.isnan(var[i]) for var in vdat):
                    # time
                    out_str = '{}'.format(tdat[i])
                    # variables
                    for var in vdat:
                        out_str += '  {:10.6g}'.format(var[i])
                    # newline
                    out_str += '\n'
                    fout.write(out_str)

def write_pfl(fnout, tdat, ddat, vdat, mask=None):
    """Write time series of profile in GOTMv5 format.

    :fnout: (str) filename of output file
    :tdat: (list) array of time
    :ddat: (list) array of depth
    :vdat: (list) array of variables
    :mask: (float, optional) value in vdat to skip
    :returns: none

    """
    nt = len(tdat[:]) # size of time
    nd = len(ddat[:]) # size of depth
    up_down = 2 # 1: data written from bottom to top (z<0 increasing)
                # otherwise: data written from top to bottom (z<0 decreasing)
    with open(fnout, 'w') as fout:
        if mask is None:
            # no mask is applied
            for i in range(nt):
                # time and dimension size
                out_str = '{}  {}  {}\n'.format(tdat[i], nd, up_down)
                fout.write(out_str)
                for j in range(nd):
                    # depth
                    out_str = '{:7.1f}'.format(ddat[j])
                    # variables
                    for var in vdat:
                        out_str += '  {:10.6g}'.format(var[i,j])
                    # newline
                    out_str += '\n'
                    fout.write(out_str)
        else:
            # skip the depth if the value of any variable matches the mask value
            # or is NaN
            for i in range(nt):
                fidx = []   # indices of filtered depth
                for j in range(nd):
                    if any(var[i,j] == mask or np.isnan(var[i,j]) for var in vdat):
                        fidx.append(j)
                nskip = len(fidx) # number of skipped depth
                if nd-nskip > 0:
                    # skip if there is no available data to write
                    # time and dimension size
                    out_str = '{}  {}  {}\n'.format(tdat[i], nd-nskip, up_down)
                    fout.write(out_str)
                    for j in range(nd):
                        if j not in fidx:
                            # depth
                            out_str = '{:7.1f}'.format(ddat[j])
                            # variables
                            for var in vdat:
                                out_str += '  {:10.6f}'.format(var[i,j])
                            # newline
                            out_str += '\n'
                            fout.write(out_str)

def write_spec(fnout, tdat, fdat, vdat):
    """Write spectra for GOMT input.

    :fnout: (str) filename of output file
    :tdat: (list) array of time
    :fdat: (list) array of frequencies
    :vdat: (list) array of variables
    :returns: none

    """
    nt = len(tdat[:]) # size of time
    nf = len(fdat[:]) # size of frequency
    up_down = 1 # 1: freq increasing with indicies
                # otherwise: freq decreasing with indicies
    with open(fnout, 'w') as fout:
        for i in range(nt):
            # time and dimension size
            out_str = '{}  {}  {}\n'.format(tdat[i], nf, up_down)
            fout.write(out_str)
            for j in range(nf):
                # frequencies
                out_str = '{:10.6g}'.format(fdat[j])
                # variables
                for var in vdat:
                    out_str += '  {:10.6g}'.format(var[i,j])
                # newline
                out_str += '\n'
                fout.write(out_str)

###############################################################################
#                               postprocessing                                #
###############################################################################

#--------------------------------
# read GOTM output
#--------------------------------

def ncread_pfl(ncvar, tidx_start=None, tidx_end=None):
    """Read in profile data [ntime, ndepth] from a netCDF variable.

    :ncvar: (netCDF variable) input variable
    :tidx_start: (int) starting index
    :tidx_end: (int) ending index
    :returns: (numpy array) profile data

    """
    # read in profile
    nsize = ncvar.ndim
    if nsize == 4:
        dat = ncvar[tidx_start:tidx_end,:,0,0]
    elif nsize == 2:
        dat = ncvar[tidx_start:tidx_end,:]
    else:
        dat = None
    return dat

def ncread_ts(ncvar, tidx_start=None, tidx_end=None):
    """Read in timeseries [ntime] from a netCDF variable.

    :ncvar: (netCDF variable) input variable
    :tidx_start: (int) starting index
    :tidx_end: (int) ending index
    :returns: (numpy array) timeseries data

    """
    # read in profile
    nsize = ncvar.ndim
    if nsize == 4:
        dat = ncvar[tidx_start:tidx_end,0,0,0]
    elif nsize == 3:
        dat = ncvar[tidx_start:tidx_end,0,0]
    elif nsize == 1:
        dat = ncvar[tidx_start:tidx_end]
    else:
        dat = None
    return dat

###############################################################################
#                                  analysis                                   #
###############################################################################

#--------------------------------
# driver
#--------------------------------

def do_analysis(name):
    """Driver to do the analysis

    :name: (str) name of the analysis
    :returns: scalar index of the analysis

    """
    switcher = {
            'mldMean': do_analysis_mld_mean
            }
    return switcher.get(name)

def do_analysis_mld_mean(infile, method='maxNsqr', tidx_start=None, tidx_end=None):
    """TODO: Docstring for do_analysis_mld_mean.

    :infile: (Dataset) netcdf file
    :method: (str) method to calculate the mixed layer depth
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (float) mean mixed layer depth

    """
    mld = get_mld(method)(infile, tidx_start=tidx_start, tidx_end=tidx_end)
    mld_mean = np.mean(mld)
    return mld_mean

#--------------------------------
# mixed layer depth
#--------------------------------

def get_mld(method):
    """Find the mixed layer depth

    :infile: (Dataset) netcdf file
    :method: (str) method to calculate the mixed layer depth
    :returns: (Numpy array) mixed layer depth

    """
    switcher = {
            'maxNsqr': get_mld_maxNsqr,
            'stratification': get_mld_maxNsqr,
            'deltaT': get_mld_deltaT,
            'temperature': get_mld_deltaT,
            'deltaR': get_mld_deltaR,
            'density': get_mld_deltaR
            }
    return switcher.get(method)

def get_mld_maxNsqr(infile, tidx_start=None, tidx_end=None):
    """Find the mixed layer depth defined as the depth where
       the stratification N^2 reaches its maximum

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (Numpy array) mixed layer depth

    """
    Nsqr = infile.variables['NN'][tidx_start:tidx_end,:,0,0]
    z = infile.variables['z'][tidx_start:tidx_end,:,0,0]
    nt = Nsqr.shape[0]
    mld = np.zeros(nt)
    # find the indices where N^2 reaches its maximum
    idx_max = np.argmax(Nsqr, 1)
    for i in np.arange(nt):
        mld[i] = z[i,idx_max[i]]
    return mld

def get_mld_deltaT(infile, deltaT=0.2, zRef=-10, tidx_start=None, tidx_end=None):
    """Find the mixed layer depth defined as the depth where
       the temperature difference from the reference level first exceed
       a threshold value

    :infile: (Dataset) netcdf file
    :deltaT: (float, optional) temperature threshold in degree C
    :zRef: (float, optional) depth of the reference level
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (Numpy array) mixed layer depth

    """
    Temp = infile.variables['temp'][tidx_start:tidx_end,:,0,0]
    z = infile.variables['z'][tidx_start:tidx_end,:,0,0]
    nt = Temp.shape[0]
    mld = np.zeros(nt)
    for i in np.arange(nt):
        idx_zref = np.argmin(np.abs(z[i,:]-zRef))
        dTemp = Temp[i,:]-Temp[i,idx_zref]
        # ignore the points above the reference level
        dTemp[idx_zref:] = 99.
        # find the maximum index (closest to the surface) where the temperature
        # difference is greater than the threshold value
        idx_min = np.max(np.where(dTemp<=-deltaT))
        mld[i] = z[i,idx_min]
    return mld

def get_mld_deltaR(infile, deltaR=0.03, zRef=-10, tidx_start=None, tidx_end=None):
    """Find the mixed layer depth defined as the depth where
       the potential density difference from the reference level first exceed
       a threshold value

    :infile: (Dataset) netcdf file
    :deltaR: (float, optional) potential density threshold in kg/m^3
    :zRef: (float, optional) depth of the reference level
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (Numpy array) mixed layer depth

    """
    Rho = infile.variables['rho'][tidx_start:tidx_end,:,0,0]
    z = infile.variables['z'][tidx_start:tidx_end,:,0,0]
    nt = Rho.shape[0]
    mld = np.zeros(nt)
    for i in np.arange(nt):
        idx_zref = np.argmin(np.abs(z[i,:]-zRef))
        dRho = Rho[i,:]-Rho[i,idx_zref]
        # ignore the points above the reference level
        dRho[idx_zref:] = -99.
        # find the maximum index (closest to the surface) where the density
        # difference is greater than the threshold value
        idx_min = np.max(np.where(dRho>=deltaR))
        mld[i] = z[i,idx_min]
    return mld

###############################################################################
#                                visualization                                #
###############################################################################

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

###############################################################################
#                                miscellaneous                                #
###############################################################################

def get_value_lat_lon(indata, lat2d, lon2d, rlat, rlon, imethod='nearest'):
    """Return the value of a variable at a given location (latitude and
    longitude).

    :lat2d: (2d Numpy array) latitude in 2d array
    :lon2d: (2d Numpy array) longitude in 2d array
    :indata: (2d/3d/4d Numpy array) the value of variables in 2d/3d/4d array
                                    last two dimensions should be consistent
                                    with lat2d and lon2d
    :rlat: (float or 1d Numpy array) output latitude
    :rlon: (float or 1d Numpy array) output longitude
    :imethod: (str) interpolation method
    :returns: (float or 1d Numpy array) value(s) of the variable at given
                                        latitude and longitude

    """
    #  TODO: no wrapping of longitude at lon=360, can be done by expanding
    #        lat2d and lon2d <07-03-18, Qing Li> #
    grid = (lat2d.flatten(),lon2d.flatten())
    nsize = indata.ndim
    rlat = np.asarray(rlat)
    rlon = np.asarray(rlon)
    if nsize == 2:
        outdata = griddata(grid, indata.flatten(), (rlat, rlon), method=imethod)
    elif nsize == 3:
        nt = indata.shape[0]
        nout = rlat.size
        outdata = np.zeros((nt, nout))
        for i in np.arange(nt):
            outdata[i,:] = griddata(grid, indata[i,:,:].flatten(), (rlat, rlon), method=imethod)
    elif nsize == 4:
        nt = indata.shape[0]
        nd = indata.shape[1]
        nout = rlat.size
        outdata = np.zeros((nt, nd, nout))
        for i in np.arange(nt):
            for j in np.arange(nd):
                outdata[i,j,:] = griddata(grid, indata[i,j,:,:].flatten(),
                        (rlat, rlon), method=imethod)
    else:
        print('The variable {} has {} dimensions, not supported'
                .format(vname, nsize))
        sys.exit(1)

    # return the interpolated data
    return outdata
