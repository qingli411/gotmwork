"""
Shared functions.

Qing Li, 20171213
"""

import datetime
import numpy as np
from netCDF4 import num2date, date2index
from scipy.interpolate import griddata

###############################################################################
#                                  constants                                  #
###############################################################################

# gravitational acceleration (m/s^2)
g = 9.81
# reference density of seawater (kg/m^3)
rho_0 = 1027.0

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

def nctime_to_datetime(nctime, tidx_start=None, tidx_end=None, tidx_arr=None):
    """Convert from nctime object to datetime object.

    :nctime: (netCDF time object) nctime object
    :tidx_start: (int) starting index
    :tidx_end: (int) ending index
    :tidx_arr: (1D numpy array of int) indices
    :returns: (datetime object) datetime object

    """
    # check if attributes exist
    t_units = nctime.units
    try:
        t_cal = nctime.calendar
    except AttributeError :
        t_cal = 'standard'
    # return sliced datetime
    if tidx_arr is None:
        return num2date(nctime[tidx_start:tidx_end], units=t_units, calendar=t_cal)
    else:
        return num2date(nctime[tidx_arr], units=t_units, calendar=t_cal)

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
    # print_dttime_range(dttime)

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
#                                miscellaneous                                #
###############################################################################

def get_value_lat_lon(indata, lat2d, lon2d, rlat, rlon, imethod='nearest'):
    """Return the value of a variable at a given location (latitude and
    longitude).

    :indata: (2d/3d/4d numpy array) the value of variables in 2d/3d/4d array
                                    last two dimensions should be consistent
                                    with lat2d and lon2d
    :lat2d: (2d numpy array) latitude in 2d array
    :lon2d: (2d numpy array) longitude in 2d array
    :rlat: (float or 1d numpy array) output latitude
    :rlon: (float or 1d numpy array) output longitude
    :imethod: (str) interpolation method
    :returns: (float or 1d numpy array) value(s) of the variable at given
                                        latitude and longitude

    """
    # make sure lon2d and rlon are compatible
    maxlon = np.max(lon2d)
    minlon = np.min(lon2d)
    maxlat = np.max(lat2d)
    minlat = np.min(lat2d)
    rlon = float(rlon)
    rlat = float(rlat)
    assert rlat <= maxlat and rlat >= minlat, \
           'Value of rlat ({}) outside the range of lat2d ({} - {}). Stop.'.format( \
           rlat, minlat, maxlat)
    if rlon > maxlon:
        rlon = rlon - 360.0
    if rlon < minlon:
        rlon = rlon + 360.0
    assert rlon <= maxlon and rlon >= minlon, \
           'Value of rlon ({}) outside the range of lon2d ({} - {}). Stop.'.format( \
           rlon, minlon, maxlon)
    if np.min(lon2d) < 0.0 and rlon > 180.0:
        rlon = rlon - 360.0
    if np.max(lon2d) > 180.0 and rlon < 0.0:
        rlon = rlon + 360.0
    nsize = indata.ndim
    if imethod == 'nearest':
        ilat, ilon = get_index_lat_lon_nearest(lat2d, lon2d, rlat, rlon)
        if nsize == 2:
            outdata = indata[ilat,ilon]
        elif nsize == 3:
            outdata = indata[:,ilat,ilon]
        elif nsize == 4:
            outdata = indata[:,:,ilat,ilon]
        else:
            print('The variable {} has {} dimensions, not supported'
                    .format(vname, nsize))
            sys.exit(1)
    else:
        grid = (lat2d.flatten(),lon2d.flatten())
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

def get_index_lat_lon_nearest(lat2d, lon2d, rlat, rlon):
    """Return the indices of the point (rlon, rlat) in the array lon2d and lat2d

    :lat2d: (2d numpy array) latitude in 2d array
    :lon2d: (2d numpy array) longitude in 2d array
    :rlat: (float or 1d numpy array) output latitude
    :rlon: (float or 1d numpy array) output longitude
    :returns: indices

    """
    # make sure lon2d and rlon are compatible
    maxlon = np.max(lon2d)
    minlon = np.min(lon2d)
    maxlat = np.max(lat2d)
    minlat = np.min(lat2d)
    rlon = float(rlon)
    rlat = float(rlat)
    assert rlat <= maxlat and rlat >= minlat, \
           'Value of rlat ({}) outside the range of lat2d ({} - {}). Stop.'.format( \
           rlat, minlat, maxlat)
    if rlon > maxlon:
        rlon = rlon - 360.0
    if rlon < minlon:
        rlon = rlon + 360.0
    assert rlon <= maxlon and rlon >= minlon, \
           'Value of rlon ({}) outside the range of lon2d ({} - {}). Stop.'.format( \
           rlon, minlon, maxlon)
    grid = (lat2d.flatten(),lon2d.flatten())
    rlon = np.asarray(rlon)
    rlat = np.asarray(rlat)
    vlon = griddata(grid, lon2d.flatten(), (rlat, rlon), method='nearest')
    vlat = griddata(grid, lat2d.flatten(), (rlat, rlon), method='nearest')
    ilat0, ilon0 = np.where(lat2d == vlat)
    ilat1, ilon1 = np.where(lon2d == vlon)
    ilat = np.intersect1d(ilat0, ilat1)
    ilon = np.intersect1d(ilon0, ilon1)

    return ilat, ilon

#--------------------------------
# UNESCO equation of state for sea water
#--------------------------------

def unesco_eos_rho(S, T, p, pcorr=True):
    """Compute the in-situ density according to the UNESCO equation of state
    for sea water. The pressure correction can be switched off by pcorr=False.
    Borrowed from GOTM source code

    :S: (float or numpy array) salinity in psu
    :T: (float or numpy array) potential temperature in degC
    :p: (float or numpy array) pressure (absolute pressure - 1.01325 bar) in bar
    :pcorr: (logical) switch to turn on pressure correction
    :returns: density in kg/m^3

    """
    T2 = T*T
    T3 = T*T2
    T4 = T2*T2
    T5 = T*T4
    S15= S**1.5
    S2 = S*S
    S3 = S*S2

    x=999.842594+6.793952e-02*T-9.09529e-03*T2+1.001685e-04*T3
    x=x-1.120083e-06*T4+6.536332e-09*T5
    x=x+S*(0.824493-4.0899e-03*T+7.6438e-05*T2-8.2467e-07*T3)
    x=x+S*5.3875e-09*T4
    x=x+np.sqrt(S3)*(-5.72466e-03+1.0227e-04*T-1.6546e-06*T2)
    x=x+4.8314e-04*S2

    if pcorr:
        if np.any(p<0):
            raise ValueError ('p should be greater than 0')
        p2 = p*p
        K = 19652.21 \
            +148.4206     *T          -2.327105    *T2 \
            +  1.360477e-2*T3         -5.155288e-5 *T4 \
            +  3.239908      *p       +1.43713e-3  *T *p \
            +  1.16092e-4 *T2*p       -5.77905e-7  *T3*p \
            +  8.50935e-5    *p2      -6.12293e-6  *T *p2 \
            +  5.2787e-8  *T2*p2 \
            + 54.6746             *S  -0.603459    *T    *S \
            +  1.09987e-2 *T2     *S  -6.1670e-5   *T3   *S \
            +  7.944e-2           *S15+1.6483e-2   *T    *S15 \
            -  5.3009e-4  *T2     *S15+2.2838e-3      *p *S \
            -  1.0981e-5  *T *p   *S  -1.6078e-6   *T2*p *S \
            +  1.91075e-4    *p   *S15-9.9348e-7      *p2*S \
            +  2.0816e-8  *T *p2*S    +9.1697e-10  *T2*p2*S
        x=x/(1.-p/K)
    return x

def unesco_eos_alpha(S, T, p, pcorr):
    """Compute the thermal expansion coefficient

    :S: (float or numpy array) salinity in psu
    :T: (float or numpy array) potential temperature in degC
    :p: (float or numpy array) pressure (absolute pressure - 1.01325 bar) in bar
    :pcorr: (logical) switch to turn on pressure correction
    :returns: (float or numpy array) thermal expansion coefficient in 1/degC

    """
    delta = 0.01
    rho_a = unesco_eos_rho(S, T+0.5*delta, p, pcorr=pcorr)
    rho_b = unesco_eos_rho(S, T-0.5*delta, p, pcorr=pcorr)
    alpha = - (rho_a - rho_b) / (rho_0*delta)
    return alpha

def unesco_eos_beta(S, T, p, pcorr):
    """Compute the saline contraction coefficient

    :S: (float or numpy array) salinity in psu
    :T: (float or numpy array) potential temperature in degC
    :p: (float or numpy array) pressure (absolute pressure - 1.01325 bar) in bar
    :pcorr: (logical) switch to turn on pressure correction
    :returns: (float or numpy array) saline contraction coefficient in 1/psu

    """
    delta = 0.01
    rho_a = unesco_eos_rho(S+0.5*delta, T, p, pcorr=pcorr)
    rho_b = unesco_eos_rho(S-0.5*delta, T, p, pcorr=pcorr)
    beta = (rho_a - rho_b) / (rho_0*delta)
    return beta
