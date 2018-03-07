"""
Shared functions.

Qing Li, 20171213
"""

import datetime
import numpy as np
from netCDF4 import num2date, date2index

def nctime_indices(nctime, date_start, date_end):
    """Returns indices corresponding to the starting and ending date & time.

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
        tidx_start, tidx_end = date2index([dt_start, dt_end], nctime,
                calendar=None, select='nearest')
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
    # get time slice
    if tidx_start is None or tidx_end is None:
        istart = 0
        iend = len(nctime[:])
    else:
        istart = tidx_start
        iend = tidx_end+1
    # check if attribute exist
    t_units = nctime.units
    try:
        t_cal = nctime.calendar
    except AttributeError : # gregorian if attribute doesn't exist
        t_cal = 'gregorian'
    # return sliced datetime
    return num2date(nctime[istart:iend], units=t_units, calendar=t_cal)

def print_dttime_range(dttime):
    """Print the range of dttime.

    :dttime: (datetime object)
    :returns: none

    """
    dt_format = '%Y%m%d'
    str_start = datetime.datetime.strftime(dttime[0], dt_format)
    str_end   = datetime.datetime.strftime(dttime[-1],   dt_format)
    print('Timeseries from {} to {}...'.format(str_start, str_end))

def ncread_pfl(ncvar, tidx_start=None, tidx_end=None):
    """Read in profile data [ntime, ndepth] from a netCDF variable.

    :ncvar: (netCDF variable) input variable
    :tidx_start: (int) starting index
    :tidx_end: (int) ending index
    :returns: (numpy array) profile data

    """
    # get time slice
    if tidx_start is None or tidx_end is None:
        istart = 0
        iend = -1
    else:
        istart = tidx_start
        iend = tidx_end+1
    # read in profile
    nsize = ncvar.ndim
    if nsize == 4:
        dat = ncvar[istart:iend,:,0,0]
    elif nsize == 2:
        dat = ncvar[istart:iend,:]
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
    # get time slice
    if tidx_start is None or tidx_end is None:
        istart = 0
        iend = -1
    else:
        istart = tidx_start
        iend = tidx_end+1
    # read in profile
    nsize = ncvar.ndim
    if nsize == 4:
        dat = ncvar[istart:iend,0,0,0]
    elif nsize == 3:
        dat = ncvar[istart:iend,0,0]
    elif nsize == 1:
        dat = ncvar[istart:iend]
    else:
        dat = None
    return dat

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
                    out_str += '  {:6.2f}'.format(var[i])
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
                        out_str += '  {:6.2f}'.format(var[i])
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
                        out_str += '  {:10.6f}'.format(var[i,j])
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

