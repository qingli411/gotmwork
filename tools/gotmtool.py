"""
Shared functions.

Qing Li, 20171213
"""

import datetime
from netCDF4 import num2date, date2index

def get_start_end_indices(nctime, date_start, date_end):
    """Returns indices corresponding to the starting and ending date & time.

    :nctime: (netCDF time object) nctime object
    :date_start: (str) starting date
    :date_end: (str) ending date
    :return: ([int,int]) starting and ending indices

    Returns 1 for the starting index if date_start is earlier than the first
    date & time in nctime, and dim_time for the ending index if date_end is
    later than the last date & time in nctime.

    """

    dim_time = len(nctime[:])
    # get time range indices
    if date_start and date_end:
        dt_start = datetime.datetime.strptime(date_start, '%Y%m%d')
        dt_end = datetime.datetime.strptime(date_end, '%Y%m%d')
        dt_in_start = num2date(nctime[1], units=nctime.units,
                calendar=nctime.calendar)
        dt_in_end = num2date(nctime[dim_time-1], units=nctime.units,
                calendar=nctime.calendar)
        # check starting time
        if dt_start < dt_in_start:
            tidx_start = 1  # skip the first data point
            date_start_str = datetime.datetime.strftime(dt_in_start, '%Y%m%d')
        else:
            tidx_start = date2index(dt_start, nctime)
            date_start_str = date_start
        # check ending time
        if dt_end > dt_in_end:
            tidx_end = dim_time
            date_end_str = datetime.datetime.strftime(dt_in_end, '%Y%m%d')
        else:
            tidx_end = date2index(dt_end, nctime)
            date_end_str = date_end
        # print some message
        print('Plotting time series from {} to {}.'
                .format(date_start_str, date_end_str))
    else:
        print('Plotting the full time series.')
        tidx_start = 1  # skip the first data point
        tidx_end = dim_time
    # return the indices
    return [tidx_start, tidx_end]
