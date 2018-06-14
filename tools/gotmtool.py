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
#                                  constants                                  #
###############################################################################

# gravitational acceleration (m/s^2)
g = 9.81
# specific heat of seawater (J/kg/degC)
cp = 3985.0
# Von Karman constant
kappa = 0.4
# reference density of seawater (kg/m^3)
rho_0 = 1027.0
# reference salinity (psu)
S_0 = 35.0
# reference temperature (degC)
T_0 = 10.0
# constant thermal expansion coefficient (1/degC)
alpha_0 = 1.65531e-4
# constant saline contraction coefficient (1/psu)
beta_0 = 7.59494e-4

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
#                                  analysis                                   #
###############################################################################

#--------------------------------
# driver
#--------------------------------

def do_analysis(name):
    """Driver to do the analysis

    :name: (str) name of the analysis
    :returns: (funciton) corresponding function

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

    :method: (str) method to calculate the mixed layer depth
    :returns: (function) corresponding function

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
    :returns: (numpy array) mixed layer depth

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
    :returns: (numpy array) mixed layer depth

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
        # ignore nan
        dTemp[np.isnan(dTemp)] = 99.
        # find the maximum index (closest to the surface) where the temperature
        # difference is greater than the threshold value
        idxlist = np.where(dTemp<=-deltaT)
        if idxlist[0].size>0:
            idx_min = np.max(idxlist)
            mld[i] = z[i,idx_min]
        else:
            mld[i] = np.min(z[i,:])
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
    :returns: (numpy array) mixed layer depth

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
        # ignore nan
        dRho[np.isnan(dRho)] = -99.
        # find the maximum index (closest to the surface) where the density
        # difference is greater than the threshold value
        idxlist = np.where(dRho>=deltaR)
        if idxlist[0].size>0:
            idx_min = np.max(idxlist)
            mld[i] = z[i,idx_min]
        else:
            mld[i] = np.min(z[i,:])
    return mld

#--------------------------------
# derived variables
#--------------------------------

def get_variable(method):
    """Find the derived variable

    :method: (str) variable name
    :returns: (function) corresponding function

    """
    switcher = {
            'LaTurb': get_la_turb,
            'LaSL': get_la_sl,
            'hoLmo': get_h_over_lmo,
            'buoyancy': get_buoyancy,
            'spice': get_spice,
            'dPEdt': get_dpedt
            }
    return switcher.get(method)

def get_la_turb(infile, tidx_start=None, tidx_end=None):
    """Find the turbulent Langmuir number defined as
       La_t = sqrt{u^*/u^S}

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) Langmuir number

    """
    # friction velocity
    ustar = infile.variables['u_taus'][tidx_start:tidx_end,0,0]
    # surface Stokes drift
    usx = infile.variables['u0_stokes'][tidx_start:tidx_end,0,0]
    usy = infile.variables['v0_stokes'][tidx_start:tidx_end,0,0]
    us = np.sqrt(usx**2.+usy**2.)
    # calculate Langmuir number
    la = np.sqrt(ustar/us)
    return la

def get_la_sl(infile, tidx_start=None, tidx_end=None):
    """Find the surface layer averaged Langmuir number defined as
       La_{SL} = sqrt{u^*/<u^S>_{SL}}, where
       <u^S>_{SL} = int_{-0.2h_b}^0 u^S(z) dz

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) Langmuir number

    """
    # friction velocity
    ustar = infile.variables['u_taus'][tidx_start:tidx_end,0,0]
    # Stokes drift profiles
    ustokes = infile.variables['u_stokes'][tidx_start:tidx_end,:,0,0]
    vstokes = infile.variables['v_stokes'][tidx_start:tidx_end,:,0,0]
    zi = infile.variables['zi'][tidx_start:tidx_end,:,0,0]
    h = infile.variables['h'][tidx_start:tidx_end,:,0,0]
    # boundary layer depth
    hbl = get_mld('deltaR')(infile, tidx_start=tidx_start, tidx_end=tidx_end)
    # surface layer: upper 20% of the boundary layer
    hsl = 0.2*hbl
    # loop over time to calculate the surface layer averaged Stokes drift
    # note that zi has indices 0:nlev whereas z has indices 0:nlev-1, this is
    # different from the indices in GOTM
    nt = ustar.shape[0]
    ussl = np.zeros(nt)
    vssl = np.zeros(nt)
    for i in np.arange(nt):
        ihsl = np.argmin(zi[i,:]<hsl[i])
        dzb = zi[i,ihsl]-hsl[i]
        ussl[i] = np.sum(ustokes[i,ihsl:]*h[i,ihsl:])+ustokes[i,ihsl-1]*dzb
        vssl[i] = np.sum(vstokes[i,ihsl:]*h[i,ihsl:])+vstokes[i,ihsl-1]*dzb
    ussl = ussl/np.abs(hsl)
    vssl = vssl/np.abs(hsl)
    # surface layer averaged Langmuir number
    la = np.sqrt(ustar/np.sqrt(ussl**2.+vssl**2.))
    return la

def get_h_over_lmo(infile, tidx_start=None, tidx_end=None):
    """Find the stability parameter defined as h/L
       where h is the boundary layer depth
       and L is the Monin–Obukhov length

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) stability parameter

    """
    # get boundary layer depth
    hbl   = get_mld('deltaR')(infile, tidx_start=tidx_start, tidx_end=tidx_end)
    # friction velocity
    ustar = infile.variables['u_taus'][tidx_start:tidx_end,0,0]
    # surface temperature and salinity
    temp0 = infile.variables['temp'][tidx_start:tidx_end,-1,0,0]
    salt0 = infile.variables['salt'][tidx_start:tidx_end,-1,0,0]
    # surface temperature flux
    tflux = infile.variables['heat'][tidx_start:tidx_end,0,0]/cp/rho_0
    # correction for solar radiation
    rad   = infile.variables['rad'][tidx_start:tidx_end,:,0,0]
    z     = infile.variables['z'][tidx_start:tidx_end,:,0,0]
    nt    = ustar.shape[0]
    rflux = np.zeros(nt)
    for i in np.arange(nt):
        ihbl = np.argmin(np.abs(z[i,:]-hbl[i]))
        rflux[i] = rad[i,-1]-rad[i,ihbl]
    tflux = tflux+rflux/cp/rho_0
    # surface salinity flux
    sflux = -(infile.variables['precip'][tidx_start:tidx_end,0,0]
          + infile.variables['evap'][tidx_start:tidx_end,0,0])*salt0
    # surface buoyancy flux (positive for stable condition)
    bflux = g*alpha_0*tflux-g*beta_0*sflux
    # Monin-Obukhov length
    Lmo = ustar**3.0/kappa/bflux
    # filter out zeros
    Lmo = np.ma.array(Lmo, mask=(Lmo==0))
    # h over L
    hoL = -abs(hbl)/Lmo
    return hoL

def get_buoyancy(infile, tidx_start=None, tidx_end=None):
    """Calculate the buoyancy from temperature and salinity
    assuming linear equation of state

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) buoyancy

    """
    # temperature and salinity
    temp  = infile.variables['temp'][tidx_start:tidx_end,:,0,0]
    salt  = infile.variables['salt'][tidx_start:tidx_end,:,0,0]
    # buoyancy
    buoy  = g*alpha_0*(temp-T_0)-g*beta_0*(salt-S_0)
    return buoy

def get_spice(infile, tidx_start=None, tidx_end=None):
    """Calculate the spice from temperature and salinity
    assuming linear equation of state

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) spice

    """
    # temperature and salinity
    temp  = infile.variables['temp'][tidx_start:tidx_end,:,0,0]
    salt  = infile.variables['salt'][tidx_start:tidx_end,:,0,0]
    # spice
    spice  = g*alpha_0*(temp-T_0)+g*beta_0*(salt-S_0)
    return spice

def get_dpedt(infile, tidx_start=None, tidx_end=None):
    """Calculate the rate of change in the total potential energy (PE)

    :infile: (Dataset) netcdf file
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) rate of change in PE

    """
    # time series of potential energy
    epot = infile.variables['Epot'][tidx_start:tidx_end,0,0]
    # time (sec)
    time = infile.variables['time'][tidx_start:tidx_end]
    # get the time derivative
    nt = time.shape[0]
    dpedt = np.zeros(nt)
    dpedt[1:-1] = (epot[2:]-epot[0:-2])/(time[2:]-time[0:-2])
    dpedt[0] = (epot[1]-epot[0])(time[1]-time[0])
    dpedt[-1] = (epot[-1]-epot[-2])/(time[-1]-time[-2])
    return dpedt

###############################################################################
#                                visualization                                #
###############################################################################

def plot_map_scatter(rlon, rlat, dat, figname, vmax=None, vmin=None, cmap='rainbow'):
    """Plot scatters on a map

    :rlon: (numpy array) 1D array of longitude
    :rlat: (numpy array) 1D array of latitude
    :dat: (numpy array) 1D array of data to plot
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

    :lat2d: (2d numpy array) latitude in 2d array
    :lon2d: (2d numpy array) longitude in 2d array
    :indata: (2d/3d/4d numpy array) the value of variables in 2d/3d/4d array
                                    last two dimensions should be consistent
                                    with lat2d and lon2d
    :rlat: (float or 1d numpy array) output latitude
    :rlon: (float or 1d numpy array) output longitude
    :imethod: (str) interpolation method
    :returns: (float or 1d numpy array) value(s) of the variable at given
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
