import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap

#--------------------------------
# Constants
#--------------------------------

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

#--------------------------------
# GOTMProfile
#--------------------------------
class GOTMProfile(object):

    """GOTMProfile object"""

    def __init__(self, time=None, z=None, data=None, name=None):
        """Initialize GOTMProfile

        :time: (1D numpy array/datetime object) time
        :z: (1D numpy array) vertical coordinate
        :data: (2D numpy array) data at each time and z
        :name: (str) name of variable

        """
        self.time = time
        self.z = z
        self.data = data
        self.name = name
        if data:
            self.data_mean = np.mean(data, axis=0)
        else:
            self.data_mean = None

    def plot(self, axis=None, xlim=None, ylim=None,
                   xlabel=None, ylabel=None, title=None,
                   ptype='contourf', **kwargs):
        """Plot the Hovmoller diagram (time - z)

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'Time' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :ptype: (str, optional) plot type, valid values: contourf (default), pcolor
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.contourf() or
                                       matplotlib.pyplot.pcolor() depending on ptype
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if not axis:
            axis = plt.gca()
        # plot type
        if ptype == 'contourf':
            fig = axis.contourf(self.time, self.z, np.transpose(self.data), **kwargs)
        elif ptype == 'pcolor':
            fig = axis.pcolor(self.time, self.z, np.transpose(self.data), **kwargs)
        else:
            raise ValueError('Plot type (ptype) should be \'contourf\' or \'pcolor\', got {}.'.format(ptype))
        # x- and y-label, turn off by passing in 'off'
        if not xlabel:
            axis.set_xlabel('Time')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if not ylabel:
            axis.set_ylabel('Depth (m)')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim:
            axis.set_xlim(xlim)
        if ylim:
            axis.set_ylim(ylim)
        # return figure
        return fig

    def plot_mean(self, axis=None, xlim=None, ylim=None,
                        xlabel=None, ylabel=None, title=None, **kwargs):
        """Plot the mean profile

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'self.name' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'Depth (m)' by default, 'off' to turn it off
        :title: (str, optional) title
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.plot()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if not axis:
            axis = plt.gca()
        # plot figure
        fig = plt.plot(self.data_mean, self.z, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if not xlabel:
            axis.set_xlabel(self.name)
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if not ylabel:
            axis.set_ylabel('Depth (m)')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim:
            axis.set_xlim(xlim)
        if ylim:
            axis.set_ylim(ylim)
        # return figure
        return fig


#--------------------------------
# GOTMTimeseries
#--------------------------------
class GOTMTimeseries(object):

    """GOTMTimeseries object"""

    def __init__(self, time=None, data=None, name=None):
        """Initialize GOTMTimeseries

        :time: (1D numpy array/datetime object) time
        :data: (1D numpy array) data at each location
        :name: (str) name of variable

        """

        self.time = time
        self.data = data
        self.name = name
        self.mean_data = np.mean(self.data, axis=0)

    def plot(self, axis=None, xlim=None, ylim=None,
                   xlabel=None, ylabel=None, title=None, **kwargs):
        """Plot timeseries

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :xlim: ([float, float], optional) upper and lower limits of the x-axis
        :ylim: ([float, float], optional) upper and lower limits of the y-axis
        :xlabel: (str, optional) x-label, 'Time' by default, 'off' to turn it off
        :ylabel: (str, optional) y-label, 'self.name' by default, 'off' to turn it off
        :title: (str, optional) title
        :**kwargs: (keyword arguments) to be passed to matplotlib.pyplot.plot()
        :returns: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if not axis:
            axis = plt.gca()
        # plot figure
        fig = axis.plot(self.time, self.data, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if not xlabel:
            axis.set_xlabel('Time')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if not ylabel:
            axis.set_ylabel(self.name)
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim:
            axis.set_xlim(xlim)
        if ylim:
            axis.set_ylim(ylim)
        # return figure
        return fig


#--------------------------------
# GOMTMap
#--------------------------------

class GOTMMap(object):

    """GOTMMap object"""

    def __init__(self, data=None, lon=None, lat=None, name=None, units=None):
        """Initialize GOTMMap

        :data: (1D numpy array) data at each location
        :lon: (1D numpy array) longitude
        :lat: (1D numpy array) latitude
        :name: (str) name of variable
        :units: (str) units of variable

        """
        self.data = data
        self.lon = lon
        self.lat = lat
        self.name = name
        self.units = units

    def save(self, path):
        """Save GOTMMap object

        :path: (str) path of file to save
        :returns: none

        """
        np.savez(path, data=self.data, lon=self.lon, lat=self.lat, name=self.name, units=self.units)

    def load(self, path):
        """Load data to GOTMMap object

        :path: (str) path of file to load
        :returns: (GOTMMap object)

        """
        dat = np.load(path)
        self.__init__(data=dat['data'], lon=dat['lon'], lat=dat['lat'],
                name=str(dat['name']), units=str(dat['units']))
        return self

    def plot(self, axis=None, levels=None, add_colorbar=True, cmap='rainbow', **kwargs):
        """Plot scatters on a map

        :axis: (matplotlib.axes, optional) axis to plot figure on
        :leveles: (list, optional) list of levels
        :add_colorbar: (bool) do not add colorbar if False
        :cmap: (str, optional) colormap
        :**kwargs: (keyword arguments) to be passed to mpl_toolkits.basemap.scatter()
        :return: (matplotlib figure object) figure

        """
        # use curret axis if not specified
        if not axis:
            axis = plt.gca()
        # plot map
        m = Basemap(projection='cyl', llcrnrlat=-72, urcrnrlat=72, llcrnrlon=0, urcrnrlon=360, ax=axis)
        # m = Basemap(projection='cyl', llcrnrlat=-72, urcrnrlat=72, llcrnrlon=0, urcrnrlon=360)
        # plot coastlines, draw label meridians and parallels.
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.fillcontinents(color='gray',lake_color='lightgray')
        m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,1,1])
        m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,1,1])
        x, y = m(self.lon, self.lat)
        # manually mapping levels to the colormap if levels is passed in,
        # otherwise linear mapping
        if levels:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            fig = m.scatter(x, y, marker='.', s=32, c=self.data, norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            fig = m.scatter(x, y, marker='.', s=32, c=self.data, cmap=plt.cm.get_cmap(cmap), **kwargs)
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-2, 2))
            cb.update_ticks()
            if self.units:
                cb.ax.set_title(self.units)
        return fig


#--------------------------------
# GOTMOutputData
#--------------------------------

class GOTMOutputData(object):

    """GOTM Output data object"""

    def __init__(self, path, init_time_location=True):
        """Initialize the data

        :path: (str) path to the GOTM output netCDF file
        :init_time_location: (bool) do not initialize time and lon etc. if
                                    False, which is used when constructing large
                                    dataset to reduce initialization time.

        """
        # path of data
        self._path = path
        # open the dataset
        self.open()
        # time and location if requested
        if init_time_location:
            # latitude
            self.lat = self.dataset.variables['lat'][:]
            # Coriolis parameter
            self.f = 4.*np.pi/86400*np.sin(self.lat*np.pi/180.)
            # longitude
            self.lon = self.dataset.variables['lon'][:]
            # time
            time = self.dataset.variables['time']
            self.time = time[:]
            self.time_units = time.units
            self.time_calendar = time.calendar
        # list of dimensions
        self.list_dimensions = list(self.dataset.dimensions.keys())
        # list of variables
        list_variables = list(self.dataset.variables.keys())
        self.list_variables = [x for x in list_variables if x not in self.list_dimensions]
        # list of time series and profiles
        self.list_timeseries = []
        self.list_profile = []
        for var in self.list_variables:
            ndim = self.dataset.variables[var].ndim
            if ndim == 4:
                self.list_profile.append(var)
            elif ndim == 3:
                self.list_timeseries.append(var)
            else:
                pass
        # update list of variables to include the derived ones
        self.list_profile = self.list_profile + self._get_derived_profile(list_keys=True)
        self.list_timeseries = self.list_timeseries + self._get_derived_timeseries(list_keys=True)
        # close dataset
        self.close()

    def open(self):
        """Open netCDF dataset

        """
        # netCDF4 Dataset
        self.dataset = Dataset(self._path, 'r')

    def close(self):
        """Close netCDF datset

        """
        if self.dataset.isopen():
            self.dataset.close()

    def read_profile(self, var, tidx_start=None, tidx_end=None, ignore_time=False, **kwargs):
        """Return profile variable and z (fixed in time)

        :var: (str) variable name
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :ignore_time: (bool) ignore time dimension if True
        :**kwargs: (keyword arguments) to be passed to _get_derived_profile()
        :returns: (GOTMProfile object) profile

        """
        # open dataset
        self.open()
        # read variable
        if var in self.list_variables:
            dat = self.dataset.variables[var][tidx_start:tidx_end,:,0,0]
            coord = self.dataset.variables[var].coordinates
            if 'zi' in coord:
                z = self.dataset.variables['zi'][0,:,0,0]
            elif 'z' in coord:
                z = self.dataset.variables['z'][0,:,0,0]
            else:
                raise AttributeError('z coordinate not fould.')
        else:
            dat, z = self._get_derived_profile(var)(tidx_start=tidx_start, tidx_end=tidx_end, **kwargs)
        # create GOTMProfile object
        if not ignore_time:
            time = num2date(self.time[tidx_start:tidx_end], units=self.time_units, calendar=self.time_calendar)
        else:
            time = None
        out = GOTMProfile(time=time, z=z, data=dat, name=var)
        # close data
        self.close()
        return out

    def read_timeseries(self, var, tidx_start=None, tidx_end=None, ignore_time=False, **kwargs):
        """Return timeseries of variable [var]

        :var: (str) variable name
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :ignore_time: (bool) ignore time dimension if True
        :**kwargs: (keyword arguments) to be passed to _get_derived_timeseries()
        :returns: (GOTMTimeseries object) time series

        """
        # open dataset
        self.open()
        # read variable
        if var in self.list_variables:
            dat = self.dataset.variables[var][tidx_start:tidx_end,0,0]
        else:
            dat = self._get_derived_timeseries(var)(tidx_start=tidx_start, tidx_end=tidx_end, **kwargs)
        # create GOTMTimeseries object
        if not ignore_time:
            time = num2date(self.time[tidx_start:tidx_end], units=self.time_units, calendar=self.time_calendar)
        else:
            time = None
        out = GOTMTimeseries(time=time, data=dat, name=var)
        # close data
        self.close()
        return out

    def _get_derived_profile(self, name=None, list_keys=False):
        """Find the derived profile variable

        :name: (str) variable name
        :list_keys: (bool) simply return the list of keys if True
        :returns: (function/str) corresponding function or list of keys

        """
        switcher = {
                'buoy': self._get_buoyancy,
                'buoyancy': self._get_buoyancy,
                'spice': self._get_spice
                }
        if list_keys:
            return list(switcher.keys())
        elif name in switcher.keys():
            return switcher.get(name)
        else:
            raise ValueError('Variable \'{}\' not found.'.format(name))

    def _get_buoyancy(self, tidx_start=None, tidx_end=None):
        """Calculate the buoyancy from temperature and salinity
        assuming linear equation of state

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) buoyancy

        """
        # temperature and salinity
        temp  = self.dataset.variables['temp'][tidx_start:tidx_end,:,0,0]
        salt  = self.dataset.variables['salt'][tidx_start:tidx_end,:,0,0]
        # buoyancy
        buoy  = g*alpha_0*(temp-T_0)-g*beta_0*(salt-S_0)
        # z
        z     = self.dataset.variables['z'][0,:,0,0]
        return buoy, z

    def _get_spice(self, tidx_start=None, tidx_end=None):
        """Calculate the spice from temperature and salinity
        assuming linear equation of state

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) spice

        """
        # temperature and salinity
        temp  = self.dataset.variables['temp'][tidx_start:tidx_end,:,0,0]
        salt  = self.dataset.variables['salt'][tidx_start:tidx_end,:,0,0]
        # spice
        spice  = g*alpha_0*(temp-T_0)+g*beta_0*(salt-S_0)
        # z
        z     = self.dataset.variables['z'][0,:,0,0]
        return spice, z

    def _get_derived_timeseries(self, name=None, list_keys=False):
        """Find the derived variable

        :name: (str) variable name
        :list_keys: (bool) simply return the list of keys if True
        :returns: (function/str) corresponding function or list of keys

        """
        switcher = {
                'LaTurb': self._get_la_turb,
                'LaSL': self._get_la_sl,
                'hoLmo': self._get_h_over_lmo,
                'bflux': self._get_bflux,
                'dPEdt': self._get_dpedt,
                'mixEf1': self._get_dpedt_over_ustar3,
                'mld_maxNsqr': self._get_mld_maxNsqr,
                'mld_deltaT': self._get_mld_deltaT,
                'mld_deltaR': self._get_mld_deltaR,
                'PE': self._get_pez
                }
        if list_keys:
            return list(switcher.keys())
        elif name in switcher.keys():
            return switcher.get(name)
        else:
            raise ValueError('Variable \'{}\' not found.'.format(name))

    def _get_la_turb(self, tidx_start=None, tidx_end=None):
        """Find the turbulent Langmuir number defined as
           La_t = sqrt{u^*/u^S}

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) Langmuir number

        """
        # friction velocity
        ustar = self.dataset.variables['u_taus'][tidx_start:tidx_end,0,0]
        # surface Stokes drift
        usx = self.dataset.variables['u0_stokes'][tidx_start:tidx_end,0,0]
        usy = self.dataset.variables['v0_stokes'][tidx_start:tidx_end,0,0]
        us = np.sqrt(usx**2.+usy**2.)
        # calculate Langmuir number
        la = np.sqrt(ustar/us)
        return la

    def _get_la_sl(self, tidx_start=None, tidx_end=None):
        """Find the surface layer averaged Langmuir number defined as
           La_{SL} = sqrt{u^*/<u^S>_{SL}}, where
           <u^S>_{SL} = int_{-0.2h_b}^0 u^S(z) dz

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) Langmuir number

        """
        # friction velocity
        ustar = self.dataset.variables['u_taus'][tidx_start:tidx_end,0,0]
        # Stokes drift profiles
        ustokes = self.dataset.variables['u_stokes'][tidx_start:tidx_end,:,0,0]
        vstokes = self.dataset.variables['v_stokes'][tidx_start:tidx_end,:,0,0]
        zi = self.dataset.variables['zi'][tidx_start:tidx_end,:,0,0]
        h = self.dataset.variables['h'][tidx_start:tidx_end,:,0,0]
        # boundary layer depth
        hbl = self._get_derived_timeseries('mld_deltaR')(tidx_start=tidx_start, tidx_end=tidx_end)
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

    def _get_bflux(self, tidx_start=None, tidx_end=None):
        """Find the surface buoyancy flux

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) surface buoyancy flux

        """
        # get boundary layer depth
        hbl   = self._get_derived_timeseries('mld_deltaR')(tidx_start=tidx_start, tidx_end=tidx_end)
        # surface temperature and salinity
        temp0 = self.dataset.variables['temp'][tidx_start:tidx_end,-1,0,0]
        salt0 = self.dataset.variables['salt'][tidx_start:tidx_end,-1,0,0]
        # surface temperature flux
        tflux = self.dataset.variables['heat'][tidx_start:tidx_end,0,0]/cp/rho_0
        # correction for solar radiation
        rad   = self.dataset.variables['rad'][tidx_start:tidx_end,:,0,0]
        z     = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
        nt    = temp0.shape[0]
        rflux = np.zeros(nt)
        for i in np.arange(nt):
            ihbl = np.argmin(np.abs(z[i,:]-hbl[i]))
            rflux[i] = rad[i,-1]-rad[i,ihbl]
        tflux = tflux+rflux/cp/rho_0
        # surface salinity flux
        sflux = -(self.dataset.variables['precip'][tidx_start:tidx_end,0,0]
                + self.dataset.variables['evap'][tidx_start:tidx_end,0,0])*salt0
        # surface buoyancy flux (positive for stable condition)
        bflux = g*alpha_0*tflux-g*beta_0*sflux
        return bflux

    def _get_h_over_lmo(self, tidx_start=None, tidx_end=None):
        """Find the stability parameter defined as h/L
           where h is the boundary layer depth
           and L is the Monin–Obukhov length

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) stability parameter

        """
        # get boundary layer depth
        hbl   = self._get_derived_timeseries('mld_deltaR')(tidx_start=tidx_start, tidx_end=tidx_end)
        # get surface buoyancy flux
        bflux = self._get_bflux(tidx_start=tidx_start, tidx_end=tidx_end)
        # friction velocity
        ustar = self.dataset.variables['u_taus'][tidx_start:tidx_end,0,0]
        # Monin-Obukhov length
        Lmo = ustar**3.0/kappa/bflux
        # filter out zeros
        Lmo = np.ma.array(Lmo, mask=(Lmo==0))
        # h over L
        hoL = -abs(hbl)/Lmo
        return hoL

    def _get_dpedt(self, tidx_start=None, tidx_end=None):
        """Calculate the rate of change in the total potential energy (PE)

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) rate of change in PE

        """
        # time series of potential energy
        epot = self.dataset.variables['Epot'][tidx_start:tidx_end,0,0]
        # time (sec)
        time = self.dataset.variables['time'][tidx_start:tidx_end]
        # get the time derivative
        nt = time.shape[0]
        dpedt = np.zeros(nt)
        dpedt[1:-1] = (epot[2:]-epot[0:-2])/(time[2:]-time[0:-2])
        dpedt[0] = (epot[1]-epot[0])/(time[1]-time[0])
        dpedt[-1] = (epot[-1]-epot[-2])/(time[-1]-time[-2])
        return dpedt

    def _get_dpedt_over_ustar3(self, tidx_start=None, tidx_end=None):
        """Calculate the bulk rate of change in the total potential
           energy (PE), normalized by firction

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) rate of change in PE

        """
        # rate of potential energy
        dpedt = self._get_dpedt(tidx_start=tidx_start, tidx_end=tidx_end)
        # friction velocity
        ustar = self.dataset.variables['u_taus'][tidx_start:tidx_end,0,0]
        # normalize by rho_0*ustar**3
        bulk_dpedt = dpedt/ustar**3/rho_0
        return bulk_dpedt

    def _get_mld_maxNsqr(self, tidx_start=None, tidx_end=None):
        """Find the mixed layer depth defined as the depth where
           the stratification N^2 reaches its maximum

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) mixed layer depth

        """
        Nsqr = self.dataset.variables['NN'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
        nt = Nsqr.shape[0]
        mld = np.zeros(nt)
        # find the indices where N^2 reaches its maximum
        idx_max = np.argmax(Nsqr, 1)
        for i in np.arange(nt):
            mld[i] = z[i,idx_max[i]]
        return mld

    def _get_mld_deltaT(self, tidx_start=None, tidx_end=None, deltaT=0.2, zRef=-10):
        """Find the mixed layer depth defined as the depth where
           the temperature difference from the reference level first exceed
           a threshold value

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :deltaT: (float, optional) temperature threshold in degree C
        :zRef: (float, optional) depth of the reference level
        :returns: (numpy array) mixed layer depth

        """
        Temp = self.dataset.variables['temp'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
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

    def _get_mld_deltaR(self, tidx_start=None, tidx_end=None, deltaR=0.03, zRef=-10):
        """Find the mixed layer depth defined as the depth where
           the potential density difference from the reference level first exceed
           a threshold value

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :deltaR: (float, optional) potential density threshold in kg/m^3
        :zRef: (float, optional) depth of the reference level
        :returns: (numpy array) mixed layer depth

        """
        Rho = self.dataset.variables['rho'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
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

    def _get_pez(self, tidx_start=None, tidx_end=None, depth=None):
        """Calculate the total potential energy (PE) integrated from
           surface to a certain depth

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :depth: (float, optional) index of depth
        :returns: (numpy array) rate of change in PE

        """

        # find layer thickness
        if 'h' in self.list_variables:
            h = self.dataset.variables['h'][tidx_start:tidx_end,:,0,0]
        else:
            zi = self.dataset.variables['zi'][tidx_start:tidx_end,:,0,0]
            h = zi[:,1:] - zi[:,0:-1]
        # time series of buoyancy
        prfl = self.read_profile('buoy', tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True)
        # find zidx
        if depth:
            zidx = np.argmin(np.abs(prfl.z-depth))
        else:
            zidx = 0
        # compute potential energy
        tmp = prfl.data[:,zidx:]*h[:,zidx:]*(prfl.z[zidx]-prfl.z[zidx:])*rho_0
        pe = np.sum(tmp, axis=1)
        pe = pe-pe[0]
        # get the time derivative
        return pe


#--------------------------------
# GOTMOutputDataSet
#--------------------------------

class GOTMOutputDataSet(object):

    """GOTMOutputDataSet object, a set of GOTMOutputData object"""

    def __init__(self, paths, keys=None):
        """Initialize GOTMOutputDataSet

        :paths: (list of str) paths to the GOTM output netCDF file
        :keys: (list of str, optional) keys to refer to each of the cases in the dataset

        """
        # paths of dataset
        self._paths = paths
        # number of cases
        self.ncase = len(self._paths)
        # set case name
        if keys:
            nkeys = len(keys)
            if nkeys != self.ncase:
                raise ValueError('The number of keys is not equal to the number of cases.')
            self.casename = keys
        else:
            self.casename = [str(i) for i in range(self.ncase)]
        # construct dictionary for cases
        self.cases = dict([(self.casename[i], GOTMOutputData(self._paths[i])) for i in range(self.ncase)])
        # case name of the reference case
        self.ref_casename = self.casename[0]

    def diff_profile(self, var, cname, ref_cname=None, tidx_start=None, tidx_end=None):
        """The difference of profile between two cases

        :var: (str) variable name
        :cname: (str) case name
        :ref_cname: (str) reference case name
        :tidx_start: (int, optional) starting index for time
        :tidx_end: (int, optional) ending index for time
        :returns: (GOTMProfile object) difference in [var]

        """
        # reference case
        if not ref_cname:
            ref_cname = self.ref_casename
        # read profiles
        prfl  = self.cases[ref_cname].read_profile(var, tidx_start=tidx_start, tidx_end=tidx_end)
        data0 = prfl.data
        time0 = prfl.time
        z0    = prfl.z
        prfl  = self.cases[cname].read_profile(var, tidx_start=tidx_start, tidx_end=tidx_end)
        data1 = prfl.data
        time1 = prfl.time
        z1    = prfl.z
        # check consistency of coordinates
        if any(time1[i] != time0[i] for i in range(len(time0))):
            raise ValueError('Length of time not consistent between two cases.')
        if len(z1) != len(z0) or any(z1[i] != z0[i] for i in range(len(z0))):
            print('z-coordinates not consistent between two cases, interpolating to that of the reference case.')
            data1 = np.array([np.interp(z0, z1, data1[i,:]) for i in range(len(time0))])
        # create GOTMProfile object
        out = GOTMProfile(time=time0, z=z0, data=data1-data0, name=var)
        return out


#--------------------------------
# GOTMOutputDataMap
#--------------------------------

class GOTMOutputDataMap(object):

    """GOTMOutputDataMap object, a large set of GOTMOutputData object that share
       the same time- and z-dimensions but at different locations"""

    def __init__(self, paths):
        """Initialize the GOTMOutputDataMap

        :paths: (list of str) paths to the GOTM output netCDF file

        """
        # paths of dataset
        self._paths = paths
        # number of cases
        self.ncase = len(self._paths)
        # loading meta data from the first data file, assuming the same meta
        # data across all data files in the dataset
        firstdata = GOTMOutputData(self._paths[0])
        # time
        self.time = firstdata.time
        self.ntime = self.time.size
        self.time_units = firstdata.time_units
        self.time_calendar = firstdata.time_calendar
        # list of dimensions
        self.list_dimensions = firstdata.list_dimensions
        # list of variables
        self.list_variables = firstdata.list_variables
        # list of time series and profiles
        self.list_timeseries = firstdata.list_timeseries
        self.list_profile = firstdata.list_profile
        # latitude and longitude
        self.lat, self.lon = self._init_latlon()

    def _init_latlon(self):
        """Initialize arrays for lat and lon.

        :returns: (numpy array, numpy array) latitude and longitude

        """
        lat = np.zeros(self.ncase)
        lon = np.zeros(self.ncase)
        for i in range(self.ncase):
            dataset = Dataset(self._paths[i], 'r')
            lat[i] = dataset.variables['lat'][:]
            lon[i] = dataset.variables['lon'][:]
        return lat, lon

    def mean_state_profile(self, var, tidx_start=None, tidx_end=None, zidx_start=None, zidx_end=None):
        """Return the mean state of a profile variable

        :var: (str) variable name
        :tidx_start: (int, optional) starting index for time
        :tidx_end: (int, optional) ending index for time
        :zidx_start: (int, optional) starting index for z
        :zidx_end: (int, optional) ending index z
        :returns: (GOTMMap object) mean state

        """
        mdat = np.zeros(self.ncase)
        for i in range(self.ncase):
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            prfl = tmp.read_profile(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True)
            mdat[i] = np.mean(np.mean(prfl.data, axis=0)[zidx_start:zidx_end], 0)
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        return out

    def mean_state_timeseries(self, var, tidx_start=None, tidx_end=None):
        """Return the mean state of a timeseries variable

        :var: (str) variable name
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (GOTMMap object) mean state

        """
        mdat = np.zeros(self.ncase)
        for i in range(self.ncase):
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            ts = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True)
            mdat[i] = ts.mean_data
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        return out

    def delta_timeseries(self, var, tidx_start=None, tidx_end=None):
        """Return the net change of a timeseries variable

        :var: (str) variable name
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (GOTMMap object) mean state

        """
        mdat = np.zeros(self.ncase)
        for i in range(self.ncase):
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            ts = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True)
            mdat[i] = ts.data[-1] - ts.data[0]
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        return out

