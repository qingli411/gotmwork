import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.colors import from_levels_and_colors
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid

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
        try:
            self.data_mean = np.mean(data, axis=0)
        except:
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
        if axis is None:
            axis = plt.gca()
        # plot type
        if ptype == 'contourf':
            fig = axis.contourf(self.time, self.z, np.transpose(self.data), **kwargs)
        elif ptype == 'pcolor':
            fig = axis.pcolor(self.time, self.z, np.transpose(self.data), **kwargs)
        else:
            raise ValueError('Plot type (ptype) should be \'contourf\' or \'pcolor\', got {}.'.format(ptype))
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel('Time')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel('Depth (m)')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
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
        if axis is None:
            axis = plt.gca()
        # plot figure
        fig = plt.plot(self.data_mean, self.z, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel(self.name)
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel('Depth (m)')
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
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
        try:
            self.data_mean = np.mean(data, axis=0)
        except TypeError:
            self.data_mean = None

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
        if axis is None:
            axis = plt.gca()
        # plot figure
        fig = axis.plot(self.time, self.data, **kwargs)
        # x- and y-label, turn off by passing in 'off'
        if xlabel is None:
            axis.set_xlabel('Time')
        else:
            if xlabel != 'off':
                axis.set_xlabel(xlabel)
        if ylabel is None:
            axis.set_ylabel(self.name)
        else:
            if ylabel != 'off':
                axis.set_ylabel(ylabel)
        # x- and y-limits
        if xlim is not None:
            axis.set_xlim(xlim)
        if ylim is not None:
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

    def __neg__(self):
        """Return the negated GOTMMap object

        """
        out = GOTMMap(data=-self.data, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        return out

    def __add__(self, other):
        """Add 'other' to a GOTMMap object

        :other: (float, int or GOTMMap object) object to be added

        """
        if isinstance(other, float) or isinstance(other, int):
            out = GOTMMap(data=self.data+other, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        elif isinstance(other, GOTMMap):
            assert self.name == other.name, "GOTMMap object to be added has a different name."
            assert self.units == other.units, "GOTMMap object to be added has a different unit."
            data = np.zeros(self.data.size)
            loc_self = list(zip(self.lon, self.lat))
            loc_other = list(zip(other.lon, other.lat))
            for idx, val in enumerate(loc_self):
                if val in loc_other:
                    idx_other = loc_other.index(val)
                    data[idx] = self.data[idx] + other.data[idx_other]
                else:
                    data[idx] = np.nan
            out = GOTMMap(data=data, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        else:
            raise TypeError('Addition is not defined between a GOTMMap object and a {} object'.format(type(other)))
        return out

    def __sub__(self, other):
        """Subtract 'other' from a GOTMMap object

        :other: (float, int, or GOTMMap object) object to be subtracted

        """
        if isinstance(other, float) or isinstance(other, int):
            out = GOTMMap(data=self.data-other, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        elif isinstance(other, GOTMMap):
            assert self.name == other.name, "GOTMMap object to be subtracted has a different name."
            assert self.units == other.units, "GOTMMap object to be subtracted has a different unit."
            data = np.zeros(self.data.size)
            loc_self = list(zip(self.lon, self.lat))
            loc_other = list(zip(other.lon, other.lat))
            for idx, val in enumerate(loc_self):
                if val in loc_other:
                    idx_other = loc_other.index(val)
                    data[idx] = self.data[idx] - other.data[idx_other]
                else:
                    data[idx] = np.nan
            out = GOTMMap(data=data, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        else:
            raise TypeError('Subtraction is not defined between a GOTMMap object and a {} object'.format(type(other)))
        return out

    def __mul__(self, other):
        """Multiply a GOTMMap object by 'other'

        :other: (float, int, or GOTMMap object) object to be multiplied

        """
        if isinstance(other, float) or isinstance(other, int):
            out = GOTMMap(data=self.data*other, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        elif isinstance(other, GOTMMap):
            data = np.zeros(self.data.size)
            loc_self = list(zip(self.lon, self.lat))
            loc_other = list(zip(other.lon, other.lat))
            for idx, val in enumerate(loc_self):
                if val in loc_other:
                    idx_other = loc_other.index(val)
                    data[idx] = self.data[idx] - other.data[idx_other]
                else:
                    data[idx] = np.nan
            out = GOTMMap(data=data, lon=self.lon, lat=self.lat, name=self.name+' * '+other.name, units=self.units+' * '+other.units)
        else:
            raise TypeError('Multiplication is not defined between a GOTMMap object and a {} object'.format(type(other)))
        return out

    def __truediv__(self, other):
        """Divide a GOTMMap object by 'other'

        :other: (float, int, or GOTMMap object) object to be divided

        """
        if isinstance(other, float) or isinstance(other, int):
            out = GOTMMap(data=self.data/other, lon=self.lon, lat=self.lat, name=self.name, units=self.units)
        elif isinstance(other, GOTMMap):
            data = np.zeros(self.data.size)
            loc_self = list(zip(self.lon, self.lat))
            loc_other = list(zip(other.lon, other.lat))
            for idx, val in enumerate(loc_self):
                if val in loc_other:
                    idx_other = loc_other.index(val)
                    data[idx] = self.data[idx] / other.data[idx_other]
                else:
                    data[idx] = np.nan
            if other.name == self.name:
                name = 'Ratio of '+self.name
                units = 'None'
            else:
                name = self.name+' / '+other.name
                units = self.units+' / '+other.units
            out = GOTMMap(data=data, lon=self.lon, lat=self.lat, name=name, units=units)
        else:
            raise TypeError('Division is not defined between a GOTMMap object and a {} object'.format(type(other)))
        return out

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
        dat = np.load(path, allow_pickle=True)
        self.__init__(data=dat['data'], lon=dat['lon'], lat=dat['lat'],
                name=str(dat['name']), units=str(dat['units']))
        return self

    def masked(self, mask, mask_data=np.nan):
        """Apply mask to GOTMMap object. The mask should also be a GOTMMap object,
           with 1 for valid and 0 for invalid.

        :mask: (GOTMMap object) mask, 1 for valid, 0 for invalid
        :mask_data: (optional) values to be filled in maked points
        :return: (GOTMMap object) masked GOTMMap

        """
        if mask.data.size != self.data.size:
            raise ValueError('The dimension of mask does not match.')
        dat = self.data
        self.data = np.where(mask.data==0, mask_data, dat)

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
        if axis is None:
            axis = plt.gca()
        # plot map
        m = Basemap(projection='cyl', llcrnrlat=-72, urcrnrlat=72, llcrnrlon=20, urcrnrlon=380, ax=axis)
        # plot coastlines, draw label meridians and parallels.
        m.drawcoastlines()
        m.drawmapboundary(fill_color='lightgray')
        m.fillcontinents(color='gray',lake_color='lightgray')
        m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1])
        m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,0,1])
        data = self.data
        lat = self.lat
        lon = self.lon
        # shift longitude
        lon = np.where(lon < 20., lon+360., lon)
        x, y = m(lon, lat)
        # manually mapping levels to the colormap if levels is passed in,
        # otherwise linear mapping
        if levels is not None:
            bounds = np.array(levels)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
            fig = m.scatter(x, y, marker='.', s=32, c=data, norm=norm, cmap=plt.cm.get_cmap(cmap), **kwargs)
        else:
            fig = m.scatter(x, y, marker='.', s=32, c=data, cmap=plt.cm.get_cmap(cmap), **kwargs)
        # add colorbar
        if add_colorbar:
            cb = m.colorbar(fig, ax=axis)
            cb.formatter.set_powerlimits((-2, 2))
            cb.update_ticks()
        return fig

    def zonal_mean(self):
        """Calculate the zonal mean.

        :returns: (numpy array) array of latitude and zonal mean data

        """
        lat_all = self.lat
        lat, counts = np.unique(lat_all, return_counts=True)
        nlat = lat.size
        val = np.zeros(nlat)
        for i in np.arange(nlat):
            if counts[i] >=5:
                tmp_arr = self.data[lat_all==lat[i]]
                val[i] = np.nanmean(tmp_arr)
            else:
                val[i] = None
        return lat, val

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
        # preprocess variable name, change sign if variable name start with '-'
        if var[0] == '-':
            var = var[1:]
            change_sign = True
        else:
            change_sign = False
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
        # change sign if variable name start with '-'
        if change_sign:
            dat = -dat
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
        # preprocess variable name, change sign if variable name start with '-'
        if var[0] == '-':
            var = var[1:]
            change_sign = True
        else:
            change_sign = False
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
        # change sign if variable name start with '-'
        if change_sign:
            dat = -dat
        out = GOTMTimeseries(time=time, data=dat, name=var)
        # close data
        self.close()
        return out

    def diag_forcing_regime_BG12(self, tidx_start=None, tidx_end=None, cfrac=0.25):
        """Diagnose the forcing regime according to the dissipation based
           definition in Belcher et al., 2012

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :cfrac: (float, optional) critical fraction, 0<cfrac<1/3
        :returns: (int) forcing regime flag

        """
        # read data
        ts_laturb = self.read_timeseries('La_Turb', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_ustar = self.read_timeseries('u_taus', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_hbl = self.read_timeseries('bld_nuh', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_obj = self.read_timeseries('bflux', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True)
        ts_bflux = ts_obj.data[1:]
        m_bflux = ts_obj.data_mean
        # stable or unstable on average
        if m_bflux > 0:
            unstable = -1
        else:
            unstable = 1
        # forcing regime in unstable condition
        mask1 = ts_bflux < 0
        mask2 = ts_ustar > 1e-10
        mask3 = ts_laturb > 1e-10
        bmask = [all(imask) for imask in zip(mask1, mask2, mask3)]
        if sum(bmask) < ts_bflux.size/3:
            # not enough data
            forcing_regime = 8
        else:
            comp_ST = 2.0*(1.0-np.exp(-0.5*ts_laturb[bmask]))
            comp_LT = 0.22/ts_laturb[bmask]**2
            comp_CT = -0.3*ts_bflux[bmask]*ts_hbl[bmask]/ts_ustar[bmask]**3
            comp_total = comp_ST + comp_LT + comp_CT
            frac_ST = comp_ST/comp_total
            frac_LT = comp_LT/comp_total
            frac_CT = comp_CT/comp_total
            mfrac_ST = np.mean(frac_ST)
            mfrac_LT = np.mean(frac_LT)
            mfrac_CT = np.mean(frac_CT)
            if mfrac_LT < cfrac and mfrac_CT < cfrac:
                # ST dominant
                forcing_regime = 1
            elif mfrac_ST < cfrac and mfrac_CT < cfrac:
                # LT dominant
                forcing_regime = 2
            elif mfrac_ST < cfrac and mfrac_LT < cfrac:
                # CT dominant
                forcing_regime = 3
            elif mfrac_ST >= cfrac and mfrac_LT >= cfrac and mfrac_CT < cfrac:
                # combined ST and LT
                forcing_regime = 4
            elif mfrac_ST >= cfrac and mfrac_CT >= cfrac and mfrac_LT < cfrac:
                # combined ST and CT
                forcing_regime = 5
            elif mfrac_LT >= cfrac and mfrac_CT >= cfrac and mfrac_ST < cfrac:
                # combined LT and CT
                forcing_regime = 6
            else:
                # combined ST, LT and CT
                forcing_regime = 7
        return forcing_regime*unstable

    def diag_forcing_regime_LF17(self, tidx_start=None, tidx_end=None, cfrac=0.25):
        """Diagnose the forcing regime according to the entrainment buoyancy
           flux based definition in Li and Fox-Kemper, 2017

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :cfrac: (float, optional) critical fraction, 0<cfrac<1/3
        :returns: (int) forcing regime flag

        """
        # read data
        ts_sl = self.read_timeseries('La_SL', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_ustar = self.read_timeseries('u_taus', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_hbl = self.read_timeseries('bld_nuh', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True).data[1:]
        ts_obj = self.read_timeseries('bflux', tidx_start=tidx_start, tidx_end=tidx_end,
                ignore_time=True)
        ts_bflux = ts_obj.data[1:]
        m_bflux = ts_obj.data_mean
        # stable or unstable on average
        if m_bflux > 0:
            unstable = -1
        else:
            unstable = 1
        # forcing regime in unstable condition
        mask1 = ts_bflux < 0
        mask2 = ts_ustar > 1e-10
        mask3 = ts_sl > 1e-10
        bmask = [all(imask) for imask in zip(mask1, mask2, mask3)]
        if sum(bmask) < ts_bflux.size/3:
            # not enough data
            forcing_regime = 8
        else:
            comp_ST = 0.17
            comp_LT = 0.083/ts_sl[bmask]**2
            comp_CT = -0.15*ts_bflux[bmask]*ts_hbl[bmask]/ts_ustar[bmask]**3
            comp_total = comp_ST + comp_LT + comp_CT
            frac_ST = comp_ST/comp_total
            frac_LT = comp_LT/comp_total
            frac_CT = comp_CT/comp_total
            mfrac_ST = np.mean(frac_ST)
            mfrac_LT = np.mean(frac_LT)
            mfrac_CT = np.mean(frac_CT)
            if mfrac_LT < cfrac and mfrac_CT < cfrac:
                # ST dominant
                forcing_regime = 1
            elif mfrac_ST < cfrac and mfrac_CT < cfrac:
                # LT dominant
                forcing_regime = 2
            elif mfrac_ST < cfrac and mfrac_LT < cfrac:
                # CT dominant
                forcing_regime = 3
            elif mfrac_ST >= cfrac and mfrac_LT >= cfrac and mfrac_CT < cfrac:
                # combined ST and LT
                forcing_regime = 4
            elif mfrac_ST >= cfrac and mfrac_CT >= cfrac and mfrac_LT < cfrac:
                # combined ST and CT
                forcing_regime = 5
            elif mfrac_LT >= cfrac and mfrac_CT >= cfrac and mfrac_ST < cfrac:
                # combined LT and CT
                forcing_regime = 6
            else:
                # combined ST, LT and CT
                forcing_regime = 7
        return forcing_regime*unstable

    def _get_derived_profile(self, name=None, list_keys=False):
        """Find the derived profile variable

        :name: (str) variable name
        :list_keys: (bool) simply return the list of keys if True
        :returns: (function/str) corresponding function or list of keys

        """
        switcher = {
                'buoy': self._get_buoyancy,
                'spice': self._get_spice,
                'wt': self._get_wt,
                'ws': self._get_ws
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

    def _get_wt(self, tidx_start=None, tidx_end=None):
        """Find the vertical temperature flux (K m/s)

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) surface buoyancy flux

        """
        # temperature profiles
        temp = self.dataset.variables['temp'][tidx_start:tidx_end,:,0,0]
        # turbulent diffusivity of heat
        nuh  = self.dataset.variables['nuh'][tidx_start:tidx_end,:,0,0]
        # nonlocal temperature flux
        gamh = self.dataset.variables['gamh'][tidx_start:tidx_end,:,0,0]
        # layer thickness
        h    = self.dataset.variables['h'][tidx_start:tidx_end,:,0,0]
        # depth
        zi   = self.dataset.variables['zi'][0,1:,0,0]
        # surface temperature flux
        wt0  = -self.dataset.variables['heat'][tidx_start:tidx_end,0,0]/cp/rho_0
        # vertical temperature flux
        wt = np.zeros(temp.shape)
        wt[:,-1] = wt0
        wt[:,0:-1] = -2.0*nuh[:,1:-1]*(temp[:,1:]-temp[:,0:-1])/(h[:,1:]+h[:,0:-1]) + gamh[:,1:-1]
        return wt, zi

    def _get_ws(self, tidx_start=None, tidx_end=None):
        """Find the vertical salinity flux (g/kg m/s)

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) surface buoyancy flux

        """
        # salinity profiles
        salt = self.dataset.variables['salt'][tidx_start:tidx_end,:,0,0]
        # surface salinity
        salt0 = salt[:,-1]
        # turbulent diffusivity of salinity
        nus  = self.dataset.variables['nus'][tidx_start:tidx_end,:,0,0]
        # nonlocal salinity flux
        gams = self.dataset.variables['gams'][tidx_start:tidx_end,:,0,0]
        # layer thickness
        h    = self.dataset.variables['h'][tidx_start:tidx_end,:,0,0]
        # depth
        zi   = self.dataset.variables['zi'][0,1:,0,0]
        # surface salinity flux
        ws0  =  (self.dataset.variables['precip'][tidx_start:tidx_end,0,0]
               - self.dataset.variables['evap'][tidx_start:tidx_end,0,0])*salt0
        # vertical temperature flux
        ws = np.zeros(salt.shape)
        ws[:,-1] = ws0
        ws[:,0:-1] = -2.0*nus[:,1:-1]*(salt[:,1:]-salt[:,0:-1])/(h[:,1:]+h[:,0:-1]) + gams[:,1:-1]
        return ws, zi

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
                'bld_nuh': self._get_bld_nuh,
                'bld_tke': self._get_bld_tke,
                'bld_kpp': self._get_bld_kpp,
                'bld_epbl': self._get_bld_epbl,
                'Nsqr_mld': self._get_Nsqr_mld,
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
        hbl = self._get_derived_timeseries('bld_nuh')(tidx_start=tidx_start, tidx_end=tidx_end)
        # surface layer: upper 20% of the boundary layer
        hsl = -0.2*hbl
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

    def _get_bflux(self, tidx_start=None, tidx_end=None, radiative_heating=False):
        """Find the surface buoyancy flux

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) surface buoyancy flux

        """
        # surface temperature and salinity
        temp0 = self.dataset.variables['temp'][tidx_start:tidx_end,-1,0,0]
        salt0 = self.dataset.variables['salt'][tidx_start:tidx_end,-1,0,0]
        # surface heat flux
        hflux = self.dataset.variables['heat'][tidx_start:tidx_end,0,0]
        # correction for penetrative solar radiation if requested,
        # otherwise use the total solar radiation
        if radiative_heating:
            rad   = self.dataset.variables['rad'][tidx_start:tidx_end,:,0,0]
            # get boundary layer depth
            hbl   = self._get_derived_timeseries('bld_nuh')(tidx_start=tidx_start, tidx_end=tidx_end)
            z     = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
            nt    = temp0.shape[0]
            rflux = np.zeros(nt)
            for i in np.arange(nt):
                ihbl = np.argmin(np.abs(z[i,:]+hbl[i]))
                rflux[i] = rad[i,-1]-rad[i,ihbl]
        else:
            # total solar radiation
            rflux = self.dataset.variables['I_0'][tidx_start:tidx_end,0,0]
        # surface temperature flux
        tflux = (hflux+rflux)/cp/rho_0
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
        hbl   = self._get_derived_timeseries('bld_nuh')(tidx_start=tidx_start, tidx_end=tidx_end)
        # get surface buoyancy flux
        bflux = self._get_bflux(tidx_start=tidx_start, tidx_end=tidx_end)
        # friction velocity
        ustar = self.dataset.variables['u_taus'][tidx_start:tidx_end,0,0]
        # Monin-Obukhov length
        Lmo = ustar**3.0/kappa/bflux
        # filter out zeros
        Lmo = np.ma.array(Lmo, mask=(Lmo==0))
        # h over L
        hoL = -hbl/Lmo
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
        nz = Nsqr.shape[1]
        mld = np.zeros(nt)
        # find the indices where N^2 reaches its maximum
        # add small noise that increase with depth to find the shallowest
        # occurrence when N^2 is constant
        noise = np.array(np.arange(nz))*1e-15
        idx_max = np.argmax(Nsqr+noise, 1)
        for i in np.arange(nt):
            mld[i] = z[i,idx_max[i]]
        return np.abs(mld)

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
        return np.abs(mld)

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
        return np.abs(mld)

    def _get_bld_nuh(self, tidx_start=None, tidx_end=None, nuh_bg=1e-5):
        """Find the boundary layer depth defined as the depth where
           the turbulent diffusivity first drops to a background value

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :nuh_bg: (float, optional) background diffusivity
        :returns: (numpy array) mixed layer depth

        """
        nuh = self.dataset.variables['nuh'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['zi'][tidx_start:tidx_end,:,0,0]
        nt = nuh.shape[0]
        nz = nuh.shape[1]
        bld = np.zeros(nt)
        for i in np.arange(nt):
            idxlist = np.where(nuh[i,:]<nuh_bg)[0]
            idxlist = idxlist[idxlist<nz-1]
            if idxlist.size==0:
                bld[i] = z[i,0]
            elif np.max(idxlist)<nz-2:
                idx1 = np.max(idxlist)
                idx0 = idx1+1
                bld[i] = z[i,idx0] - (z[i,idx0]-z[i,idx1]) * \
                         (nuh[i,idx0]-nuh_bg) / (nuh[i,idx0]-nuh[i,idx1])
            else:
                bld[i] = z[i,-1]
        return np.abs(bld)

    def _get_bld_tke(self, tidx_start=None, tidx_end=None, tke_crit=1e-7):
        """Find the boundary layer depth defined as the depth where
           the turbulent kinetic energy approaches zero (equals a small
           critical value.

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :tke_crit: (float, optional) critial TKE in m^2/s^2
        :returns: (numpy array) mixed layer depth

        """
        tke = self.dataset.variables['tke'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['zi'][tidx_start:tidx_end,:,0,0]
        nt = tke.shape[0]
        nz = tke.shape[1]
        bld = np.zeros(nt)
        for i in np.arange(nt):
            idxlist = np.where(tke[i,:]<tke_crit)[0]
            idxlist = idxlist[idxlist<nz-1]
            if idxlist.size==0:
                bld[i] = z[i,0]
            elif np.max(idxlist)<nz-2:
                idx1 = np.max(idxlist)
                idx0 = idx1+1
                bld[i] = z[i,idx0] - (z[i,idx0]-z[i,idx1]) * \
                         (tke[i,idx0]-tke_crit) / (tke[i,idx0]-tke[i,idx1])
            else:
                bld[i] = z[i,-1]
        return np.abs(bld)

    def _get_bld_kpp(self, tidx_start=None, tidx_end=None):
        """Return the boundary layer depth in KPP in the output

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) mixed layer depth

        """
        bld = self.dataset.variables['KPP_OSBL'][tidx_start:tidx_end,0,0]
        return np.abs(bld)

    def _get_bld_epbl(self, tidx_start=None, tidx_end=None):
        """Return the boundary layer depth in EPBL in the output

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) mixed layer depth

        """
        bld = self.dataset.variables['ePBL_OSBL'][tidx_start:tidx_end,0,0]
        return np.abs(bld)

    def _get_Nsqr_mld(self, tidx_start=None, tidx_end=None):
        """Find the stratification N^2 at the base of the
           mixed layer (density criterion)

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (numpy array) N^2

        """
        # get mixed layer depth
        hml = self._get_derived_timeseries('mld_deltaR')(tidx_start=tidx_start, tidx_end=tidx_end)
        # read N^2
        Nsqr = self.dataset.variables['NN'][tidx_start:tidx_end,:,0,0]
        z = self.dataset.variables['z'][tidx_start:tidx_end,:,0,0]
        nt = Nsqr.shape[0]
        Nsqr_mld = np.zeros(nt)
        for i in np.arange(nt):
            ihml = np.argmin(np.abs(z[i,:]+hml[i]))
            Nsqr_mld[i] = Nsqr[i, ihml]
        return Nsqr_mld

    def _get_pez(self, tidx_start=None, tidx_end=None, depth=None):
        """Calculate the total potential energy (PE) integrated from
           surface to a certain depth

        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :depth: (float, optional) reference depth
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
        if ref_cname is None:
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
            dataset.close()
        return lat, lon

    def diagnostics(self, name=None, **kwargs):
        """Find the diagnostics

        :name: (str) name of diagnostics
        :returns: (GOTMMap object) requested diagnostics in GOTMMap object
                                          format or a list of names of the supported
                                          diagnostics

        """
        if name == "mld_deltaR_mean":
            out = self.mean_state_timeseries(var='mld_deltaR', **kwargs)
        elif name == "mld_deltaRp1_mean":
            out = self.mean_state_timeseries(var='mld_deltaR', deltaR=0.1, **kwargs)
        elif name == "PE_delta":
            out = self.delta_timeseries(var='PE', **kwargs)
        elif name == "SST_mean":
            out = self.mean_state_timeseries(var='sst', **kwargs)
        elif name == 'SSS_mean':
            out = self.mean_state_timeseries(var='sss', **kwargs)
        elif name == 'Nsqr_mld_mean':
            out = self.mean_state_timeseries(var='Nsqr_mld', **kwargs)
        elif name == 'forcing_regime_BG12':
            out = self.forcing_regime_map(forcing_reg_type='BG12', **kwargs)
        elif name == 'forcing_regime_LF17':
            out = self.forcing_regime_map(forcing_reg_type='LF17', **kwargs)
        else:
            raise ValueError('Diagnostics \'{}\' not found.'.format(name))
        return out

    def mean_state_profile(self, var, tidx_start=None, tidx_end=None, zidx_start=None, zidx_end=None, **kwargs):
        """Return the mean state of a profile variable

        :var: (str) variable name
        :tidx_start: (int, optional) starting index for time
        :tidx_end: (int, optional) ending index for time
        :zidx_start: (int, optional) starting index for z
        :zidx_end: (int, optional) ending index z
        :returns: (GOTMMap object) mean state

        """
        tmp = GOTMOutputData(self._paths[0], init_time_location=False)
        prfl_dat = tmp.read_profile(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data[:,zidx_start:zidx_end]
        nt = prfl_dat.shape[0]
        nz = prfl_dat.shape[1]
        dat = np.zeros([self.ncase, nt, nz])
        dat[0,:,:] = prfl_dat
        npcount = np.floor(self.ncase/20)
        for i in np.arange(self.ncase-1)+1:
            if np.mod(i, npcount) == 0:
                print('{:6.1f} %'.format(i/self.ncase*100.0))
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            dat[i,:,:] = tmp.read_profile(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data[:,zidx_start:zidx_end]
        mdat = np.mean(dat, axis=(1,2))
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        print('{:6.1f} %'.format(100.0))
        return out

    def mean_state_timeseries(self, var, fillvalue=None, tidx_start=None, tidx_end=None, **kwargs):
        """Return the mean state of a timeseries variable

        :var: (str) variable name
        :fillvalue: (float) invalid value
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (GOTMMap object) mean state

        """
        tmp = GOTMOutputData(self._paths[0], init_time_location=False)
        ts_dat = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data
        nt = ts_dat.shape[0]
        dat = np.zeros([self.ncase, nt])
        dat[0,:] = ts_dat
        npcount = np.floor(self.ncase/20)
        for i in np.arange(self.ncase-1)+1:
            if np.mod(i, npcount) == 0:
                print('{:6.1f} %'.format(i/self.ncase*100.0))
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            dat[i,:] = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data
        if fillvalue is not None:
            dat[dat==fillvalue] = np.nan
        mdat = dat.mean(axis=1)
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        print('{:6.1f} %'.format(100.0))
        return out

    def delta_timeseries(self, var, fillvalue=None, tidx_start=None, tidx_end=None, **kwargs):
        """Return the net change of a timeseries variable

        :var: (str) variable name
        :tidx_start: (int, optional) starting index
        :tidx_end: (int, optional) ending index
        :returns: (GOTMMap object) mean state

        """
        tmp = GOTMOutputData(self._paths[0], init_time_location=False)
        ts_dat = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data
        nt = ts_dat.shape[0]
        dat = np.zeros([self.ncase, nt])
        dat[0,:] = ts_dat
        npcount = np.floor(self.ncase/20)
        for i in np.arange(self.ncase-1)+1:
            if np.mod(i, npcount) == 0:
                print('{:6.1f} %'.format(i/self.ncase*100.0))
            tmp = GOTMOutputData(self._paths[i], init_time_location=False)
            dat[i,:] = tmp.read_timeseries(var, tidx_start=tidx_start, tidx_end=tidx_end, ignore_time=True, **kwargs).data
        if fillvalue is not None:
            dat[dat==fillvalue] = np.nan
        mdat = dat[:,-1] - dat[:,0]
        # create GOTMMap object
        out = GOTMMap(data=mdat, lon=self.lon, lat=self.lat, name=var)
        print('{:6.1f} %'.format(100.0))
        return out

    def forcing_regime_map(self, forcing_reg_type='BG12', fillvalue=None, **kwargs):
        """Return indices for the forcing regime defined in diag_forcing_regime_*.

        :forcing_reg_type: (str) type of forcing regime: BG12 or LF17
        :fillvalue: (float) invalid value, not for now
        :returns: (GOTMMap object) forcing regime

        """
        forcing_regime = np.zeros(self.ncase)
        npcount = np.floor(self.ncase/20)
        if forcing_reg_type == 'BG12':
            for i in np.arange(self.ncase):
                if np.mod(i, npcount) == 0:
                    print('{:6.1f} %'.format(i/self.ncase*100.0))
                tmp = GOTMOutputData(self._paths[i], init_time_location=False)
                forcing_regime[i] = tmp.diag_forcing_regime_BG12(**kwargs)
        elif forcing_reg_type == 'LF17':
            for i in np.arange(self.ncase):
                if np.mod(i, npcount) == 0:
                    print('{:6.1f} %'.format(i/self.ncase*100.0))
                tmp = GOTMOutputData(self._paths[i], init_time_location=False)
                forcing_regime[i] = tmp.diag_forcing_regime_LF17(**kwargs)
        else:
            raise ValueError('forcing regime name \'{}\' not supported.')
        # create GOTMMap object
        out = GOTMMap(data=forcing_regime, lon=self.lon, lat=self.lat, name='forcing_regime_'+forcing_reg_type)
        print('{:6.1f} %'.format(100.0))
        return out

#--------------------------------
# Functions
#--------------------------------

def plot_dist_3p(hst, xi, yi, axis=None, filled=False, fcolors=None, **kwargs):
    """Plot bi-dimensional histogram. Show the contours of the
       histogram which enclose the highest 30%, 60%, and 90%
       centered distribution.

    :his: (2D numpy array) bi-dimensional histogram
    :xi: (1D numpy array) centers of x dimension
    :yi: (1D numpy array) centers of y dimension
    :axis: (matplotlib.axes, optional) axis to plot figure on
    :filled: (bool) filled contour if True
    :fcolors: (list, optional) color string or sequence of colors, optional)
    :return: (matplotlib figure object) figure

    """
    vl = [0.3, 0.6, 0.9]
    fig = plot_dist_xp(hst, xi, yi, axis=axis, levels=vl, filled=filled, fcolors=fcolors, **kwargs)
    return fig

def plot_dist_4p(hst, xi, yi, axis=None, filled=False, fcolors=None, **kwargs):
    """Plot bi-dimensional histogram. Show the contours of the
       histogram which enclose the highest 30%, 60%, 90% and 99%
       centered distribution.

    :his: (2D numpy array) bi-dimensional histogram
    :xi: (1D numpy array) centers of x dimension
    :yi: (1D numpy array) centers of y dimension
    :axis: (matplotlib.axes, optional) axis to plot figure on
    :filled: (bool) filled contour if True
    :fcolors: (list, optional) color string or sequence of colors, optional)
    :return: (matplotlib figure object) figure

    """
    vl = [0.3, 0.6, 0.9, 0.99]
    fig = plot_dist_xp(hst, xi, yi, axis=axis, levels=vl, filled=filled, fcolors=fcolors, **kwargs)
    return fig

def plot_dist_xp(hst, xi, yi, axis=None, levels=None, filled=False, fcolors=None, **kwargs):
    """Plot bi-dimensional histogram. Show the contours of the
       histogram which enclose the highest p1%, p2%, ... and pN%
       centered distribution.

    :his: (2D numpy array) bi-dimensional histogram
    :xi: (1D numpy array) centers of x dimension
    :yi: (1D numpy array) centers of y dimension
    :axis: (matplotlib.axes, optional) axis to plot figure on
    :levels: (list of float, optional) contour levels, 0.0-1.0
    :filled: (bool) filled contour if True
    :fcolors: (list, optional) color string or sequence of colors
    :return: (matplotlib figure object) figure

    """
    # use curret axis if not specified
    if axis is None:
        axis = plt.gca()
    hsum = np.sum(hst)
    hlist = -np.sort(-hst.flatten())/hsum
    hcum = np.cumsum(hlist)
    vl = levels
    nv = len(vl)
    vlev = np.zeros(nv)
    for i in np.arange(nv):
        ind = np.argmin(abs(hcum-vl[i]))
        vlev[i] = hlist[ind]
    pdfData = hst/hsum
    pdfData[pdfData==0] = 1e-12
    if not filled:
        fig = axis.contour(xi, yi, np.log10(np.transpose(pdfData)), levels=np.log10(vlev[::-1]), **kwargs)
    else:
        if fcolors is None:
            cmap = cm.get_cmap('bone')
            fcolors = cmap(np.linspace(1.0, 0.0, 11)[0:nv+1])
        else:
            nfc = len(fcolors)
            if nfc != nv+1:
                raise ValueError('Length of fcolors should equal to number of levels + 1.')
        fig = axis.contourf(xi, yi, np.log10(np.transpose(pdfData)), levels=np.log10(vlev[::-1]),
                            colors=fcolors, extend='both', **kwargs)
    return fig


def plot_regime_diagram_background_BG12(axis=None):
    """Plot the background of the regime diagram
       following Fig. 3 of Belcher et al., 2012

    :axis: (matplotlib.axes, optional) axis to plot figure on

    """
    if axis is None:
        axis = plt.gca()

    # range of power
    xpr = [-1, 1]
    ypr = [-3, 3]
    # range
    xlims = [10**i for i in xpr]
    ylims = [10**i for i in ypr]
    # size of x and y
    nx = 500
    ny = 500
    xx = np.logspace(xpr[0], xpr[1], nx)
    yy = np.logspace(ypr[0], ypr[1], ny)
    zz1 = np.zeros([nx, ny])
    zz2 = np.zeros([nx, ny])
    zz3 = np.zeros([nx, ny])
    for i in np.arange(nx):
        for j in np.arange(ny):
            zz1[i,j] = 2*(1-np.exp(-0.5*xx[i]))
            zz2[i,j] = 0.22*xx[i]**(-2)
            zz3[i,j] = 0.3*xx[i]**(-2)*yy[j]
    zz = zz1 + zz2 + zz3
    axis.contourf(xx, yy, np.transpose(np.log10(zz)),
                  levels=[-0.1, 0, 0.1, 0.25, 0.5, 1, 2, 3, 4],
                  cmap='summer', extend='both')
    axis.contour(xx, yy, np.transpose(np.log10(zz)),
                  levels=[-0.1, 0, 0.1, 0.25, 0.5, 1, 2, 3, 4],
                  colors='darkgray')
    axis.contour(xx, yy, np.transpose(zz1/zz), levels=0.9, colors='k',
                linestyles='-', linewidths=2)
    axis.contour(xx, yy, np.transpose(zz2/zz), levels=0.9, colors='k',
                linestyles='-', linewidths=2)
    axis.contour(xx, yy, np.transpose(zz3/zz), levels=0.9, colors='k',
                linestyles='-', linewidths=2)
    axis.set_xlim(xlims)
    axis.set_ylim(ylims)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel('La$_t$')
    axis.set_ylabel('$h/L_L$')
    axis.set_aspect(aspect=1/3)
    axis.text(0.85, 3e-3, '0', color='k', fontsize=8, rotation=-90)
    axis.text(1.6, 1e-2, '0.1', color='k', fontsize=8, rotation=-90)
    axis.text(3.8, 1e-1, '0.25', color='k', fontsize=8, rotation=-90)
    axis.text(4, 1e2, '0.5', color='k', fontsize=8, rotation=33)
    axis.text(3.2, 3e2, '1', color='k', fontsize=8, rotation=36)
    axis.text(0.53, 1e2, '2', color='k', fontsize=8, rotation=38)
    axis.text(0.3, 3.1e2, '3', color='k', fontsize=8, rotation=39)
    axis.text(0.12, 5e2, '4', color='k', fontsize=8, rotation=40)
    axis.text(0.11, 4e-3, 'Langmuir', bbox=dict(boxstyle="square",ec='k',fc='w'))
    axis.text(3, 4e-3, 'Shear', bbox=dict(boxstyle="square",ec='k',fc='w'))
    axis.text(0.13, 1e2, 'Convection', bbox=dict(boxstyle="square",ec='k',fc='w'))


def plot_regime_diagram_background_L19(axis=None):
    """Plot the background of the reegime diagram
       in Li et al., 2019

    :axis: (matplotlib.axes, optional) axis to plot figure on

    """
    if axis is None:
        axis = plt.gca()
    # range of power
    xpr = [-1, 1]
    ypr = [-3, 3]
    # range
    xlims = [10**i for i in xpr]
    ylims = [10**i for i in ypr]
    # background following Fig. 3 of Belcher et al., 2012
    nx = 500
    ny = 500
    xx = np.logspace(xpr[0], xpr[1], nx)
    yy = np.logspace(ypr[0], ypr[1], ny)
    zz1 = np.zeros([nx, ny])
    zz2 = np.zeros([nx, ny])
    zz3 = np.zeros([nx, ny])
    for i in np.arange(nx):
        for j in np.arange(ny):
            zz1[i,j] = 2*(1-np.exp(-0.5*xx[i]))
            zz2[i,j] = 0.22*xx[i]**(-2)
            zz3[i,j] = 0.3*xx[i]**(-2)*yy[j]
    zz = zz1 + zz2 + zz3

    rz_ST = zz1/zz
    rz_LT = zz2/zz
    rz_CT = zz3/zz
    fr = np.ones(zz.shape) * 7
    cfrac = 0.25
    fr[(rz_LT<cfrac) & (rz_CT<cfrac)] = 1
    fr[(rz_ST<cfrac) & (rz_CT<cfrac)] = 2
    fr[(rz_ST<cfrac) & (rz_LT<cfrac)] = 3
    fr[(rz_ST>=cfrac) & (rz_LT>=cfrac) & (rz_CT<cfrac)] = 4
    fr[(rz_ST>=cfrac) & (rz_CT>=cfrac) & (rz_LT<cfrac)] = 5
    fr[(rz_LT>=cfrac) & (rz_CT>=cfrac) & (rz_ST<cfrac)] = 6

    color_list = ['firebrick','forestgreen','royalblue','gold','orchid','turquoise','w']
    cb_ticks = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
    cmap, norm = from_levels_and_colors(cb_ticks, color_list)
    axis.contourf(xx, yy, np.transpose(fr), cmap=cmap, norm=norm)
    axis.contour(xx, yy, np.transpose(fr), colors='darkgray')
    axis.set_xlim(xlims)
    axis.set_ylim(ylims)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel('La$_t$')
    axis.set_ylabel('$h/L_L$')
    axis.set_aspect(aspect=1/3)
    axis.text(0.11, 4e-3, 'Langmuir', bbox=dict(boxstyle="square",ec='k',fc='w'))
    axis.text(3, 4e-3, 'Shear', bbox=dict(boxstyle="square",ec='k',fc='w'))
    axis.text(0.13, 1e2, 'Convection', bbox=dict(boxstyle="square",ec='k',fc='w'))


def plot_forcing_regime_map(f_regime, axis=None, add_colorbar=True, **kwargs):
    """Plot forcing regime in scatters on a map

    :f_regime: (GOTMMap object) forcing regime
    :axis: (matplotlib.axes, optional) axis to plot figure on
    :add_colorbar: (bool) do not add colorbar if False
    :**kwargs: (keyword arguments) to be passed to mpl_toolkits.basemap.scatter()
    :return: (matplotlib figure object) figure

    """
    # use curret axis if not specified
    if axis is None:
        axis = plt.gca()
    # plot map
    m = Basemap(projection='cyl', llcrnrlat=-72, urcrnrlat=72, llcrnrlon=20, urcrnrlon=380, ax=axis)
    # plot coastlines, draw label meridians and parallels.
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='gray',lake_color='lightgray')
    m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1])
    m.drawmeridians(np.arange(-180.,181.,60.), labels=[1,0,0,1])
    # plot forcing regime
    data = np.abs(f_regime.data)
    lat = f_regime.lat
    lon = f_regime.lon
    # shift longitude
    lon = np.where(lon < 20., lon+360., lon)
    x, y = m(lon, lat)
    # levels
    levels = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
    cb_ticks = [1, 2, 3, 4, 5, 6, 7, 8]
    cb_ticks_labels = ['S', 'L', 'C', 'SL', 'SC', 'LC', 'SLC', 'NA']
    color_list = ['firebrick','forestgreen','royalblue','gold','orchid','turquoise','w','gray']
    cmap = colors.LinearSegmentedColormap.from_list('rgb', color_list, len(color_list))
    bounds = np.array(levels)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=len(color_list))
    fig = m.scatter(x, y, marker='.', s=32, c=data, norm=norm, cmap=cmap, **kwargs)
    # mark stable regions
    freg_data = np.where(np.isnan(f_regime.data), 999, f_regime.data)
    smask = freg_data<0
    lat_s = lat[smask]
    lon_s = lon[smask]
    # shift longitude
    lon_s = np.where(lon_s < 20., lon_s+360., lon_s)
    x_s, y_s = m(lon_s, lat_s)
    # fig_s = m.scatter(x_s, y_s, marker='*', s=6, c='black', linewidth=0.1, **kwargs)
    fig_s = m.scatter(x_s, y_s, marker='x', s=7, c='black', linewidth=0.4, **kwargs)
    # add colorbar
    if add_colorbar:
        cb = m.colorbar(fig, ax=axis, ticks=cb_ticks)
        cb.ax.set_yticklabels(cb_ticks_labels)
    return fig
