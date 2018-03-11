# gotmwork #

This is my [GOTM](http://gotm.net) working directory.

Use `setup.sh` to set up the source code, test cases and tools from Github. Use `build_src.sh` to compile the source code.


## tools/ ##

A set of tools to preprocess observational data for input and postprocess GOTM output data.

- Change the entry value of a namelist.

  `nmlchange`

- Query the value of en entry in a namelist.

  `nmlquery`

- Convert input observational data from netCDF to formatted text file.

  `nc2dat`

  `nc2dat_latlon` Convert CESM/CORE-II data to text file for given lat and lon

  `nc2dat_cdip_spec` Convert wave spectrum data of CDIP to text file.

  `nc2dat_core2ww3_usp` Convert partitioned Stokes drift of WW3 to text file.

- Check is the point given by latitude and longitude is a sea point.

  `is_sea`

- Plot time series from GOTM output. Accept multiple variables.

  `plotts`

- Plot time series of profile from GOTM.

  `plotpfl`

- Preprocess the input data for GOTM and modify the namelist accordingly.

  `case_preproc`

- Show timeseries in the observation.

  `obs_show`

- Shared functions
  `gotmtool.py`

- NCL script to prepare the net heat flux (excluding shortwave) data from longwave, sensible and latent heat fluxes for GOTM.

  `ocs_heatflux.ncl`

- NCL script to prepare the surface fluxes data and temperature and salinity profiles from CESM output.

  `cesm_prep_fluxes.ncl`

  `cesm_prep_profiles.ncl`

- NCL script to prepare meteorology data from CORE-II.

  `core2_prep_meteo.ncl`

## namelist/ ##

Directory for all namelist


## data/ ##

Directory for data description file in XML


## cases/ ##

Test cases.
