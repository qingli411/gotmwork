# GOTM work

This is my work directory to setup, test and run  [GOTM](http://gotm.net) and is continually growing. A set of scripts and tools are described in the following.

## ./
The script `build_src.sh` is used to compile the source code. Change `${blddir}` (temporary directory for builing the code) and `${srcdir}` (directory of the source code) accordingly before running it.

The script `set_tools.sh` is used by `case_run` (in case/) to set paths and tools.

Note: The script `setup.sh` downloads and sets up the **original** GOTM source code, test cases and tools from Github. It was only used when I initially set up the codebase.   

## ./tools/

This folder contains a set of tools to preprocess observational data for input and postprocess GOTM output data. Most of them are written in Python3 and can be run from the command line. The arguments are managed using `argparse`.

- Shared Python3 functions.

  `gotmtool.py`

- Change the entry value of a namelist.

  `nmlchange`

- Query the value of en entry in a namelist.

  `nmlquery`

- Convert netCDF data to formatted text file.

  `nc2dat` Observational data from OCS etc.

  `nc2dat_argo` Argo profile data

  `nc2dat_cdip_spec` CDIP wave spectrum data.

  `nc2dat_core2ww3` WW3 wave variables and partitioned surface Stokes drift data.

  `nc2dat_latlon` Select CESM/CORE-II data at given lat and lon.

- Check if the point given by latitude and longitude is a sea point.

  `is_sea`

- Plot time series from GOTM output. Accept multiple variables.

  `plotts`

- Plot time series of profile from GOTM.

  `plotpfl`

- Preprocess the input data for GOTM and modify the namelist accordingly.

  `case_preproc`

- Show timeseries in the observation.

  `obs_show`

- NCL script to prepare the net heat flux (excluding shortwave) data from longwave, sensible and latent heat fluxes for GOTM.

  `ocs_heatflux.ncl`

- NCL script to prepare the surface fluxes data and temperature and salinity profiles from CESM output.

  `cesm_prep_fluxes.ncl`

  `cesm_prep_profiles.ncl`

- NCL script to prepare meteorology data from CORE-II.

  `core2_prep_meteo.ncl`

- Matlab script to convert Argo profile data from MAT to netCDF.

  `argo_mat2nc.m`
  
- Script to postprocess a single run, used by `case_run`.

  `case_postproc.sh`
  
- Script to set the namelist according to the turbulent methods.

  `case_turbmethod.sh`

## ./namelist/

Directory for all namelist


## ./data/

This directory contains preprocessed input data for test cases:

- OCSPapa_20120101-20131204
- OSMOSIS_winter
- OSMOSIS_spring
- COREII_LAT-54_LON254_20080115-20091231
- COREII_LAT10_LON86_20080115-20091231
- COREII_LAT2_LON234_20080115-20091231
- Idealized
- Idealized_Hurricane

Also included in this directory are the data description files in XML format, which are used by `case_preproc` to preprocess the input data.


## ./cases/

Test cases. In each case, `case_test` sets up the namelist, preprocess the input data and run the simulation, whereas `case_run` (if exists) runs the simulation from preprocessed input data.

- COREII

  - `case_test_multi` sets up multiple runs under CORE-II forcing, currently one run for each 4 by 4 degree box globally.
  - `case_run_multi` is similar to `case_test_multi`, but uses preprocessed CORE-II data and is therefore significantly faster.
  - `preproc_data` preprocesses the CORE-II data. 
  - `do_parallel` manually distributes jobs to 8 cores on a Mac Pro.
  - `kill_all` kills all the jobs.

- OCSKEO

- OCSPapa

- OSMOSIS

- TEST_RES
Sensitivity test of different boundary layer schemes to different vertical resolutions and time steps.
   - `case_loop.sh` loops over different turbulent methods, vertical resolutions and time steps.
   - `OCSPapa` runs test case using OCS Papa data
   - `OSMOSIS` runs test case using OSMOSIS data

- Idealized_Tests
- Idealized_Hurricane