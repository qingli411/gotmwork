# GOTM work

A work directory of [GOTM](http://gotm.net), which contains a set of scripts and tools to preprocess the input data, set up runs, and analyze and visualize the output data.

---
## Quick Start

### Installation

Create a local copy of `gotmwork` by

```
git clone https://github.com/qingli411/gotmwork
```

### Setup

The script `setup_gotmwork.sh` provides a step-by-step guide to set up the work environment for GOTM. These steps include:

1. Set up environment variables.
2. Download source code of [CVMix](http://cvmix.github.io).
3. Download source code of [GOTM](https://github.com/qingli411/gotm).
4. Compile CVMix.
5. Compile GOTM.
6. Set up Python environment for GOTM via Conda.

To use it, chanage directory to `gotmwork`, type in
```
./setup_gotmwork.sh
```
and follow the instructions.

Here is a brief description of each step:

- Step 1 generates a file `.gotmwork_env.sh` in the HOME directory which saves all the necessary environment variables. You may optionally add `source $HOME/.gotmwork_env.sh` in the file `$HOME/.bashrc` (bash shell) to automatically set up these environment variables when opening a new terminal.
- Step 2 downloads the source code of CVMix from Github and checkout the required tag. See `scripts/get_cvmix.sh`.
- Step 3 downloads the source code of GOTM from Github and checkout the required tag. See `scripts/get_gotm.sh`.
- Step 4 compiles CVMix library using `build_cvmix.sh`. You will be asked to input a Fortran compiler (e.g., gfortran) and the path of netCDF when running it for the first time. See [CVMix Homepage](http://cvmix.github.io) for more information on how to build CVMix.
- Step 5 compiles GOTM using `build_gotm.sh`. CVMix has to be built before building GOTM.
- Step 6 sets up Python environment for [preprocess data ](#Preprocess-data) and [analysis & visualization](#Analysis-and-Visualization) using [Conda](https://docs.conda.io/en/latest/). It will create a new conda environment `gotm` from `gotm_env_*.yml` and activate it using `set_conda_env.sh`.


### Compile CVMix and GOTM

In cases of compiling steps are skipped when [setting up](#Setup) the gotmwork environment or source code of CVMix and GOTM are changes, CVMix and GOTM can be compiled by
```
./build_cvmix.sh [-clean] [-build]
```
and
```
./build_gotm.sh [-clean] [-build]
```
Use optional arguments `-build` to build, `-clean` to clean the old build, and `-clean -build` to do a clean build.

### Preprocess data
`tools/` and `scripts/` contain some [tools](#A-List-of-Tools) and [scripts](#A-List-of-Scripts) to preprocess data for GOTM.
`tools/gotmtool.py` provides functions (e.g., `write_ts()`, `write_pfl()` and `write_spec()`) to write data to file in the format that GOTM5 requires.

### Run a case
Change directory to a test case (see [Test Cases](#Test-Cases) for more detail) and run
```
./case_run
```

### Analysis and Visualization
`tools/gotmanalysis.py` provides classes and functions for data analysis and visualization. Some examples using these classes and functions are in `visualization/examples`.


---
## Test Cases

Test cases. In each case, `case_test` sets up the namelist, preprocess the input data and run the simulation, whereas `case_run` (if exists) runs the simulation from preprocessed input data. `case_postproc.sh` contains steps to postprocess output data and is used by `case_run` and `case_test`, either for visulaization or manipulation of the data.

- __COREII__
Run one set of simulations in each 4 by 4 degree box to cover the global ocean, foced by CORE-II.
  - `case_test_multi` sets up multiple runs under CORE-II forcing.
  - `case_run_multi` is similar to `case_test_multi`, but uses preprocessed CORE-II data and is therefore significantly faster.
  - `do_parallel` automatically submit parallel jobs to multiple cores.
  - `kill_all` kills all the jobs.
  - `preproc_data` preprocesses the CORE-II data.

- __JRA55do__
Run one set of simulations in each 4 by 4 degree box to cover the global ocean, foced by JRA55-do.
  - `case_run_multi` sets up multiple runs under JRA55-do forcing using preprocessed data.
  - `do_parallel` automatically submit parallel jobs to multiple cores.
  - `kill_all` kills all the jobs.
  - `preproc_data` preprocesses the JRA55-do data.

- __OCSKEO__

- __OCSPapa__
Single site simulation forced by [Ocean Station Papa](https://www.pmel.noaa.gov/ocs/Papa) data.
  - `case_postproc.sh` is a Bash script to postprocess a single run, used by `case_run`.
- __OSMOSIS__
Single site simulation forced by [OSMOSIS](https://www.bodc.ac.uk/projects/data_management/uk/osmosis/) data.
- __TEST_RES__
Sensitivity test of different boundary layer schemes to different vertical resolutions and time steps.
   - `case_loop.sh` loops over different turbulent methods, vertical resolutions and time steps.
   - `OCSPapa` runs test case using OCS Papa data
   - `OSMOSIS` runs test case using OSMOSIS data
   - `COREII` runs test case using selected COREII data
- __Idealized_Tests__
Set up cases with constant surface forcing and idealzied initial conditions.
- __Idealized_Hurricane__
Set up the idealized hurricane cases of [Reichl et al., 2016](https://doi.org/10.1175/MWR-D-16-0074.1).
- __Idealized_Tests_LF17__
Set up idealized cases using the initial conditions and surface forcing conditions of Case S-L1 and Case S-B in [Li and Fox-Kemper, 2017](https://doi.org/10.1175/JPO-D-17-0085.1) (see their Table~1).
- __Idealized_Tests_MSM97__
Set up idealized cases using the initial conditions and surface forcing conditions of [McWilliams et al., 1997](https://doi.org/10.1017/S0022112096004375).


### Preprocessed Data

The preprocessed input data and namelists for [Test Cases](#Test-Cases) are in the directrory `data/`:

- OCSPapa_20120101-20131204
- OSMOSIS_winter
- OSMOSIS_spring
- COREII_LAT-54_LON254_20080601-20091231
- COREII_LAT10_LON86_20080601-20091231
- COREII_LAT2_LON234_20080601-20091231
- Idealized
- Idealized_Hurricane
- Idealized_Tests_LF17
- Idealized_Tests_MSM97

In each directory the tool `update_nml` can be used to update the namelist from `data/namelist/` in the case where new entries are added.

Also included in this directory are the data description files in XML format, which are used by `case_preproc` to preprocess the input data.

### Namelist

All namelists are in the directory `data/namelist/`. Use `init_namelist` to generate namelist from schemas in the source code according to the type of turbulence closure. It requires Python2 environment and the tool `editscenario`, which can be installed using the script `scripts/install_python_tools.sh`.

---
## A List of Tools

A list of tools in the directory `tools/`.
Most of the tools listed below are written in Python3, some in Bash script. The file `gotmtool.py` contains some shared Python3 functions used by many of the tools. Option `-h` can be used with all tools to get the usage.

| Tool name                  | Description |
| -------------------------- |:----------- |
| `argo_mld`                 | Read Argo temperature and salinity profiles in GOTM input data format and return mixed layer depth based on density threshold|
| `case_preproc`             | Preprocess the input data for GOTM and modify the namelist according to the input xml file. |
| `gotm_archive_data`        | Compress and archive GOTM output data. |
| `gotm_extract_data`        | Extract data from archive generated by `gotm_archive_data`. |
| `gotm_map_quality_control` | Remove runs with NaNs in the output data. |
| `is_sea`                   | Check if the point given by latitude and longitude is a sea point. |
| `nc2dat`                   | Convert observational data (OCS etc.) in netCDF format to formatted text file for GOTM input. |
| `nc2dat_argo`              | Convert Argo profile data in netCDF format to formatted text file for GOTM input. |
| `nc2dat_cdip_spec`         | Convert CDIP wave spectrum data in netCDF format to formatted text file for GOTM input. |
| `nc2dat_core2swr`          | Read the daily maximum shortwave radiation from COREII in netCDF format, add an idealized diurnal cycle, and output the hourly shortwave radiation into formatted text file for GOTM input. |
|  `nc2dat_latlon`           | Select CORE-II/JRA55-do/CESM data (in netCDF format) at given latitude and longitude and output into formatted text file for GOTM input. |
| `nc2dat_ww3`               | Convert WW3 wave variables and partitioned surface Stokes drift data in netCDF format to formatted text file for GOTM input. |
| `nmlchange`                | Change the entry value of a namelist. |
| `nmlquery`                 | Query the value of en entry in a namelist. |
| `plotpfl`                  | Plot Hovmoller diagram (time-depth) from GOTM output.
| `plotts`                   | Plot time series from GOTM output. Accept multiple variables. |

### Other Tools

- `gotmanalysis.py` contains classes and functions for data analysis and visualization.
- `windwave.py` contains tools to estimate and test Stokes drift computed from empirical wind wave spectrum.

---
## A List of Scripts

A list of scripts in the directory `scripts/`.
Bash scripts to setup the tools and runs, and Matlab and NCL scripts to preprocess the input data.

| Script name                | Description |
| -------------------------- |:----------- |
| `argo_mat2nc.m`            | Matlab script to convert Argo profile data from MAT to netCDF. |
| `case_turbmethod.sh`       | Bash script to set the namelist according to the turbulent methods. |
| `cesm_prep_fluxes.ncl`     | NCL script to prepare the surface fluxes data from CESM output. |
| `cesm_prep_profiles.ncl`   | NCL script to prepare the temperature and salinity profiles from CESM output. |
| `cesm_ww3a_to_gx1v6.ncl`   | NCL script to interpolate the WW3 output data onto POP grid gx1v6. |
| `core2_prep_meteo.ncl`     | NCL script to prepare meteorology data from CORE-II. |
| `get_cvmix.sh`             | Bash script to download CVMix source code from Github. |
| `get_gotm.sh`              | Bash script to download GOTM source code from Github. |
| `install_python_tools.sh`  | Bash script to download and install Python tools for GOTM from Github, including `editscenario`, `xmlstore`, `xmlplot` and `gotmgui`. |
| `jar55do_prep_meteo.ncl`   | NCL script to prepare meteorology data from JRA55-do. |
| `ocs_heatflux.ncl`         | NCL script to prepare the net heat flux (excluding shortwave) data from longwave, sensible and latent heat fluxes for GOTM.
| `roms_dz.m`                | Matlab script to generate ROMS style stretching vertical grid. |
| `set_gotmwork_env.sh`      | Bash script to set up environment variables for GOTM. |

### Other Scripts

Scripts in the root directory:

- `build_cvmix.sh` is a Bash script to build GOTM.
- `build_gotm.sh` is a Bash script to build GOTM.
- `set_tools.sh` is a Bash script to set up paths and tools, used by `case_run`.
- `setup_gotmwork.sh` is a Bash script to set up work environment for GOTM.



