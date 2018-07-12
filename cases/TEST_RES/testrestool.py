"""
Shared functions

Qing Li, 20180614
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
tooldir_default = os.environ['GOTMWORK_ROOT']+'/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import get_variable

#--------------------------------
# common variables
#--------------------------------
# input directory
dir_in = os.environ['GOTMRUN_ROOT']+'/TEST_RES/'
# output directory
dir_out = os.environ['GOTMFIG_ROOT']+'/TEST_RES/'
# list of cases
case_list = ['OSMOSIS_winter',
             'OSMOSIS_spring',
             'OCSPapa_20130621-20131201',
             'COREII_LAT2_LON234_20080601-20081231',
             'COREII_LAT10_LON86_20080601-20081231',
             'COREII_LAT-54_LON254_20080901-20090831']
# list of depths
depth_list = np.array([-200, -240, -120, -150, -120, -400])
# list of turbulent methods
turbmethod_list = ['KPP-CVMix',
                   'KPPLT-EFACTOR',
                   'KPPLT-ENTR',
                   # 'KPPLT-RWHGK',
                   'OSMOSIS',
                   'EPBL',
                   'EPBL-LT',
                   'K-EPSILON-SG',
                   'SMC',
                   'SMCLT']
# list of legend for turbulent methods
legend_list = ['KPP-CVMix',
               'KPPLT-VR12',
               'KPPLT-LF17',
               # 'KPPLT-RWHGK16',
               'OSMOSIS',
               'ePBL',
               'ePBL-LT',
               'k-epsilon',
               'SMC-KC94',
               'SMCLT-H15']
# list of profiles
pfl_list = ['temp', 'salt', 'rho', 'buoy', 'spice']
# max and min values for colorbar ([ncase, nprofile])
pfl_cmax_list = np.array([[16, 20, 18, 27, 29, 7],
                      [35.9, 35.9, 33.8, 35.1, 34.6, 34.2],
                      [1027.3, 1027.2, 1026.8, 1026.4, 1026.2, 1027.1],
                      [8.e-3, 8.e-3, 3.2e-2, 2.8e-2, 4.2e-2, 2.3e-3],
                      [1.7e-2, 2.2e-2, -0.5e-2, 2.8e-2, 2.4e-2, -1.2e-2]])
pfl_cmin_list = np.array([[12, 12, 4, 12, 14, 4],
                      [35.7, 35.6, 32.2, 34.7, 33.3, 33.9],
                      [1026.0, 1025.3, 1023.2, 1022.7, 1020.9, 1026.7],
                      [-4.e-3, -4.e-3, 0, 4.e-3, 6.e-3, -3.2e-3],
                      [0.9e-2, 0.6e-2, -2.7e-2, 0, 1.2e-2, -1.8e-2]])
pfl_dmax_list = np.array([[1, 1, 1, 1, 1, 1],
                      [0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
                      [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
                      [5.e-4, 1.e-3, 2.e-3, 2.e-3, 1.e-3, 1.e-3],
                      [5.e-4, 1.e-3, 2.e-3, 2.e-3, 1.e-3, 1.e-3]])

def read_pfl(infile, var, tidx_start=None, tidx_end=None):
    """Read profile variable and z (fixed in time) from GOTM file

    :infile: (netCDF4 Dataset) input netCDF file
    :var: (str) variable name
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) profile, z

    """
    varlist = infile.variables.keys()
    if var in varlist:
        dat = infile.variables[var][tidx_start:tidx_end,:,0,0]
        coord = infile.variables[var].coordinates
        if 'zi' in coord:
            z = infile.variables['zi'][0,:,0,0]
        elif 'z' in coord:
            z = infile.variables['z'][0,:,0,0]
        else:
            raise AttributeError('z coordinate not fould.')
    else:
        dat = get_variable(var)(infile, tidx_start=tidx_start, tidx_end=tidx_end)
        #  TODO: better handle of z for derived variables <14-06-18, Qing Li> #
        z = infile.variables['z'][0,:,0,0]
    return dat, z

def read_ts(infile, var, tidx_start=None, tidx_end=None):
    """Read time series of variable from GOTM file

    :infile: (netCDF4 Dataset) input netCDF file
    :var: (str) variable name
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) time series

    """
    varlist = infile.variables.keys()
    if var in varlist:
        dat = infile.variables[var][tidx_start:tidx_end,0,0]
    else:
        dat = get_variable(var)(infile, tidx_start=tidx_start, tidx_end=tidx_end)
    return dat
