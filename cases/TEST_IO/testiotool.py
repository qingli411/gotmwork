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
from gotmtool import gotm_read_pfl, gotm_read_ts

#--------------------------------
# common variables
#--------------------------------
# input directory
dir_in = os.environ['GOTMRUN_ROOT']+'/TEST_IO/'
# output directory
dir_out = os.environ['GOTMFIG_ROOT']+'/TEST_IO/'
# list of cases
case_list = ['OSMOSIS_winter',
             'OSMOSIS_spring',
             'OCSPapa_20130621-20131201',
             'COREII_LAT2_LON234_20080615-20081231',
             'COREII_LAT10_LON86_20080615-20081231',
             'COREII_LAT-54_LON254_20080915-20090915']
# list of depths
depth_list = np.array([-200, -240, -120, -150, -120, -400])
# list of turbulent methods
turbmethod_list = ['KPP-CVMix',
                   'KPPLT-EFACTOR',
                   'KPPLT-ENTR',
                   'KPPLT-RWHGK',
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
               'KPPLT-RWHGK16',
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

