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
from gotmanalysis import *

#--------------------------------
# common variables
#--------------------------------
# input directory
# dir_in = os.environ['GOTMRUN_ROOT']+'/TEST_IO/'
dir_in = '/Volumes/Qing_Work/work/gotmrun/TEST_IO/'
# output directory
dir_out = os.environ['GOTMFIG_ROOT']+'/TEST_IO/'
# list of cases
case_list = ['OSMOSIS_winter',
             'OSMOSIS_spring',
             'OCSPapa_20130621-20131201',
             'COREII_LAT2_LON234_20080601-20081231',
             'COREII_LAT10_LON86_20080601-20081231',
             'COREII_LAT-54_LON254_20090101-20090731']
# list of titles
title_list = ['OSMOSIS Winter',
              'OSMOSIS Spring',
              'OCSPapa',
              'COREII LAT2 LON234',
              'COREII LAT10 LON86',
              'COREII LAT-54 LON254']
# list of depths
depth_list = np.array([-200, -480, -100, -150, -120, -400])
# list of turbulent methods
turbmethod_list = ['KPP-CVMix',
                   'KPP-ROMS',
                   'KPPLT-EFACTOR',
                   'KPPLT-ENTR',
                   'KPPLT-RWHGK',
                   'EPBL',
                   'EPBL-LT',
                   'SMC',
                   'SMCLT',
                   'K-EPSILON-SG',
                   'OSMOSIS']
# list of legend for turbulent methods
legend_list = ['KPP-CVMix',
               'KPP-ROMS',
               'KPPLT-VR12',
               'KPPLT-LF17',
               'KPPLT-RWHGK16',
               'ePBL',
               'ePBL-LT',
               'SMC-KC94',
               'SMCLT-H15',
               'k-epsilon',
               'OSMOSIS']
tm_color = ['black',
            'blue',
            'red',
            'orange',
            'purple',
            'skyblue',
            'steelblue',
            'limegreen',
            'green',
            'mediumvioletred',
            'darkgoldenrod']
dzdt_list = ['VR1m_DT60s',
             'VR1m_DT600s',
             'VR1m_DT1800s',
             'VR1m_DT3600s',
             'VR5m_DT60s',
             'VR5m_DT600s',
             'VR5m_DT1800s',
             'VR5m_DT3600s',
             'VR10m_DT60s',
             'VR10m_DT600s',
             'VR10m_DT1800s',
             'VR10m_DT3600s']
l_interp = [False, False, False, False,
            True, True, True, True,
            True, True, True, True]
# list of location
irow_2col = [1, 2, 0, 1, 2, 3, 3, 4, 4, 5, 5]
icol_2col = [0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1]
labels_2col = ['(b)', '(c)', '(g)', '(h)', '(i)', '(d)', '(j)', '(e)', '(k)','(f)','(l)']
# list of profiles
pfl_list = ['temp', 'salt', 'rho', 'buoy', 'spice']
# max and min values for colorbar ([ncase, nprofile])
pfl_cmax_list = np.array([[16, 20, 18, 27, 29, 7],
                      [35.9, 35.9, 33.8, 35.1, 35.0, 34.2],
                      [1027.3, 1027.2, 1026.8, 1026.4, 1024.2, 1027.1],
                      [8.e-3, 8.e-3, 3.2e-2, 2.8e-2, 5.0e-2, 2.3e-3],
                      [1.7e-2, 2.2e-2, -0.5e-2, 2.8e-2, 2.4e-2, -1.2e-2]])
pfl_cmin_list = np.array([[12, 12, 4, 12, 20, 4],
                      [35.7, 35.6, 32.2, 34.7, 33.0, 33.9],
                      [1026.0, 1025.3, 1023.2, 1022.7, 1021.5, 1026.7],
                      [-3.5e-3, -2.5e-3, 8e-3, 4.e-3, 1.5e-2, -3.2e-3],
                      [0.9e-2, 0.6e-2, -2.7e-2, 0, 1.2e-2, -1.8e-2]])
pfl_dmax_list = np.array([[1, 1, 5, 2, 1, 1],
                      [0.05, 0.05, 0.3, 0.05, 0.5, 0.05],
                      [0.1, 0.1, 1.0, 0.1, 0.5, 0.1],
                      [2.e-3, 2.e-3, 1.e-2, 2.e-3, 1.e-2, 1.e-3],
                      [2.e-3, 2.e-3, 1.e-2, 2.e-3, 1.e-2, 1.e-3]])

nc = len(case_list)
nv = len(pfl_list)
nm = len(turbmethod_list)
nzt = len(dzdt_list)
