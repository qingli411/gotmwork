#!/usr/bin/env python3
"""
Qing Li, 20180516
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import datetime
from netCDF4 import Dataset, num2date
tooldir_default = os.environ['GOTMWORK_ROOT']+'/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import *

def main():

    # case
    case_list = ['OSMOSIS_winter',
                 'OSMOSIS_spring',
                 'OCSPapa_20130621-20131201',
                 'COREII_LAT2_LON234_20080615-20081231',
                 'COREII_LAT10_LON86_20080615-20081231',
                 'COREII_LAT-54_LON254_20080915-20090915']
    var_list = ['La']
    ylabel = ['$La_t$']

    # loop over all cases
    nc = len(case_list)
    nv = len(var_list)
    for i in np.arange(nc):
        case = case_list[i]
        print(case)
        for j in np.arange(nv):
            var = var_list[j]
            plot_ts_var(case, var, ylabel[j])

def plot_ts_var(case, var, ylabel):

    # input data directory
    dataroot = os.environ['GOTMRUN_ROOT']+'/TEST_RES/'+case

    # output figure name
    figdir = os.environ['GOTMFIG_ROOT']+'/TEST_RES/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/TS_LaTurb_'+var+'.png'

    # use the first in the list as the reference case
    data0 = dataroot+'/KPP-CVMix_VR1m_DT60s/gotm_out.nc'

    # read data
    if var == 'La':
        infile0 = Dataset(data0, 'r')
        # fld01x = infile0.variables['u0_stokes'][:,0,0]
        # fld01y = infile0.variables['v0_stokes'][:,0,0]
        # fld01 = np.sqrt(fld01x**2.+fld01y**2.)
        # fld02 = infile0.variables['u_taus'][:,0,0]
        # fld0 = np.sqrt(fld02/fld01)
        fld0 = get_la('LaTurb')(infile0)
        nctime0 = infile0.variables['time']
        t_cal = 'standard'
        dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure size
    f = plt.gcf()
    fig_width = 6
    fig_height = 2.5
    f.set_size_inches(fig_width, fig_height)

    # plot figure
    ax = plt.gca()
    ax.plot(dttime0, fld0, '-k', linewidth=1.5)
    ax.set_ylabel(ylabel)
    ax.set_ylim([0, 1])
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # auto adjust the x-axis label
    plt.gcf().autofmt_xdate()

    # reduce margin
    plt.tight_layout()

    # save figure
    plt.savefig(figname, dpi = 300)

    # close figure
    plt.close()

if __name__ == "__main__":
    main()
