#!/usr/bin/env python3
"""
Qing Li, 20180516
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset, num2date
from testrestool import *

def main():

    # case
    var_list = ['LaTurb', 'LaSL', 'hoLmo']
    ylabel = ['$La_t$', '$La_{SL}$', '$-h_b/L_{MO}$']

    # test
    l_test = True
    l_test = False
    i_test = 1
    j_test = 2

    if l_test:
        plot_ts_var(case_list[i_test], var_list[j_test], ylabel[j_test])
    else:
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
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/TS_'+var+'.png'

    # use the first in the list as the reference case
    data0 = dataroot+'/KPP-CVMix_VR1m_DT60s/gotm_out.nc'

    # read data
    tidx_start = 1
    tidx_end = None
    infile0 = Dataset(data0, 'r')
    fld0 = read_ts(infile0, var, tidx_start=tidx_start, tidx_end=tidx_end)
    # remove negative values
    fld0 = np.ma.array(fld0, mask=(fld0<=0))

    # time
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[tidx_start:tidx_end], units=nctime0.units, calendar=t_cal)

    # figure size
    f = plt.gcf()
    fig_width = 6
    fig_height = 2.5
    f.set_size_inches(fig_width, fig_height)

    # plot figure
    ax = plt.gca()
    ax.plot(dttime0, fld0, '-k', linewidth=1.5)
    ax.set_ylabel(ylabel)
    # ax.set_ylim([-500, 500])
    ax.set_yscale('log')
    # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

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
