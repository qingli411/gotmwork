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

    var_list = ['Epot', 'dPEdt']
    ylabel = ['$E_P$', '$dE_P/dt$']
    color = ['k', 'b', 'r', 'g', 'm', 'c', 'y']

    # loop over all cases
    nc = len(case_list)
    nv = len(var_list)
    nm = len(turbmethod_list)
    for i in np.arange(nc):
        case = case_list[i]
        print(case)
        for j in np.arange(nv):
            var = var_list[j]
            plot_ts_cmp_turbmethods(turbmethod_list, case, var, ylabel[j], color, legend_list)

def plot_ts_cmp_turbmethods(turbmethod_list, case, var, ylabel, color, legend_list):

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/TS_cmp_turbmethods_'+var+'.png'

    # number of turbulence methods
    nm = len(turbmethod_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod_list[0]+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0 = read_ts(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure size
    f, axarr = plt.subplots(2, sharex=True)
    fig_width = 6
    fig_height = 4
    f.set_size_inches(fig_width, fig_height)

    # plot figure
    axarr[0].plot(dttime0, fld0, '-k', linewidth=1.5, label=legend_list[0])
    axarr[0].set_ylabel(ylabel)
    axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # loop over other turbmethods
    for i in np.arange(nm-1):
        j = i+1
        data1 = dataroot+'/'+turbmethod_list[j]+'_VR1m_DT60s/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        fld1 = read_ts(infile1, var)
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)
        axarr[1].plot(dttime1, fld1-fld0, '-', color=color[j], linewidth=1.5, label=legend_list[j])

    # ylabel
    axarr[1].set_ylabel(ylabel+' Diff.')
    axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # horizontal reference line
    axarr[1].axhline(y=0, color='k', linestyle='-', linewidth=0.75)

    # legend
    axarr[0].legend(fontsize='small')
    axarr[1].legend(fontsize='small', loc=3)

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
