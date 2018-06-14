#!/usr/bin/env python3
"""
Qing Li, 20180508
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset, num2date
from testrestool import *

def main():

    # test
    l_test = True
    l_test = False
    i_test = 5
    j_test = 4
    k_test = 0

    if l_test:
        case = case_list[i_test]
        print(case)
        depth = depth_list[i_test]
        turbmethod = turbmethod_list[k_test]
        var = pfl_list[j_test]
        c_max = pfl_cmax_list[j_test, i_test]
        c_min = pfl_cmin_list[j_test, i_test]
        plot_pfl(turbmethod, case, var, depth, c_max=c_max, c_min=c_min)
    else:
        # loop over all cases
        nc = len(case_list)
        nv = len(pfl_list)
        nm = len(turbmethod_list)
        for i in np.arange(nc):
            case = case_list[i]
            depth = depth_list[i]
            print(case)
            for j in np.arange(nv):
                var = pfl_list[j]
                print(var)
                c_max = pfl_cmax_list[j,i]
                c_min = pfl_cmin_list[j,i]
                for k in np.arange(nm):
                    turbmethod = turbmethod_list[k]
                    plot_pfl(turbmethod, case, var, depth, c_max=c_max, c_min=c_min)

def plot_pfl(turbmethod, case, var, depth, c_max=None, c_min=None):

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_'+turbmethod+'_'+var+'.png'

    data0 = dataroot+'/'+turbmethod+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0, z0 = read_pfl(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure
    fig_width = 6
    fig_height = 3

    # plot figure
    f = plt.gcf()
    f.set_size_inches(fig_width, fig_height)

    # contourf plot
    if c_max and c_min:
        # contour levels
        c_int = (c_max-c_min)/20
        levels0 = np.arange(c_min, c_max+c_int, c_int)
        im0 = plt.contourf(dttime0, z0, np.transpose(fld0), levels0, extend='both', cmap='rainbow')
    else:
        im0 = plt.contourf(dttime0, z0, np.transpose(fld0), extend='both', cmap='rainbow')
    plt.ylabel('Depth (m)')
    plt.ylim([depth, 0])
    title0 = case+' '+var+' '+turbmethod
    plt.title(title0, fontsize=10)
    cb0 = plt.colorbar(im0)
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

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
