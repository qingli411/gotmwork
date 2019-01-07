#!/usr/bin/env python3
"""
Qing Li, 20180718
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset, num2date
from testiotool import *

def main():

    cmp_list = ['dampV_none', 'dampV_10d', 'dampV_1d']
    l_interp = [False, False, False]

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
            c_max = pfl_cmax_list[j,i]
            c_min = pfl_cmin_list[j,i]
            d_max = pfl_dmax_list[j,i]
            for k in np.arange(nm):
                turbmethod = turbmethod_list[k]
                plot_pfl_cmp(cmp_list, l_interp, case, turbmethod, var, c_max, c_min, d_max, depth)

def plot_pfl_cmp(cmp_list, l_interp, case, turbmethod, var, c_max, c_min, d_max, depth):

    # input data
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_'+turbmethod+'_'+var+'.png'

    # number of dz dt cases
    ncmp = len(cmp_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod+'_'+cmp_list[0]+'/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0, z0 = gotm_read_pfl(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure
    fig_width = 6
    fig_height = 2+2*(ncmp-1)

    # plot figure
    f, axarr = plt.subplots(ncmp, sharex=True)
    f.set_size_inches(fig_width, fig_height)

    # contour levels
    c_int = (c_max-c_min)/20
    levels0 = np.arange(c_min, c_max+c_int, c_int)
    d_int = d_max/10
    levels1 = np.arange(-d_max, d_max+d_int, d_int)

    # contourf plot
    im0 = axarr[0].contourf(dttime0, z0, np.transpose(fld0), levels0, extend='both', cmap='rainbow')
    axarr[0].set_ylabel('Depth (m)')
    axarr[0].set_ylim([depth, 0])
    title0 = case+' '+var+' '+cmp_list[0]
    axarr[0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over other turbmethods
    for i in np.arange(ncmp-1):
        j = i+1
        data1 = dataroot+'/'+turbmethod+'_'+cmp_list[j]+'/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        fld1_tmp, z1_tmp = gotm_read_pfl(infile1, var)
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)

        # interpolate to z0
        if l_interp[j]:
            nt = fld0.shape[0]
            fld1 = np.zeros(fld0.shape)
            for k in np.arange(nt):
                fld1[k,:] = np.interp(z0, z1_tmp, fld1_tmp[k,:])
            z1 = z0
        else:
            fld1 = fld1_tmp
            z1 = z1_tmp

        im1 = axarr[j].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        axarr[j].set_ylabel('Depth (m)')
        axarr[j].set_ylim([depth, 0])
        title1 = 'Diff. '+cmp_list[j]+' $-$ '+cmp_list[0]
        axarr[j].set_title(title1, fontsize=10)
        cb1 = plt.colorbar(im1, ax=axarr[j])
        cb1.formatter.set_powerlimits((-2, 2))
        cb1.update_ticks()

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
