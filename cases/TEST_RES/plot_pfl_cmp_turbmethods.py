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

    l_ens = False
    l_twocol = True
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
            if l_ens:
                plot_pfl_cmp_turbmethods_ens(case, var, c_max, c_min, d_max, depth)
            elif l_twocol:
                plot_pfl_cmp_turbmethods_2col(case, var, c_max, c_min, d_max, depth)
            else:
                plot_pfl_cmp_turbmethods(case, var, c_max, c_min, d_max, depth)

def plot_pfl_cmp_turbmethods(case, var, c_max, c_min, d_max, depth):

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_turbmethods_'+var+'.png'

    # number of turbulence methods
    nm = len(turbmethod_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod_list[0]+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0, z0 = gotm_read_pfl(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure
    fig_width = 6
    fig_height = 2+2*(nm-1)

    # plot figure
    f, axarr = plt.subplots(nm, sharex=True)
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
    title0 = case+' '+var+' '+turbmethod_list[0]
    axarr[0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over other turbmethods
    for i in np.arange(nm-1):
        j = i+1
        data1 = dataroot+'/'+turbmethod_list[j]+'_VR1m_DT60s/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        fld1, z1 = gotm_read_pfl(infile1, var)
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)

        im1 = axarr[j].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        axarr[j].set_ylabel('Depth (m)')
        axarr[j].set_ylim([depth, 0])
        title1 = 'Diff. '+legend_list[j]+' $-$ '+legend_list[0]
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

def plot_pfl_cmp_turbmethods_2col(case, var, c_max, c_min, d_max, depth):

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_turbmethods_'+var+'.png'

    # number of turbulence methods
    nm = len(turbmethod_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod_list[0]+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0, z0 = gotm_read_pfl(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)

    # figure
    nrow = (nm+1)//2
    fig_width = 12
    fig_height = 2+2*(nrow-1)

    # plot figure
    f, axarr = plt.subplots(nrow, 2, sharex='col', sharey='row')
    f.set_size_inches(fig_width, fig_height)

    # contour levels
    c_int = (c_max-c_min)/20
    levels0 = np.arange(c_min, c_max+c_int, c_int)
    d_int = d_max/10
    levels1 = np.arange(-d_max, d_max+d_int, d_int)

    # contourf plot
    im0 = axarr[0, 0].contourf(dttime0, z0, np.transpose(fld0), levels0, extend='both', cmap='rainbow')
    axarr[0, 0].set_ylabel('Depth (m)')
    axarr[0, 0].set_ylim([depth, 0])
    title0 = case+' '+var+' '+turbmethod_list[0]
    axarr[0, 0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0, 0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over other turbmethods
    for i in np.arange(nm-1):
        j = i+1
        n = j%2
        m = j//2
        data1 = dataroot+'/'+turbmethod_list[j]+'_VR1m_DT60s/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        fld1, z1 = gotm_read_pfl(infile1, var)
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)

        im1 = axarr[m, n].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        if n == 0:
            axarr[m, n].set_ylabel('Depth (m)')
        axarr[m, n].set_ylim([depth, 0])
        title1 = 'Diff. '+legend_list[j]+' $-$ '+legend_list[0]
        axarr[m, n].set_title(title1, fontsize=10)
        cb1 = plt.colorbar(im1, ax=axarr[m, n])
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

def plot_pfl_cmp_turbmethods_ens(case, var, c_max, c_min, d_max, depth):

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_turbmethods_ens_'+var+'.png'

    # number of turbulence methods
    nm = len(turbmethod_list)

    # use ensemble median as the reference
    for i in np.arange(nm):
        data0 = dataroot+'/'+turbmethod_list[i]+'_VR1m_DT60s/gotm_out.nc'
        # read data
        infile0 = Dataset(data0, 'r')
        if i == 0:
            nctime0 = infile0.variables['time']
            t_cal = 'standard'
            dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)
            nt = nctime0[:].shape[0]
            tmp, z0 = gotm_read_pfl(infile0, var)
            nz = z0.shape[0]
            tmp0 = np.zeros([nm, nt, nz])
            tmp0[i,:,:], z0 = gotm_read_pfl(infile0, var)
        else:
            tmp0[i,:,:], z0 = gotm_read_pfl(infile0, var)
    # get the median
    fld0 = np.median(tmp0, axis=0)

    # figure
    fig_width = 6
    fig_height = 2+2*nm

    # plot figure
    f, axarr = plt.subplots(nm+1, sharex=True)
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
    title0 = case+' '+var+' Ens. Median'
    axarr[0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over turbmethods
    for i in np.arange(nm):
        j = i+1
        data1 = dataroot+'/'+turbmethod_list[i]+'_VR1m_DT60s/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        fld1, z1 = gotm_read_pfl(infile1, var)
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)

        im1 = axarr[j].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        axarr[j].set_ylabel('Depth (m)')
        axarr[j].set_ylim([depth, 0])
        title1 = 'Diff. '+legend_list[i]+' $-$ Ens. Median'
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
