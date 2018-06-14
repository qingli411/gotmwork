#!/usr/bin/env python3
"""
Qing Li, 20180508
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset, num2date

def main():

    # case
    case_list = ['OCSPapa_20130621-20131201']
    var_list = ['temp', 'salt']
    turbmethod_list = ['KPP-CVMix',
                       'KPPLT-EFACTOR',
                       'KPPLT-ENTR',
                       'OSMOSIS',
                       'EPBL',
                       'SMC',
                       'SMCLT']
    cmax_list = np.array([[18],
                          [33.8]])
    cmin_list = np.array([[4],
                          [32.2]])
    dmax_list = np.array([[3],
                          [0.2]])

    # loop over all cases
    nc = len(case_list)
    nv = len(var_list)
    nm = len(turbmethod_list)
    for i in np.arange(nc):
        case = case_list[i]
        print(case)
        for j in np.arange(nv):
            var = var_list[j]
            c_max = cmax_list[j,i]
            c_min = cmin_list[j,i]
            d_max = dmax_list[j,i]
            plot_pfl_cmp_turbmethods_obs(turbmethod_list, case, var, c_max, c_min, d_max)

def plot_pfl_cmp_turbmethods_obs(turbmethod_list, case, var, c_max, c_min, d_max):

    # input data directory
    dataroot = '/Users/qingli/work/gotmrun/TEST_RES/'+case

    # output figure name
    figdir = '/Users/qingli/work/gotmfigures/TEST_RES/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_turbmethods_obs_'+var+'.png'

    # number of turbulence methods
    nm = len(turbmethod_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod_list[0]+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    ncvar0 = infile0.variables[var+'_obs']
    fld0 = ncvar0[:,:,0,0]
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)
    z0 = read_z(infile0, ncvar0)

    # figure
    fig_width = 6
    fig_height = 2+2*(nm-1)

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
    axarr[0].set_ylim([-120, 0])
    title0 = case+' '+var+': OBS'
    axarr[0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over other turbmethods
    for i in np.arange(nm):
        j = i+1
        data1 = dataroot+'/'+turbmethod_list[i]+'_VR1m_DT60s/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        ncvar1 = infile1.variables[var]
        fld1 = ncvar1[:,:,0,0]
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)
        z1 = read_z(infile1, ncvar1)

        im1 = axarr[j].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        axarr[j].set_ylabel('Depth (m)')
        axarr[j].set_ylim([-120, 0])
        title1 = 'Diff. from OBS: '+turbmethod_list[i]
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

def read_z(infile, ncvar):
    """Read the z coordinate of a variable in GOTM output assuming fixed layer

    :infile: (netCDF4 Dataset) input netCDF file
    :ncvar: (netCDF variable) variable
    :returns: z coordinate (negative below the surface)

    """
    # choose veritcal coordinate
    varlist = infile.variables.keys()
    # GOTM output (fixed z)
    try:
        coord = ncvar.coordinates
    except AttributeError:
        coord = 'v4'
    if 'zi' in coord:
        z = infile.variables['zi'][0,:,0,0]
    elif 'z' in coord:
        z = infile.variables['z'][0,:,0,0]
    else:
        z = infile.variables['z'][:]
    return z

if __name__ == "__main__":
    main()
