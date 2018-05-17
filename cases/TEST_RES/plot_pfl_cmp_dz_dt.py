#!/usr/bin/env python3
"""
Qing Li, 20180508
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import datetime
from netCDF4 import Dataset, num2date

def main():

    # case
    case_list = ['OSMOSIS_winter',
                 'OSMOSIS_spring',
                 'OCSPapa_20130621-20131201',
                 'COREII_LAT2_LON234_20080615-20081231',
                 'COREII_LAT10_LON86_20080615-20081231',
                 'COREII_LAT-54_LON254_20080915-20090915']
    turbmethod_list = ['KPP-CVMix',
                       'KPPLT-EFACTOR',
                       'KPPLT-ENTR',
                       'OSMOSIS',
                       'EPBL',
                       'SMC',
                       'SMCLT']
    var_list = ['temp', 'salt']
    dzdt_list = ['VR1m_DT60s', 'VR1m_DT1800s', 'VR5m_DT60s', 'VR5m_DT1800s']
    l_interp = [False, False, True, True]
    cmax_list = np.array([[16, 20, 18, 27, 29, 7],
                          [35.9, 35.9, 33.8, 35.1, 34.6, 34.2]])
    cmin_list = np.array([[12, 12, 4, 12, 14, 4],
                          [35.7, 35.6, 32.2, 34.7, 33.3, 33.9]])
    dmax_list = np.array([[1, 1, 1, 1, 1, 1],
                          [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]])

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
            for k in np.arange(nm):
                turbmethod = turbmethod_list[k]
                plot_pfl_cmp_dz_dt(dzdt_list, l_interp, case, turbmethod, var, c_max, c_min, d_max)

def plot_pfl_cmp_dz_dt(dzdt_list, l_interp, case, turbmethod, var, c_max, c_min, d_max):

    # input data
    dataroot = '/Users/qingli/work/gotmrun/TEST_RES/'+case

    # output figure name
    figdir = '/Users/qingli/work/gotmfigures/TEST_RES/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_'+turbmethod+'_'+var+'.png'

    # number of dz dt cases
    nzt = len(dzdt_list)

    # use the first in the list as the reference case
    data0 = dataroot+'/'+turbmethod+'_'+dzdt_list[0]+'/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    ncvar0 = infile0.variables[var]
    fld0 = ncvar0[:,:,0,0]
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)
    z0 = read_z(infile0, ncvar0)

    # figure
    fig_width = 6
    fig_height = 2+2*(nzt-1)

    # plot figure
    f, axarr = plt.subplots(nzt, sharex=True)
    f.set_size_inches(fig_width, fig_height)

    # contour levels
    c_int = (c_max-c_min)/20
    levels0 = np.arange(c_min, c_max+c_int, c_int)
    d_int = d_max/10
    levels1 = np.arange(-d_max, d_max+d_int, d_int)

    # contourf plot
    im0 = axarr[0].contourf(dttime0, z0, np.transpose(fld0), levels0, extend='both', cmap='rainbow')
    axarr[0].set_ylabel('Depth (m)')
    title0 = case+' '+var+' '+dzdt_list[0]
    axarr[0].set_title(title0, fontsize=10)
    cb0 = plt.colorbar(im0, ax=axarr[0])
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()

    # loop over other turbmethods
    for i in np.arange(nzt-1):
        j = i+1
        data1 = dataroot+'/'+turbmethod+'_'+dzdt_list[j]+'/gotm_out.nc'
        infile1 = Dataset(data1, 'r')
        ncvar1 = infile1.variables[var]
        fld1_tmp = ncvar1[:,:,0,0]
        nctime1 = infile1.variables['time']
        dttime1 = num2date(nctime1[:], units=nctime1.units, calendar=t_cal)
        z1_tmp = read_z(infile1, ncvar1)

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
        title1 = 'Diff. '+dzdt_list[j]+' $-$ '+dzdt_list[0]
        axarr[j].set_title(title1, fontsize=10)
        cb1 = plt.colorbar(im1, ax=axarr[j])
        cb1.formatter.set_powerlimits((-2, 2))
        cb1.update_ticks()

    # auto adjust the x-axis label
    plt.gcf().autofmt_xdate()

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
