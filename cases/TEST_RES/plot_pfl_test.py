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
tooldir_default = os.environ['HOME']+'/models/gotm/gotmwork/tools'
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
    var_list = ['temp', 'salt', 'buoy', 'buoyancy', 'spice']
    turbmethod_list = ['KPP-CVMix',
                       'KPPLT-EFACTOR',
                       'KPPLT-ENTR',
                       'OSMOSIS',
                       'EPBL',
                       'SMC',
                       'SMCLT']
    depth_list = np.array([-200, -240, -120, -150, -120, -400])

    i_test = 0
    j_test = 4
    k_test = 0
    c_max = 1.4e-2
    c_min = -0.4e-2
    case = case_list[i_test]
    depth = depth_list[i_test]
    turbmethod = turbmethod_list[k_test]
    print(case)
    var = var_list[j_test]
    plot_pfl(turbmethod, case, var, depth, c_max=c_max, c_min=c_min)

def plot_pfl(turbmethod, case, var, depth, c_max=None, c_min=None):

    # input data directory
    dataroot = '/Users/qingli/work/gotmrun/TEST_RES/'+case

    # output figure name
    figdir = '/Users/qingli/work/gotmfigures/TEST_RES/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_'+var+'.png'

    data0 = dataroot+'/'+turbmethod+'_VR1m_DT60s/gotm_out.nc'

    # read data
    infile0 = Dataset(data0, 'r')
    fld0 = read_var(infile0, var)
    nctime0 = infile0.variables['time']
    t_cal = 'standard'
    dttime0 = num2date(nctime0[:], units=nctime0.units, calendar=t_cal)
    z0 = read_z(infile0, var)

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

def read_var(infile, var, tidx_start=None, tidx_end=None):
    """Read variable from infile

    :infile: (netCDF4 Dataset) input netCDF file
    :var: (str) variable name
    :tidx_start: (int, optional) starting index
    :tidx_end: (int, optional) ending index
    :returns: (numpy array) variable

    """
    varlist = infile.variables.keys()
    if var in varlist:
        dat = infile.variables[var][tidx_start:tidx_end,:,0,0]
    else:
        dat = get_variable(var)(infile, tidx_start=tidx_start, tidx_end=tidx_end)
    return dat


def read_z(infile, var):
    """Read the z coordinate of a variable in GOTM output assuming fixed layer

    :infile: (netCDF4 Dataset) input netCDF file
    :ncvar: (str) variable name
    :returns: z coordinate (negative below the surface)

    """
    # choose veritcal coordinate
    varlist = infile.variables.keys()
    if var in varlist:
        # nc variable
        ncvar = infile.variables[var]
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
    else:
        z = infile.variables['z'][0,:,0,0]
    return z

if __name__ == "__main__":
    main()
