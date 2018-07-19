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
    i_test = 5
    j_test = 2

    case = case_list[i_test]
    coord1 = 'LaSL'
    coord2 = 'hoLmo'
    var = 'dPEdt'
    label1 = '$\mathrm{La}_{SL}$'
    label2 = '$-h_b/L_{MO}$'

    # input data directory
    dataroot = dir_in+'/'+case

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/SCT_'+var+'.png'


    # read data
    tidx_start = 1
    tidx_end = None
    data0 = dataroot+'/'+turbmethod_list[0]+'_VR1m_DT60s/gotm_out.nc'
    infile0 = Dataset(data0, 'r')
    time = infile0.variables['time'][:]
    nt = time.shape[0]
    nm = len(turbmethod_list)
    print(nt)
    print(nm)
    fld0 = np.zeros([nm,nt-1])
    xx0 = np.zeros([nm,nt-1])
    yy0 = np.zeros([nm,nt-1])
    for i in np.arange(nm):
        turbmethod = turbmethod_list[i]
        data0 = dataroot+'/'+turbmethod+'_VR1m_DT60s/gotm_out.nc'
        infile0 = Dataset(data0, 'r')
        fld0[i,:] = gotm_read_ts(infile0, var, tidx_start=tidx_start, tidx_end=tidx_end)
        xx0[i,:] = gotm_read_ts(infile0, coord1, tidx_start=tidx_start, tidx_end=tidx_end)
        yy0[i,:] = gotm_read_ts(infile0, coord2, tidx_start=tidx_start, tidx_end=tidx_end)
    # remove negative values
    # fld0 = np.ma.array(fld0, mask=(fld0<=0))

    # figure size
    f = plt.gcf()
    fig_width = 6
    fig_height = 5
    f.set_size_inches(fig_width, fig_height)

    # plot figure
    cmap = 'rainbow'
    mfld0 = np.mean(fld0, axis=0)
    dat = np.std(fld0, axis=0)/abs(mfld0)
    # dat = np.std(fld0, axis=0)/mfld0
    datn = np.ma.array(dat, mask=(mfld0>0))
    datp = np.ma.array(dat, mask=(mfld0<0))
    xx = np.mean(xx0, axis=0)
    yy = np.mean(yy0, axis=0)
    # sc = plt.scatter(xx, yy, marker='s', s=5, c=datp, cmap=plt.cm.get_cmap(cmap), vmin=0, vmax=0.5)
    sc = plt.scatter(xx, yy, marker='o', s=5, c=dat, cmap=plt.cm.get_cmap(cmap), vmin=0, vmax=0.5)
    # sc = plt.scatter(xx, yy, marker='o', s=5, c=dat, cmap=plt.cm.get_cmap(cmap))
    plt.colorbar(sc)
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e-1, 1e1])
    plt.ylim([1e-3, 1e3])
    # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # reduce margin
    plt.tight_layout()

    # save figure
    plt.savefig(figname, dpi = 300)

    # close figure
    plt.close()

if __name__ == "__main__":
    main()
