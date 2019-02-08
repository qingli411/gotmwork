#!/usr/bin/env python3
"""
Qing Li, 20181003
"""

import matplotlib.dates as mdates
from netCDF4 import Dataset, num2date
from testrestool import *

def main():

    # variable
    var = 'PE'
    # loop over all cases
    # nc = len(case_list)
    # nm = len(turbmethod_list)
    if len(sys.argv) == 1:
        args_list = np.arange(nc)
    else:
        args_list = sys.argv[1:]

    for i in args_list:
        i = int(i)
        case = case_list[i]
        depth = depth_list[i]
        print(case)
        plot_ts_cmp_dz_dt(case, var, depth)

def plot_ts_cmp_dz_dt(case, var, depth):

    dz = np.zeros(nzt)
    dt = np.zeros(nzt)
    dz_str, dt_str = dzdt_list[0].split('_')
    dz[0] = float(dz_str.replace('VR','').replace('m',''))
    dt[0] = float(dt_str.replace('DT','').replace('s',''))

    # input data directory
    dataroot = dir_in+'/'+case
    # paths of files
    paths = [dataroot+'/'+turbmethod_list[i]+'_VR1m_DT60s/gotm_out.nc' for i in range(nm)]
    # initialize dataset
    data = GOTMOutputDataSet(paths=paths, keys=turbmethod_list)

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/IPE_cmp_dzdt_'+var+'.png'

    # figure
    nrow = (nm+2)//2
    fig_width = 12
    fig_height = 3+2*(nrow-1)

    # plot figure
    height_ratios = [1]*nrow
    height_ratios.append(0.15)
    width_ratios = [1, 1, 0.05]
    f, axarr = plt.subplots(nrow, 2)
    f.set_size_inches(fig_width, fig_height)

    # panel a
    gotmdata0 = data.cases['KPP-CVMix']
    ts  = gotmdata0.read_timeseries(var, depth=depth)
    ts0 = ts.data
    dttime0 = ts.time
    par1 = axarr[0, 0].twinx()
    par1.plot(dttime0, ts0, color='lightgray', linewidth=2.5)
    par1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    par1.set_ylabel('$\Delta PE$ (J)', fontsize=12)
    for i in np.arange(nm):
        gotmdata1 = data.cases[turbmethod_list[i]]
        ts  = gotmdata1.read_timeseries(var, depth=depth)
        ts1 = ts.data
        dttime1 = ts.time
        axarr[0, 0].plot(dttime1, ts1-ts0, color=tm_color[i], linewidth=1)
    axarr[0, 0].set_ylabel('$\Delta PE -\Delta PE_r$ (J)', fontsize=12)
    axarr[0, 0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axarr[0, 0].autoscale(enable=True, axis='x', tight=True)
    axarr[0, 0].set_zorder(par1.get_zorder()+1)
    axarr[0, 0].patch.set_visible(False)
    axarr[0, 0].text(0.04, 0.92, '(a)', transform=axarr[0, 0].transAxes, fontsize=16,
                     fontweight='bold', va='top')
    axarr[0, 0].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))

    # panel b-l
    # loop over other turbmethods
    for i in np.arange(nm):
        n = icol_2col[i]
        m = irow_2col[i]
        # paths of files
        tm_paths = [dataroot+'/'+turbmethod_list[i]+'_'+dzdt_list[k]+'/gotm_out.nc'
                    for k in range(nzt)]
        # initialize dataset
        tm_data = GOTMOutputDataSet(paths=tm_paths, keys=dzdt_list)

        # base case
        gotmdata0 = tm_data.cases['VR1m_DT60s']
        fld0 = gotmdata0.read_timeseries(var, depth=depth).data
        dfld0 = fld0[-1] - fld0[0]
        error_dzdt = np.zeros(nzt)
        # loop over other cases
        for ii in np.arange(nzt-1):
            j = ii+1
            gotmdata1 = tm_data.cases[dzdt_list[j]]
            fld1 = gotmdata1.read_timeseries(var, depth=depth).data
            # compute percentage error
            error_dzdt[j] = np.sqrt(((fld1-fld0)**2).mean())/abs(dfld0)*100
            # get coordinate
            dz_str, dt_str = dzdt_list[j].split('_')
            dz[j] = float(dz_str.replace('VR','').replace('m',''))
            dt[j] = float(dt_str.replace('DT','').replace('s',''))

        # plt.plot(dz[0], error3_dzdt[0], 'ko')
        axarr[m, n].plot(dz[0:9:4], error_dzdt[0:9:4], '--kx', linewidth=1.5)
        axarr[m, n].plot(dz[1:10:4], error_dzdt[1:10:4], '--rx', linewidth=1.5)
        axarr[m, n].plot(dz[2:11:4], error_dzdt[2:11:4], '--bx', linewidth=1.5)
        axarr[m, n].plot(dz[3:12:4], error_dzdt[3:12:4], '--gx', linewidth=1.5)
        axarr[m, n].axhline(0, color='black', linewidth=0.75)
        axarr[m, n].set_xlabel('$\Delta z$ (m)', fontsize=12)
        axarr[m, n].set_ylabel('% error in $PE$', fontsize=12)
        axarr[m, n].set_xlim(0,11)
        axarr[m, n].text(0.04, 0.92, labels_2col[i], transform=axarr[m, n].transAxes,
                         fontsize=16, fontweight='bold', va='top')
        axarr[m, n].text(0.98, 1.15, legend_list[i], color=tm_color[i],
                         transform=axarr[m, n].transAxes, fontsize=13, fontweight='bold',
                         va='top', ha='right')
        axarr[m, n].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axarr[m, n].set_ylim(0, 4)

     # reduce margin
    plt.tight_layout()

    # save figure
    plt.savefig(figname, dpi = 300)

    # close figure
    plt.close()

if __name__ == "__main__":
    main()
