#!/usr/bin/env python3
"""
Qing Li, 20180508
"""

from netCDF4 import Dataset, num2date
from testrestool import *

def main():

    # loop over all cases
    # nc = len(case_list)
    # nv = len(pfl_list)
    # nm = len(turbmethod_list)
    if len(sys.argv) == 1:
        args_list = np.arange(nc)
    else:
        args_list = sys.argv[1:]

    for i in args_list:
        i = int(i)
        case = case_list[i]
        title = title_list[i]
        depth = depth_list[i]
        print(case)
        for j in np.arange(nv):
            var = pfl_list[j]
            c_max = pfl_cmax_list[j,i]
            c_min = pfl_cmin_list[j,i]
            d_max = pfl_dmax_list[j,i]
            for k in np.arange(nm):
                turbmethod = turbmethod_list[k]
                plot_pfl_cmp_dz_dt(case, title, turbmethod, var, c_max, c_min, d_max, depth)

def plot_pfl_cmp_dz_dt(case, title, turbmethod, var, c_max, c_min, d_max, depth):

    # input data
    dataroot = dir_in+'/'+case
    # paths of files
    paths = [dataroot+'/'+turbmethod+'_'+dzdt_list[i]+'/gotm_out.nc' for i in range(nzt)]
    # initialize dataset
    data = GOTMOutputDataSet(paths=paths, keys=dzdt_list)

    # output figure name
    figdir = dir_out+'/'+case
    os.makedirs(figdir, exist_ok=True)
    figname = figdir+'/Pfl_cmp_'+turbmethod+'_'+var+'.png'

    # figure
    fig_width = 12
    fig_height = 6.75

    # plot figure
    f, axarr = plt.subplots(int(nzt/3), 3, sharex=True)
    f.set_size_inches(fig_width, fig_height)

    # contour levels
    c_int = (c_max-c_min)/20
    levels0 = np.arange(c_min, c_max+c_int, c_int)
    d_int = d_max/10
    levels1 = np.arange(-d_max-0.5*d_int, d_max+d_int, d_int)
    cb_ticks = np.arange(-d_max, d_max+d_int*2, d_int*2)

    # contourf plot
    gotmdata0 = data.cases['VR1m_DT60s']
    dttime0 = num2date(gotmdata0.time, units=gotmdata0.time_units, calendar=gotmdata0.time_calendar)
    prfl = gotmdata0.read_profile(var)
    fld0 = prfl.data
    z0 = prfl.z
    im0 = axarr[0, 0].contourf(dttime0, z0, np.transpose(fld0), levels0, extend='both', cmap='rainbow')
    axarr[0, 0].set_ylabel('Depth (m)')
    axarr[0, 0].set_ylim([depth, 0])
    title0 = title+' '+var+' '+dzdt_list[0]
    axarr[0, 0].set_title(title0, fontsize=10)

    # loop over other turbmethods
    for i in np.arange(nzt-1):
        j = i+1
        n, m = divmod(j, 4)
        gotmdata1 = data.cases[dzdt_list[j]]
        prfl = gotmdata1.read_profile(var)
        fld1_tmp = prfl.data
        z1_tmp = prfl.z
        dttime1 = num2date(gotmdata1.time, units=gotmdata1.time_units, calendar=gotmdata1.time_calendar)

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

        im1 = axarr[m, n].contourf(dttime1, z1, np.transpose(fld1-fld0), levels1, extend='both', cmap='RdBu_r')
        if n == 0:
            axarr[m, n].set_ylabel('Depth (m)')
        axarr[m, n].set_ylim([depth, 0])
        title1 = 'Diff. '+dzdt_list[j]+' $-$ '+dzdt_list[0]
        axarr[m, n].set_title(title1, fontsize=10)

    # auto adjust the x-axis label
    plt.gcf().autofmt_xdate()

    # reduce margin
    plt.tight_layout()

    plt.subplots_adjust(bottom=0.095, right=0.9)
    cax0 = plt.axes([0.85, 0.55, 0.1, 0.4])
    cax0.set_visible(False)
    cb0 = plt.colorbar(im0, ax=cax0)
    cb0.formatter.set_powerlimits((-2, 2))
    cb0.update_ticks()
    cax1 = plt.axes([0.85, 0.08, 0.1, 0.4])
    cax1.set_visible(False)
    cb1 = plt.colorbar(im1, ax=cax1, ticks=cb_ticks)
    cb1.formatter.set_powerlimits((-2, 2))
    cb1.update_ticks()


    # save figure
    plt.savefig(figname, dpi = 300)

    # close figure
    plt.close()


if __name__ == "__main__":
    main()
