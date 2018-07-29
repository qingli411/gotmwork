"""
Shared functions.

Qing Li, 20180517
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
tooldir_default = os.environ['GOTMWORK_ROOT']+'/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import *

#--------------------------------
# common variables
#--------------------------------
# forcing scheme
forc_scheme = 'COREII'
# input directory
dir_in = os.environ['GOTMFIG_ROOT']+'/'+forc_scheme
# output directory
dir_out = os.environ['GOTMFIG_ROOT']+'/'+forc_scheme

#--------------------------------
# plot the field in a map
#--------------------------------
def plot_map(turb_scheme0, month, analysis, method, units, vmax=None, vmin=None, levels=None):

    # starting and ending dates
    date_start_analy, date_end_analy = get_analy_dates(month)

    # input data
    dat0name = forc_scheme+'_'+turb_scheme0+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.npz'

    # load data
    d0 = np.load(dir_in+'/'+dat0name)
    dat0 = d0['dat']
    dat = dat0
    # for i in np.arange(len(dat)):
        # print(dat[i])
    dat[dat==-499.5] = np.nan
    lon = d0['lon']
    lat = d0['lat']

    # output figure name
    figname = forc_scheme+'_'+turb_scheme0+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.png'

    # plot figures
    plt.figure()
    outfig = dir_out+'/'+figname
    plot_map_scatter(lon, lat, dat, units=units, levels=levels, vmax=vmax, vmin=vmin)
    # save figure
    plt.savefig(outfig, dpi=300)
    plt.close()

#--------------------------------
# plot difference
#--------------------------------
def plot_map_diff(turb_scheme0, turb_scheme1, month, analysis, method, units=None, vmax=None, vmin=None):

    # starting and ending dates
    date_start_analy, date_end_analy = get_analy_dates(month)

    # input data
    dat0name = forc_scheme+'_'+turb_scheme0+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.npz'
    dat1name = forc_scheme+'_'+turb_scheme1+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.npz'

    # load data
    d0 = np.load(dir_in+'/'+dat0name)
    d1 = np.load(dir_in+'/'+dat1name)
    dat0 = d0['dat']
    dat1 = d1['dat']
    if units == '%':
        dat = (dat1-dat0)/dat0*100.0
    else:
        dat = dat1-dat0
    lon = d0['lon']
    lat = d0['lat']

    # output figure name
    figname = forc_scheme+'_'+turb_scheme1+'-'+turb_scheme0+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.png'

    # plot figures
    plt.figure()
    outfig = dir_out+'/'+figname
    plot_map_scatter(lon, lat, dat, units=units, vmax=vmax, vmin=vmin, cmap='RdBu_r')
    # save figure
    plt.savefig(outfig, dpi=300)
    plt.close()

#--------------------------------
# get the starting and ending date for analyses
#--------------------------------
def get_analy_dates(month):
    adate_start = {
            "test":  "20080601",
            "Feb":   "20090201",
            "Jun":   "20080601",
            "Jul":   "20080701",
            "Aug":   "20080801",
            "Dec":   "20081201",
            }
    adate_end = {
            "test":  "20080630",
            "Feb":   "20090228",
            "Jun":   "20080630",
            "Jul":   "20080731",
            "Aug":   "20080831",
            "Dec":   "20081231",
            }
    return [adate_start.get(month, None), adate_end.get(month, None)]

