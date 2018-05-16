#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
tooldir_default = os.environ['HOME']+'/models/gotm/gotmwork/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import *

def main():
    forc_scheme = "COREII"
    turb_scheme1 = "KPP-CVMix"
    # turb_scheme2 = "KPPLT-ENTR"
    turb_scheme2 = "OSMOSIS"
    analysis = "mldMean"
    method = "deltaR"
    month = "DecY2"

    # input directory
    dir_in = os.environ['HOME']+'/work/gotmfigures/'+forc_scheme

    # output directory
    dir_out = os.environ['HOME']+'/work/gotmfigures/'+forc_scheme

    # starting and ending dates
    date_start_analy, date_end_analy = get_analy_dates(month)

    # input data
    dat1name = forc_scheme+'_'+turb_scheme1+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.npz'
    dat2name = forc_scheme+'_'+turb_scheme2+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.npz'

    # load data
    d1 = np.load(dir_in+'/'+dat1name)
    d2 = np.load(dir_in+'/'+dat2name)
    dat1 = d1['dat']
    dat2 = d2['dat']
    dat = (dat2-dat1)/dat1*100.0
    lon = d1['lon']
    lat = d1['lat']

    # output figure name
    figname = forc_scheme+'_'+turb_scheme2+'-'+turb_scheme1+'_'+analysis+'_'+method+'_'+date_start_analy+'-'+date_end_analy+'.png'

    # plot figures
    plt.figure()
    outfig = dir_out+'/'+figname
    plot_map_scatter(lon, lat, dat, outfig, vmax=40, vmin=-40, cmap='RdBu_r')
    plt.close()

def get_analy_dates(month):
    adate_start = {
            "test":  "20090101",
            "Dec":   "20081201",
            "Jun":   "20090601",
            "DecY2": "20091201"
            }
    adate_end = {
            "test":  "20090331",
            "Dec":   "20081231",
            "Jun":   "20090630",
            "DecY2": "20091231"
            }
    return [adate_start.get(month, None), adate_end.get(month, None)]

if __name__ == "__main__":
    main()
