#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import sys
tooldir_default = os.environ['GOTMWORK_ROOT']+'/tools'
sys.path.append(os.environ.get('tooldir', tooldir_default))
from gotmtool import *
from core2tool import *

def main():
    # turb_scheme0 = 'KPP-CVMix'
    # turb_scheme1_list = ['KPPLT-ENTR', 'OSMOSIS', 'EPBL', 'SMC']
    turb_scheme0 = 'KPP-CVMix'
    turb_scheme1_list = ['KPP-ROMS', 'KPP-CVMix-nDC']
    # month_list = ['Dec', 'Feb', 'Jun', 'Aug', 'DecY2']
    month_list = ['Aug']
    analysis = 'mldMean'
    method = 'deltaR'
    nturb = len(turb_scheme1_list)
    nmon = len(month_list)
    # plot scheme0
    levels = [-500, -450, -400, -350, -300, -250, -200, -180, -160, -140, -120, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0]
    print(turb_scheme0)
    for j in np.arange(nmon):
        month = month_list[j]
        print('  {}'.format(month))
        plot_map(turb_scheme0, month, analysis, method, units='m', levels=levels, vmax=0, vmin=-500)

    # plot differences between scheme1 and scheme0
    for i in np.arange(nturb):
        turb_scheme1 = turb_scheme1_list[i]
        print(turb_scheme1)
        for j in np.arange(nmon):
            month = month_list[j]
            print('  {}'.format(month))
            plot_map_diff(turb_scheme0, turb_scheme1, month, analysis, method, units='%', vmax=40, vmin=-40)

if __name__ == "__main__":
    main()
