#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Qing Li, 20171216

import numpy as np
import datetime

def main():
    fname = '~/data/CDIP/Papa/de166p101_201601-201612'
    nheader = 3
    tbin = np.array([100, 22, 18, 16, 14, 12, 10, 8, 6, 2])
    fbin = 1.0/tbin
    print(fbin)
    f = 0.5*(fbin[0:-1]+fbin[1:])
    print(f)
    df = fbin[1:]-fbin[0:-1]
    print(df)

    # load data in txt format
    data = np.loadtxt(fname, skiprows=nheader)
    nt = data.shape[0]
    # unpack data
    time = data[:,0]
    dt_format = '%Y%m%d%H%M'
    strtime = np.array(['{}'.format(int(time[i])) for i in range(nt)])
    dttime = np.array([datetime.datetime.strptime(strtime[i], dt_format)
        for i in range(nt)])
    hs = data[:,1]*1e-2 # cm -> m
    tp = data[:,2]
    spec = data[:,3:]*1e-4 # cm**2 -> m**2
    print(spec.shape)
    hs0 = 4*np.sqrt(np.sum(spec,axis=1))
    print(hs0)

if __name__ == "__main__":
    main()
