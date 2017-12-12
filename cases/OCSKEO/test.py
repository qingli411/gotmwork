#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import datetime

def main():
    dateend = '20160411'
    dend = datetime.datetime.strptime(dateend, '%Y%m%d')

    oname = 'sst.dat'
    keyre_datetime = r'^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}'
    # keyre_datetime = r'^[0-9]{4}-[0-9]{2}-[0-9]{2}'
    with open(oname) as fdat:
        for line in reversed(fdat.readlines()):
            if re.search(keyre_datetime, line):
                print(line)
                break
    date, time, val = line.split()
    print(date)
    print(time)
    dtime_new = datetime.datetime.strptime(date+' '+time, '%Y-%m-%d %H:%M:%S')
    isnew = dtime_new > dend
    print(isnew)

if __name__ == "__main__":
    main()
