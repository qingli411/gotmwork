#!/usr/bin/env bash

# echo "Test 1: set the starting and ending dates"
# ../tools/nc2dat_argo -i "/Volumes/Qing_Work/data/obs/Argo/WOD13_PFL_2000.nc" -lat 25.7 -lon 180.2 -r 2 -ds 20090225 -de 20090425 -ot "tprof_test1.dat" -os "sprof_test1.dat"

echo "Test 2: set the starting date and date range"
../tools/nc2dat_argo -i "/Volumes/Qing_Work/data/obs/Argo/WOD13_PFL_2000.nc" -lat 26 -lon -42 -r 2 -ds 20080601 -dr 30 -ot "tprof_test2.dat" -os "sprof_test2.dat" -dry


