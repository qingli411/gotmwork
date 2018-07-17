#!/usr/bin/env bash

# setup paths and tools
source "../set_tools.sh"

date_str="20080115"
lat1=-10
lon1=2
lat2=10
lon2=2
rdeg=2
rdate=20
echo "Find Argo profiles near [lon=${lon1}, lat=${lat1}] on ${date_str}"

echo "Test 1a: dry run"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat1} -lon ${lon1} -r ${rdeg} -d ${date_str} -dr ${rdate} -dry

echo "Test 1b: save output"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat1} -lon ${lon1} -r ${rdeg} -d ${date_str} -dr ${rdate} -ot "tprof_test1.dat" -os "sprof_test1.dat"

echo "Test 2a: ignore year, dry run"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat1} -lon ${lon1} -r ${rdeg} -d ${date_str} -dr ${rdate} -dry -iy

echo "Test 2b: ignore year, save output"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat1} -lon ${lon1} -r ${rdeg} -d ${date_str} -dr ${rdate} -ot "tprof_test2.dat" -os "sprof_test2.dat" -iy

echo ""
echo "Find Argo profiles near [lon=${lon2}, lat=${lat2}] on ${date_str}"

echo "Test 3a: ignore year, no profile available, dry run"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat2} -lon ${lon2} -r ${rdeg} -d ${date_str} -dr ${rdate} -dry -iy

echo "Test 3b: ignore year, no profile available, save output"
../tools/nc2dat_argo -i "${GOTMDATA_ROOT}/Argo/WOD13_PFL_2000.nc" -lat ${lat2} -lon ${lon2} -r ${rdeg} -d ${date_str} -dr ${rdate} -ot "tprof_test3.dat" -os "sprof_test3.dat" -iy

