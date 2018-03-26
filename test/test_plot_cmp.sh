#!/bin/sh

datadir="/Users/qingli/work/gotmrun"
case1="OCSPapa_KPP-CVMix_20130701-20140701"
case2="OCSPapa_KPP-Original_20130701-20140701"
# case2="COREII_LAT50_LON215_KPP-CVMix_20090701-20091231"
file1="gotm_out.nc"
file2="gotm_out.nc"

echo "Testing plotts: dataset 1"
../tools/plotts -f ${datadir}/${case1}/${file1} -v tx ty heat I_0 -o gotm_ts_surface_dset1.pdf

echo "Testing plotts: dataset 2"
../tools/plotts -f ${datadir}/${case2}/${file2} -v tx ty heat I_0 -o gotm_ts_surface_dset2.pdf

echo "Testing plotts: comparison between two datasets"
../tools/plotts -f ${datadir}/${case1}/${file1} -f2 ${datadir}/${case2}/${file2} -v tx ty heat I_0 -o gotm_ts_surface_cmp.pdf

# for ptype in contourf pcolor scatter; do
    ptype="contourf"
    echo ${ptype}
echo "Testing plotpfl: dataset 1"
../tools/plotpfl -f ${datadir}/${case1}/${file1} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_dset1.pdf

echo "Testing plotpfl: dataset 2"
../tools/plotpfl -f ${datadir}/${case2}/${file2} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_dset2.pdf

echo "Testing plotpfl: comparison between two datasets"
../tools/plotpfl -f ${datadir}/${case1}/${file1} -f2 ${datadir}/${case2}/${file2} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_cmp.pdf
# done
