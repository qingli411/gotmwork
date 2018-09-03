#!/bin/sh

###################
#  GOTM vs CVMix  #
###################

# datestr="20130101-20140701"
# datestr="20130701-20140701"
# casename1="OCSPapa_KPP-CVMix"
# casename2="OCSPapa_KPP-Original"

#######################
#  COREII vs OCSPapa  #
#######################

# datestr="20090701-20091231"
# casename1="OCSPapa_KPP-CVMix"
# casename2="COREII_LAT50_LON215_KPP-CVMix"

######################
#  COREII vs OCSKEO  #
######################

# datestr="20080701-20081231"
datestr="20080102-20081231"
casename1="OCSKEO_KPP-CVMix"
casename2="COREII_LAT32_LON145_KPP-CVMix"

#######################################################################
#                              set paths                              #
#######################################################################

case1="${casename1}_${datestr}"
case2="${casename2}_${datestr}"
file1="gotm_out.nc"
file2="gotm_out.nc"

#######################################################################
#                            plot figures                             #
#######################################################################

echo "Testing plotts: dataset 1"
${cmd_plotts} -f ${GOTMRUN_ROOT}/${case1}/${file1} -v tx ty heat I_0 -o gotm_ts_surface_dset1.png

echo "Testing plotts: dataset 2"
${cmd_plotts} -f ${GOTMRUN_ROOT}/${case2}/${file2} -v tx ty heat I_0 -o gotm_ts_surface_dset2.png

echo "Testing plotts: comparison between two datasets"
${cmd_plotts} -f ${GOTMRUN_ROOT}/${case1}/${file1} -f2 ${GOTMRUN_ROOT}/${case2}/${file2} -v tx ty heat I_0 -o gotm_ts_surface_cmp.png

# for ptype in contourf pcolor scatter; do
    ptype="contourf"
    echo ${ptype}
echo "Testing plotpfl: dataset 1"
${cmd_plotpfl} -f ${GOTMRUN_ROOT}/${case1}/${file1} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_dset1.png

echo "Testing plotpfl: dataset 2"
${cmd_plotpfl} -f ${GOTMRUN_ROOT}/${case2}/${file2} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_dset2.png

echo "Testing plotpfl: obs in dataset 1"
${cmd_plotpfl} -f ${GOTMRUN_ROOT}/${case1}/${file1} -v temp_obs -ptype ${ptype} -o gotm_pfl_temp_${ptype}_obs.png

echo "Testing plotpfl: comparison between two datasets"
${cmd_plotpfl} -f ${GOTMRUN_ROOT}/${case1}/${file1} -f2 ${GOTMRUN_ROOT}/${case2}/${file2} -v temp -ptype ${ptype} -o gotm_pfl_temp_${ptype}_cmp.png

echo "Testing plotpfl: comparison between dataset 1 and obs"
${cmd_plotpfl} -f ${GOTMRUN_ROOT}/${case1}/${file1} -f2 ${GOTMRUN_ROOT}/${case1}/${file1} -v temp -v2 temp_obs -ptype ${ptype} -o gotm_pfl_temp_${ptype}_cmp_obs.png

# done
