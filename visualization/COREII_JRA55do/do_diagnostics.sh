#!/bin/bash
# Do the diagnostics for each month

source ../../set_tools.sh

ttag_list=("20080601-20080630" "20080701-20080731" "20080801-20080831" "20080901-20080930" "20081001-20081031" "20081101-20081130" "20081201-20081231" "20090101-20090131" "20090201-20090228" "20090301-20090331" "20090401-20090430" "20090501-20090531")
# ttag_list=( "20090501-20090531")
# diag_list=( "mld_deltaR_mean" "SST_mean" "PE_delta" "Nsqr_mld_mean" )
diag_list=( "mld_deltaR_mean" )
cname="JRA55-do_Global_dampV5d"
for ttag in ${ttag_list[@]}; do
    for diag in ${diag_list[@]}; do
        echo ${cname}_${ttag}_${diag}
        # step 1
        # ./plot_map_diagnostics -c ${cname} -t ${ttag} -d ${diag} -U -P
        # step 2
        # ./gen_map_mask -c ${cname} -t ${ttag} -d ${diag}
        # step 3
        ./plot_map_diagnostics -c ${cname} -t ${ttag} -d ${diag} -M -P
    done
done
