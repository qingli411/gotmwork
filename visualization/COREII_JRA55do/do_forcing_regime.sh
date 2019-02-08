#!/bin/bash
# Compute and plot maps of forcing regimes for each month

source ../../set_tools.sh

ttag_list=("20080601-20080630" "20080701-20080731" "20080801-20080831" "20080901-20080930" "20081001-20081031" "20081101-20081130" "20081201-20081231" "20090101-20090131" "20090201-20090228" "20090301-20090331" "20090401-20090430" "20090501-20090531")
cname="JRA55-do_Global_dampV5d"
diag="BG12"
for ttag in ${ttag_list[@]}; do
    echo ${cname}_${ttag}_${diag}
    ./plot_map_forcing_regime -c ${cname} -t ${ttag} -d ${diag} -P -M
done
