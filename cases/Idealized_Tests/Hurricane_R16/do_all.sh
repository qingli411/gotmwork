#!/bin/bash
# run all cases

# indices for turbulence models
# - start from 0
# - leave empty if run all models
imodel_list=( )

# list of turbulence models
model_list=( "KPP-CVMix" "KPP-ROMS" "KPPLT-EFACTOR" "KPPLT-ENTR" "KPPLT-RWHGK" "SMC" "SMCLT" "EPBL" "EPBL-LT" "OSMOSIS" "K-EPSILON-SG" "EPBL-RH18" "EPBL-RL19" "SMC-C01A" )

# run all the models if $imodel_list is empty
if [[ ${#imodel_list[@]} -eq 0 ]]; then
    imodel_list=( ${!model_list[@]} )
fi

# number of models
nmodel=${#imodel_list[@]}

# get pool of indices for each job
for ((m=0; m<nmodel; m++)); do
    jj=${imodel_list[m]}
    echo "Job $m: ./case_run -m ${model_list[jj]}"
    ./case_run -m ${model_list[jj]}
done
