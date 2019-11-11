#!/bin/bash
# Do the diagnostics (ts) for each month

source ../../set_tools.sh

# flag for updating data
update_data=True

# number of parallel jobs
njob=2

# case name
cname="COREII_Global_dampV5d_3h"

# indices for time tags
# - start from 0
# - leave empty if run all time tags
isetup_list=(  )

# indices for diagnostics
# - start from 0
# - leave empty if run all diagnostics
imodel_list=(  )

# list of time tags
setup_list=( \
            "VR1m_DT600s_20080601-20080630" \
            "VR1m_DT600s_20080701-20080731" \
            "VR1m_DT600s_20080801-20080831" \
            "VR1m_DT600s_20080901-20080930" \
            "VR1m_DT600s_20081001-20081031" \
            "VR1m_DT600s_20081101-20081130" \
            "VR1m_DT600s_20081201-20081231" \
            "VR1m_DT600s_20090101-20090131" \
            "VR1m_DT600s_20090201-20090228" \
            "VR1m_DT600s_20090301-20090331" \
            "VR1m_DT600s_20090401-20090430" \
            "VR1m_DT600s_20090501-20090531" \
          )

# list of turbulence models
model_list=( "KPP-CVMix" "KPP-ROMS" "KPPLT-EFACTOR" "KPPLT-ENTR" "KPPLT-RWHGK" "SMC" "SMCLT" "EPBL" "EPBL-LT" "OSMOSIS" "K-EPSILON-SG" "EPBL-RH18" "EPBL-RL19" "SMC-C01A" )

# run all the dates if $isetup_list is empty
if [[ ${#isetup_list[@]} -eq 0 ]]; then
    isetup_list=( ${!setup_list[@]} )
fi

# number of dates
nsetup=${#isetup_list[@]}

# run all the models if $idiag_list is empty
if [[ ${#imodel_list[@]} -eq 0 ]]; then
    imodel_list=( ${!model_list[@]} )
fi

# number of models
nmodel=${#imodel_list[@]}

# number of total runs
nrun=$((nsetup*nmodel))
# number of cases in a job
if [[ $((nrun % njob)) -eq 0 ]]; then
    njcase=$((nrun/njob))
else
    njcase=$((nrun/njob+1))
fi
# flags
if [[ ${update_data} == "True" ]]; then
    s1_flag="-U"
else
    s1_flag=""
fi
# print a summary
echo "--------"
echo "Number of parallel jobs: ${njob}"
echo "Number of total runs: ${nrun}"
echo "Maximum number of runs in each job: ${njcase}"
echo "--------"
echo "List of runs:"
# get pool of indices for each job
j_list=()
k_list=()
j=0
k=0
for ((m=0; m<njob; m++)); do
    for ((n=0; n<njcase; n++)); do
        if [[ ${k} -lt ${nsetup} ]]; then
            j_list+=(${imodel_list[j]})
            k_list+=(${isetup_list[k]})
            # print out the command of runs
            jj=${imodel_list[j]}
            kk=${isetup_list[k]}
            echo "./diagnostics_mapts -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag}"
            echo "./diagnostics_mapts_langmuir -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag}"
            echo "./diagnostics_mapts_airsea -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag}"
            j=$((j+1))
            if [[ ${j} -eq ${nmodel} ]]; then
                k=$((k+1))
                j=0
            fi
        fi
    done
done
echo "--------"
# submit jobs
for ((m=0; m<njob; m++)); do
    {
    for ((n=0; n<njcase; n++)); do
        ii=$((m*njcase+n))
        if [[ ${ii} -lt $((nsetup*nmodel)) ]]; then
            jj=${j_list[ii]}
            kk=${k_list[ii]}
            ./diagnostics_mapts -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag} > log.${ii};
            ./diagnostics_mapts_langmuir -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag} > log.${ii};
            ./diagnostics_mapts_airsea -c ${cname} -s ${setup_list[kk]} -m ${model_list[jj]} ${s1_flag} > log.${ii};
        fi
    done
    } &
done
