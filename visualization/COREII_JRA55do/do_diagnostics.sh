#!/bin/bash
# Do the diagnostics for each month

source ../../set_tools.sh

# flag for updating data
update_data=False

# number of parallel jobs
njob=4

# case name
cname="JRA55-do_Global_dampV5d"

# indices for time tags
# - start from 0
# - leave empty if run all time tags
ittag_list=( )

# indices for diagnostics
# - start from 0
# - leave empty if run all diagnostics
idiag_list=( )

# list of time tags
ttag_list=("20080601-20080630" "20080701-20080731" "20080801-20080831" "20080901-20080930" "20081001-20081031" "20081101-20081130" "20081201-20081231" "20090101-20090131" "20090201-20090228" "20090301-20090331" "20090401-20090430" "20090501-20090531")

# list of diagnostics
diag_list=( "mld_deltaR_mean" "mld_deltaRp1_mean" "SST_mean" "PE_delta" "Nsqr_mld_mean" )

# run all the dates if $ittag_list is empty
if [[ ${#ittag_list[@]} -eq 0 ]]; then
    ittag_list=( ${!ttag_list[@]} )
fi

# number of dates
nttag=${#ittag_list[@]}

# run all the models if $idiag_list is empty
if [[ ${#idiag_list[@]} -eq 0 ]]; then
    idiag_list=( ${!diag_list[@]} )
fi

# number of models
ndiag=${#idiag_list[@]}

# number of total runs
nrun=$((nttag*ndiag))
# number of cases in a job
if [[ $((nrun % njob)) -eq 0 ]]; then
    njcase=$((nrun/njob))
else
    njcase=$((nrun/njob+1))
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
        if [[ ${k} -lt ${nttag} ]]; then
            j_list+=(${idiag_list[j]})
            k_list+=(${ittag_list[k]})
            # print out the command of runs
            jj=${idiag_list[j]}
            kk=${ittag_list[k]}
            echo "./plot_map_diagnostics -c ${cname} -t ${ttag_list[kk]} -d ${diag_list[jj]}"
            j=$((j+1))
            if [[ ${j} -eq ${ndiag} ]]; then
                k=$((k+1))
                j=0
            fi
        fi
    done
done
echo "--------"
# # submit jobs
if [[ ${update_data} == "True" ]]; then
    s1_flag="-U"
else
    s1_flag=""
fi
for ((m=0; m<njob; m++)); do
    {
    for ((n=0; n<njcase; n++)); do
        ii=$((m*njcase+n))
        if [[ ${ii} -lt $((nttag*ndiag)) ]]; then
            jj=${j_list[ii]}
            kk=${k_list[ii]}
            # step 1
            ./plot_map_diagnostics -c ${cname} -t ${ttag_list[kk]} -d ${diag_list[jj]} ${s1_flag} > logs1.${ii};
            # step 2
            ./gen_map_mask -c ${cname} -t ${ttag_list[kk]} -d ${diag_list[jj]} > logs2.${ii};
            # step 3
            ./plot_map_diagnostics -c ${cname} -t ${ttag_list[kk]} -d ${diag_list[jj]} -M -P > logs3.${ii};
        fi
    done
    } &
done
