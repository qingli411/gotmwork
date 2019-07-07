#!/bin/bash
#
# This script sets up working environment for GOTM
#
# Qing Li, 20190706

function print_hline() {
    printf "\n"
    printf -- '-%.0s' {1..36}
    printf "\n"
}

function check_error() {
    local error=$1
    local proc=$2
    if [[ ${error} != 0 ]]; then
        echo -e "Error in ${proc}. Stop.\n"
        exit 1
    fi
}

function inquire_yes_no() {
    local msg=$1
    local dft_val=$2
    echo "${msg} (${dft_val}):"
    read -r
    while [[ ! -z ${REPLY} ]] && [[ ${REPLY} != "yes" ]] && [[ ${REPLY} != "no" ]] ; do
        echo -e "Not a valid answer. Please type \"yes\" or \"no\"."
        read -r
    done
}

function get_inquire() {
    local dft_val=$1
    if [[ ! -z ${REPLY} ]]; then
        val=${REPLY}
    else
        val=${dft_val}
    fi
    echo ${val}
}


# set up gotmwork environment variables
print_hline
./scripts/set_gotmwork_env.sh
check_error $? "set_gotmwork_env.sh"

# source environment variables
source ~/.gotmwork_env.sh

# get CVMix
print_hline
inquire_yes_no "Download CVMix from Github?" "yes"
get_cvmix=$(get_inquire "yes")
if [[ ${get_cvmix} == "yes" ]]; then
    ./scripts/get_cvmix.sh
    check_error $? "get_cvmix.sh"
fi

# get GOTM
print_hline
inquire_yes_no "Download GOTM from Github?" "yes"
get_gotm=$(get_inquire "yes")
if [[ ${get_gotm} == "yes" ]]; then
    ./scripts/get_gotm.sh
    check_error $? "get_gotm.sh"
fi

# build CVMix
print_hline
inquire_yes_no "Compile CVMix?" "no"
bld_cvmix=$(get_inquire "no")
if [[ ${bld_cvmix} == "yes" ]]; then
    ./build_cvmix.sh -clean -build
    check_error $? "build_cvmix"
fi

# build GOTM
print_hline
inquire_yes_no "Compile GOTM?" "no"
bld_gotm=$(get_inquire "no")
if [[ ${bld_gotm} == "yes" ]]; then
    ./build_gotm.sh -clean -build
    check_error $? "build_gotm"
fi

# check Python environment

# done
print_hline
echo -e "** Done\n"

