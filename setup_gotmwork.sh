#!/bin/bash
#
# This script sets up working environment for GOTM
#
# Qing Li, 20190706

function print_hline() {
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

# set up gotmwork environment variables
print_hline
./scripts/set_gotmwork_env.sh
check_error $? "set_gotmwork_env.sh"

# source environment variables
source ~/.gotmwork_env.sh

# acquire CVMix
print_hline
./scripts/get_cvmix.sh
check_error $? "get_cvmix.sh"

# acquire GOTM
print_hline
./scripts/get_gotm.sh
check_error $? "get_gotm.sh"

# build CVMix


# build GOTM

# check Python environment

# done
print_hline
echo -e "** Done\n"

