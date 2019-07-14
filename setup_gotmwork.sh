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

# help information
print_hline
echo -e "** This script guides you through all the following steps to set up work"
echo -e "   environment for GOTM."
echo -e "     1. set up environment variables."
echo -e "     2. download source code of CVMix."
echo -e "     3. download source code of GOTM."
echo -e "     4. compile CVMix."
echo -e "     5. compile GOTM."
echo -e "     6. set up Python environment for GOTM via Conda."
echo -e "   Please follow the instructions and type in 'yes' or 'no' at each step.\n"
echo -e "** Note:"
echo -e "   A working fortran compiler (e.g., gfortran) and netCDF C and"
echo -e "   Fortran libraries are required for compiling CVMix and GOTM."
echo -e "   'nc-config' and 'nf-config' are used to determine the paths of"
echo -e "   netCDF libraries. Open another terminal and type in the following"
echo -e "   commands to check if they are installed:\n"
echo -e "   > which gfortran"
echo -e "   > which nc-config"
echo -e "   > which nf-config\n"

# a note for Mac user
os=$(uname)
if [[ ${os} == "Darwin" ]]; then
    echo -e "** A note for Mac users:"
    echo -e "   If netCDF libraries are installed using Homebrew,"
    echo -e "   'nf-config' may not work correctly and give an error"
    echo -e "   'nf-config not yet implemented for cmake builds'."
    echo -e "   MacPort version seems fine. A workaround is to use your"
    echo -e "   own 'nf-config' instead of the Homebrew version. An"
    echo -e "   example (NetCDF v4.5.0 installed in '/usr/local' with"
    echo -e "   Fortran compiler 'gfortran') is given in"
    echo -e "   https://github.com/qingli411/gotm/blob/master/scripts/nf-config\n"
fi

# continue
inquire_yes_no "Continue?" "yes"
run_setup=$(get_inquire "yes")
if [[ ${run_setup} == "no" ]]; then
    echo -e "** Stop setting up work environment for GOTM."
    exit 0
fi

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
    check_error $? "build_cvmix.sh"
else
    echo -e "** Skip compiling CVMix. You may use ./build_cvmix.sh to compile CVMix."
    echo -e "   Note that CVMix needs to be compiled prior to compiling GOTM."
fi

# build GOTM
print_hline
inquire_yes_no "Compile GOTM?" "no"
bld_gotm=$(get_inquire "no")
if [[ ${bld_gotm} == "yes" ]]; then
    ./build_gotm.sh -clean -build
    check_error $? "build_gotm.sh"
else
    echo -e "** Skip compiling GOTM. You may use ./build_gotm.sh to compile GOTM."
fi

# check Python environment
print_hline
echo -e "** Python 3 is needed by many tools in ./tools/ to preprocess input"
echo -e "   data for GOTM, and by the scripts and Jupyter notebooks in"
echo -e "   ./visualization/ to postprocess and visualize the results. Some"
echo -e "   tools and scripts may have dependencies on specific versions of"
echo -e "   the required Python packages. Therefore, it is recommended to create"
echo -e "   a conda environment for GOTM from ./gotm_env.yml to ensure compatibility."
echo -e "   However, it may not be necessary if you already have these packages"
echo -e "   installed. Check ./gotm_env.yml for required Python packages.\n"
inquire_yes_no "Proceed to create a conda environment namded 'gotm' and activate it?" "no"
set_conda=$(get_inquire "no")
if [[ ${set_conda} == "yes" ]]; then
    ./set_conda_env.sh
    check_error $? "set_conda_env.sh"
else
    echo -e "** Skip setting conda environment. You may use"
    echo -e "   ${GOTMWORK_ROOT}/set_conda_env.sh"
    echo -e "   to set up conda environment."
fi

# done
print_hline
echo -e "** Done\n"

