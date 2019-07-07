#!/bin/bash
#
# This script sets up environments for building GOTM source code,
# running cases, and postprocessing.
#
# Qing Li, 20180626

function get_inquire() {
    local dft_val=$1
    if [[ ! -z ${REPLY} ]]; then
        val=${REPLY}
    else
        val=${dft_val}
    fi
    echo ${val}
}

function make_dir() {
    local dft_val=$1
    dir=$(get_inquire ${dft_val})
    if [[ ${dir} != /* ]]; then
        echo "1"
    else
        mkdir -p ${dir}
        echo $?
    fi
}

function inquire_dir() {
    local msg=$1
    local dft_val=$2
    echo "${msg} (${dft_val}):"
    read -r
    while [[ $(make_dir ${dft_val}) != 0 ]]; do
        echo -e "Not a valid full path of a directory. Please try agian."
        read -r
    done
}

# gotmwork environment file
gotmwork_env_file="${HOME}/.gotmwork_env.sh"

if [[ ! -f ${gotmwork_env_file} ]]; then

    # set up the environment file
    echo -e "** Setting up GOTMWORK environment file"
    echo -e "   ${gotmwork_env_file}\n"

    # gotmwork root directory
    gotmwork_root=$(pwd)

    # default directories
    dft_gotm_root=$(dirname ${gotmwork_root})
    dft_cvmix_root="${dft_gotm_root}/CVMix-src"
    dft_gotmcode_root="${dft_gotm_root}/code"
    dft_gotmexe_root="${dft_gotm_root}/exe"
    dft_gotmbuild_root="${dft_gotm_root}/build"
    dft_gotmrun_root="${dft_gotm_root}/run"
    dft_gotmdata_root="${dft_gotm_root}/data"
    dft_gotmfig_root="${dft_gotm_root}/fig"
    dft_gotmarchive_root="${dft_gotm_root}/archive"

    # instruction
    echo -e "Type in the full path of the directories. Leave it"
    echo -e "empty to use the default values in parentheses.\n"

    # inquire environment variables
    inquire_dir "Root directory of CVMix source code" ${dft_cvmix_root}
    cvmix_root=$(get_inquire ${dft_cvmix_root})
    inquire_dir "Root directory of GOTM source code" ${dft_gotmcode_root}
    gotmcode_root=$(get_inquire ${dft_gotmcode_root})
    inquire_dir "Root directory of GOTM input data" ${dft_gotmdata_root}
    gotmdata_root=$(get_inquire ${dft_gotmdata_root})
    inquire_dir "Directory to build GOTM" ${dft_gotmbuild_root}
    gotmbuild_root=$(get_inquire ${dft_gotmbuild_root})
    inquire_dir "Directory of GOTM executable" ${dft_gotmexe_root}
    gotmexe_root=$(get_inquire ${dft_gotmexe_root})
    inquire_dir "Directory to run GOTM" ${dft_gotmrun_root}
    gotmrun_root=$(get_inquire ${dft_gotmrun_root})
    inquire_dir "Directory for visualizations of the results" ${dft_gotmfig_root}
    gotmfig_root=$(get_inquire ${dft_gotmfig_root})
    inquire_dir "Directory to archive GOTM output data" ${dft_gotmarchive_root}
    gotmarchive_root=$(get_inquire ${dft_gotmarchive_root})

    # write to the environment file
    echo -e "\n** Write environment variables to file:"
    echo -e "   ${gotmwork_env_file}\n"

    cat > ${gotmwork_env_file} << GOTMWORK_ENV
# Environment variables for gotmwork
#
export GOTMWORK_ROOT=${gotmwork_root}
export GOTMCODE_ROOT=${gotmcode_root}
export GOTMDATA_ROOT=${gotmdata_root}
export GOTMBUILD_ROOT=${gotmbuild_root}
export GOTMEXE_ROOT=${gotmexe_root}
export GOTMRUN_ROOT=${gotmrun_root}
export GOTMFIG_ROOT=${gotmfig_root}
export GOTMARCHIVE_ROOT=${gotmarchive_root}
export CVMIX_ROOT=${cvmix_root}

GOTMWORK_ENV

else
    echo -e "** Environment variables already set in:"
    echo -e "   ${gotmwork_env_file}"
    echo -e "   Remove it first to set up new environment.\n"
fi

# print out the environment file for confirmation
echo -e "** Current environment variables:\n"
cat ${gotmwork_env_file}
