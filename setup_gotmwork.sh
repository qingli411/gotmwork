#!/bin/bash
#
# This script sets up environments for building GOTM source code,
# running cases, and postprocessing.
#
# Qing Li, 20180626

function python_version() {
    # Return 2 if in Python 2.x environment, 3 if in Python 3.x environment
    python -c "import sys; print(sys.version_info[0])"
}

function if_python_module() {
    # Return 1 if the module (name as argument) is installed in the current
    # Python environment, 0 otherwise
    python -c "import pkgutil; print(1 if pkgutil.find_loader(\"$1\") else 0)"
}

function if_python3_module() {
    # Return 1 if the module (name as argument) is installed in the current
    # Python environment, 0 otherwise
    python3 -c "import pkgutil; print(1 if pkgutil.find_loader(\"$1\") else 0)"
}

function inquire_dir() {
    local msg=$1
    local dft_val=$2
    echo -e "${msg} (${dft_val}):"
    read -r
    while [ ! -z "${REPLY}" ] && [ ! -d "${REPLY}" ]; do
        echo "Not a valid directory. Please try agian."
        read -r
    done
}

function get_inquire() {
    local dft_val=$1
    if [ ! -z "${REPLY}" ]; then
        val=${REPLY}
    else
        val=${dft_val}
    fi
    echo ${val}
}

function print_hline() {
    printf -- '-%.0s' {1..36}
    printf "\n"
}

# check if in python3 environment
if [ -x "$(command -v python3)" ]; then
    cmd_if_python_module="if_python3_module"
else
    cmd_if_python_module="if_python_module"
    if [[ $(python_version) == 2 ]]; then
        echo "Require Python 3.x environment. Stop."
        exit 1
    fi
fi

# check if required python modules are installed
python_module_list="argparse matplotlib mpl_toolkits.basemap netCDF4 numpy scipy subprocess xml.etree.ElementTree"
pm_counter=0
for pm in ${python_module_list}; do
    if [[ $(${cmd_if_python_module} ${pm}) == 0 ]]; then
        echo "Require Python module ${pm}. Not found."
        pm_counter=$((pm_counter+1))
    fi
done
if [[ ${pm_counter} > 0 ]]; then
    exit 1
fi

# gotmwork environment file
gotmwork_env_file="${HOME}/.gotmwork_env.sh"

if [ ! -f ${gotmwork_env_file} ]; then

    # set up the environment file
    print_hline
    echo "  Setting up GOTMWORK environment file"
    echo "  ${gotmwork_env_file}"
    print_hline

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

    # write to the environment file
    print_hline
    echo "  Write environment variables to file:"
    echo "  ${gotmwork_env_file}"
    print_hline
    echo ""

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
export CVMIX_ROOT=${cvmix_root}

GOTMWORK_ENV

    # print out the environment file for confirmation
    cat ${gotmwork_env_file}

    print_hline
    echo "  Done"
    print_hline
    echo ""
fi

if [[ $? == 0 ]]; then
    print_hline
    echo "Use 'source ${gotmwork_env_file}'"
    echo "to set up GOTMWORK environment"
    print_hline
fi

