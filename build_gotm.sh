#!/bin/bash
# This script builds and installs GOTM code
#
# Qing Li, 20190706

function usage() {
    echo -e "\nBuild GOTM source code\n"
    echo -e "Usage:\n$0 [arguments]\n"
    echo -e "Build the source code without cleaning the old build"
    echo -e "if no arguments are used\n"
    echo -e "Arguments:"
    echo -e "  -build: build source code"
    echo -e "  -clean: clean old build\n"
}

# set paths
# gotmwork environment file
gotmwork_env_file="${HOME}/.gotmwork_env.sh"
if [[ -f ${gotmwork_env_file} ]]; then
    source ${gotmwork_env_file}
else
    echo "** GOTMWORK environment not set. Use set_gotmwork_env.sh to set it up."
    exit 1
fi

blddir="${GOTMBUILD_ROOT}"
srcdir="${GOTMCODE_ROOT}/src"
exedir="${GOTMEXE_ROOT}"
cvmixdir="${CVMIX_ROOT}"

# flags to build GOTM with CVMix
l_cvmix=true

# flags to build and clean build
l_build=false
l_clean=false

# check input arguments
if [[ $# == 0 ]]; then
    # build the source code by defaut
    l_build=true
else
    for param in "$@"; do
        shift
        case $param in
            -build)
                l_build=true
                ;;
            -clean)
                l_clean=true
                ;;
            *)
                usage
                exit 1
                ;;
        esac
    done
fi

if [[ "${l_clean}" == "true" ]]; then
    echo "** Cleanning old GOTM build..."
    # clean the build directory
    rm -rf ${blddir}
fi

if [[ "${l_build}" == "true" ]]; then
    echo "** Building GOTM..."
    # create directory for build
    mkdir -p ${blddir}
    cd ${blddir}
    # create directory for exe
    mkdir -p ${exedir}

    # build
    if [[ "${l_cvmix}" == "true" ]]; then
        cmake ${srcdir} -DGOTM_USE_FABM=false \
            -DCMAKE_INSTALL_PREFIX:PATH=${exedir} \
            -DGOTM_USE_CVMix=true
    else
        cmake ${srcdir} -DGOTM_USE_FABM=false \
            -DCMAKE_INSTALL_PREFIX:PATH=${exedir} \
            -DGOTM_USE_CVMix=false
    fi

    # install
    make install
fi

