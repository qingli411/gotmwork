#!/bin/bash
# This script builds CVMix code
# Python2 is required to build CVMix. Starting from v0.95-beta, Python3 is supported.
#
# Qing Li, 20190706

function usage() {
    echo -e "\nBuild CVMix source code\n"
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

srcdir="${CVMIX_ROOT}/src"

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

# clean build
if [[ "${l_clean}" == "true" ]]; then
    echo "** Cleanning old CVMix build..."
    cd ${srcdir}
    make allclean
    cd - &> /dev/null
fi

# build
if [[ "${l_build}" == "true" ]]; then
    echo "** Building CVMix..."
    cd ${srcdir}
    make
    cd - &> /dev/null

fi

