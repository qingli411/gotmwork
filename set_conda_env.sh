#!/bin/bash
# This script sets up conda environment for GOTM
#
# Qing Li, 20190707

# gotmwork environment file
gotmwork_env_file="${HOME}/.gotmwork_env.sh"
if [[ -f ${gotmwork_env_file} ]]; then
    source ${gotmwork_env_file}
else
    echo "** GOTMWORK environment not set. Use set_gotmwork_env.sh to set it up."
    exit 1
fi

# check if conda exists
conda_exe=$(which conda)
if [[ -z ${conda_exe} ]]; then
    echo -e "** Conda not found. Stop.\n"
    echo -e "** Follow instructions on"
    echo -e "   https://docs.conda.io/en/latest/miniconda.html"
    echo -e "   or"
    echo -e "   https://www.anaconda.com/distribution/"
    echo -e "   to install conda."
    exit 1
fi

# check os
os=$(uname)

# check if conda environment "gotm" exists.
# if not, create it.
conda_envs=$(conda env list | awk '{print $1}')
conda_envs=${conda_envs//\#/}
if [[ ${conda_envs} =~ "gotm" ]]; then
    echo -e "** Conda environment \"gotm\" exists. Skip.\n"
    exit 0
else
    echo -e "** Creating conda environment from ${GOTMWORK_ROOT}/gotm_env.yml..."
    conda env create -f ./gotm_env_${os}.yml
    # activate gotm environment
    if [[ $? == 0 ]]; then
        echo -e "** Activating conda environment \"gotm\"..."
        conda activate gotm
    fi
fi

