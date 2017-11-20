#!/bin/bash
# This script builds and installs GOTM code

blddir="${HOME}/Downloads/gotmbuild"
srcdir="${HOME}/models/gotm/code/src"

# create directory for build
mkdir -p ${blddir}
cd ${blddir}

# build
cmake ${srcdir} -DGOTM_USE_FABM=false

# install
make install
