#!/usr/bin/env bash

# setup paths and tools
source "../set_tools.sh"

${cmd_nc2dat_cdip} -i "${GOTMDATA_ROOT}/CDIP/Papa/166p1_rt.nc" -o "spec_file.dat" -ds 20160325 -de 20160411


