#!/usr/bin/env bash

# setup paths and tools
source "../set_tools.sh"

${cmd_nc2dat_core2swr} -i "${GOTMDATA_ROOT}/COREII_IAF/ncar_rad.2008-2009.23OCT2012.nc" -o "swr_file.dat" -lat 25.7 -lon 180.2 -ds 20090325 -de 20090411



