#!/usr/bin/env bash

# setup paths and tools
source "../set_tools.sh"

# ${cmd_nc2dat_latlon} -i "${GOTMDATA_ROOT}/COREII_IAF/u_10.2009.23OCT2012.nc" -v "U_10_MOD" -o "u10_file.dat" -lat 25 -lon 180 -ds 20090101 -de 20090111 -maxd 300 -ignore_year

${cmd_nc2dat_latlon} -i "${GOTMDATA_ROOT}/GPCP/GPCP_1996_2015.nc" -v "precip" -o "precip_file.dat" -lat 25 -lon 180 -ds 20090101 -de 20090211



