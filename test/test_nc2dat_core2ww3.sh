#!/usr/bin/env bash

# setup paths and tools
source "../set_tools.sh"

${cmd_nc2dat_core2ww3} -i "${GOTMDATA_ROOT}/WAVEWATCH_COREII/ww3.2008-2009_usp.nc" -o "usp_file.dat" -lat 25.7 -lon 180.2 -ds 20090325 -de 20090411 -usp

${cmd_nc2dat_core2ww3} -i "${GOTMDATA_ROOT}/WAVEWATCH_COREII/ww3.2008-2009.nc" -o "wave_file.dat" -lat 25.7 -lon 180.2 -ds 20090325 -de 20090411 -v hs fp dir

${cmd_nc2dat_core2ww3} -i "${GOTMDATA_ROOT}/WAVEWATCH_COREII/ww3.2008-2009.nc" -o "wave_file.dat" -lat 25.7 -lon 180.2 -ds 20090325 -de 20090411



