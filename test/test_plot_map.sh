#!/bin/bash

# setup paths and tools
source "../set_tools.sh"

turbmethod="KPPLT-ENTR"
${cmd_plot_map} -f COREII -t ${turbmethod} -a mldMean -m deltaR -vmax 0 -vmin -500 -ds 20080801 -de 20080831 -ds_a 20080801 -de_a 20080831 -dir_in ${GOTMRUN_ROOT}/COREII_GLOBAL/VR1m_DT600s_20080801-20080831/${turbmethod}
