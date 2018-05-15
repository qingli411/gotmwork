#!/bin/bash

turbmethod="KPP-CVMix"
../tools/plot_map -f COREII -t ${turbmethod} -a mldMean -m deltaR -vmax 0 -vmin -200 -ds 20080615 -de 20091231 -ds_a 20090101 -de_a 20090331 -dir_in /Users/qingli/scratch/VR5m_DT60s/${turbmethod}
