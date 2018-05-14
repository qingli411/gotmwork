#!/bin/bash

../tools/plot_map -f COREII -t OSMOSIS -a mldMean -m deltaR -vmax 0 -vmin -200 -ds 20080615 -de 20091231 -ds_a 20090101 -de_a 20090331
../tools/plot_map -f COREII -t OSMOSIS -a mldMean -m maxNsqr -vmax 0 -vmin -200 -ds 20080615 -de 20091231 -ds_a 20090101 -de_a 20090331
