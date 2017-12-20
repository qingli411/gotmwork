#!/bin/sh

echo "Testing plot_obs"
../tools/obs_show -xml ../data/OCSPapa.xml -root ${HOME}/models/gotm/gotmwork \
    -data ${HOME}/data/OCS/Papa/2016 -out ${HOME}/models/gotm/gotmwork/test/OCSPapa \
    -ds 20160325 -de 20160411
