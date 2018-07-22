#!/bin/sh

# setup paths and tools
source "../set_tools.sh"

echo "Testing plot_ts with GOTM data"
${cmd_plotts} -f ocspapa.nc -v tx ty heat I_0 -o gotm_ts_surface.pdf -ds 20160320 -de 20160411

echo "Testing plot_pfl with GOTM data"
${cmd_plotpfl} -f ocspapa.nc -v temp -ptype pcolor -ds 20160330 -de 20160405 -o gtom_pfl_temp.pdf

echo "Testing plot_ts with OCSPapa data"
${cmd_plotts} -f wind_10meter_50n145w_hr.cdf -v UZS_2422 VZS_2423 -o papa_ts_wind.pdf -ds 20160201 -de 20160205

echo "Testing plot_pfl with OCSPapa data"
${cmd_plotpfl} -f t50n145w_hr.cdf -v T_20 -ptype scatter -ds 20160201 -de 20160205 -o papa_pfl_temp.pdf

echo "Testing plot_ts with SPURS data"
${cmd_plotts} -f SPURS_1_D_F1H.nc -v QS QL QB QH -o spurs_ts_flux.pdf -ds 20121201 -de 20121203

echo "Testing plot_pfl with SPURS data"
${cmd_plotpfl} -f SPURS_1_D_TS.nc -v TEMP -ptype scatter -ds 20121201 -de 20121203 -o spurs_pfl_temp.pdf
