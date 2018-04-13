#!/usr/bin/env bash

# ../tools/nc2dat_latlon -i "/Volumes/Qing_Work/data/yellowstone/b1850_f19_gx1_vr12-ma_v2/b1850_f19_gx1_vr12-ma_v2_STProfiles_monthly_0101.nc" -v "TEMP" -o "temp_file.dat" -lat 25 -lon 180 -ds 20160101 -de 20160111 -maxd 300

# ../tools/nc2dat_latlon -i "/Volumes/Qing_Work/data/yellowstone/b1850_f19_gx1_vr12-ma_v2/b1850_f19_gx1_vr12-ma_v2_parSpace_daily_0101.nc" -v "SHF_HEAT" -o "cesm_file.dat" -lat 25 -lon 180 -ds 20160101 -de 20160111 -maxd 300 -ignore_year

# ../tools/nc2dat_latlon -i "/Volumes/Qing_Work/data/COREII_IAF/u_10.2009.23OCT2012.nc" -v "U_10_MOD" -o "u10_file.dat" -lat 25 -lon 180 -ds 20090101 -de 20090111 -maxd 300 -ignore_year

../tools/nc2dat_latlon -i "/Volumes/Qing_Work/data/obs/GPCP_1996_2015.nc" -v "precip" -o "precip_file.dat" -lat 25 -lon 180 -ds 20090101 -de 20090211



