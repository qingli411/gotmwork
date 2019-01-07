# script for post-processing

# variable list
pp_vars="temp temp_obs salt salt_obs"
# plot surface fluxes
${cmd_plotts} -f ${outname}.nc -v tx ty heat I_0 precip u0_stokes v0_stokes -o "gotm_ts_surface_forcing.png"
# plot temperature and salinity
for pp_var in ${pp_vars}; do
    ${cmd_plotpfl} -f ${outname}.nc -v ${pp_var} -o "gotm_pfl_${pp_var}.png" -ptype contourf
done
