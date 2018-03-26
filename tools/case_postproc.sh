# script for post-processing

# variable list
pp_vars="temp temp_obs salt salt_obs"
# plot surface fluxes
${cmd_plotts} -f ${outname}.nc -v tx ty heat I_0 -o "gotm_ts_surface_forcing.pdf"
# plot temperature and salinity
for pp_var in ${pp_vars}; do
    ${cmd_plotpfl} -f ${outname}.nc -v ${pp_var} -o "gotm_pfl_${pp_var}.pdf" -ptype contourf
done
