# This script set up paths and tools used by case_run

# paths
tooldir="${workdir}/tools"
nmldir="${workdir}/namelist"
xmldir="${workdir}/data"

# tools
cmd_nmlchange="${tooldir}/nmlchange"
cmd_case_preproc="${tooldir}/case_preproc"
cmd_nc2dat="${tooldir}/nc2dat"
cmd_nc2dat_argo="${tooldir}/nc2dat_argo"
cmd_nc2dat_cdip="${tooldir}/nc2dat_cdip_spec"
cmd_nc2dat_core2ww3="${tooldir}/nc2dat_core2ww3"
cmd_nc2dat_latlon="${tooldir}/nc2dat_latlon"
cmd_plotpfl="${tooldir}/plotpfl"
cmd_plotts="${tooldir}/plotts"
cmd_is_sea="${tooldir}/is_sea"

# bash script
scpt_case_postproc="${tooldir}/case_postproc.sh"
