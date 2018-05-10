# script for setting up namelist depending on turbulence method

# set turbulence model
case ${turbmethod} in
    "OSMOSIS")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 98
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .true.
        ;;
    "SMC")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e tke_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e len_scale_method -v 9
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_coeff -v 3
        ;;
    "SMCLT")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e tke_method -v 4
        ${cmd_nmlchange} -f gotmturb.nml -e len_scale_method -v 11
        ${cmd_nmlchange} -f gotmturb.nml -e e6 -v 6.0
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_coeff -v 3
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .true.
        ;;
    "KPP-CVMix")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e lcvmix -v .true.
        ;;
    "KPP-GOTM")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e lcvmix -v .false.
        ;;
    "KPPLT-EFACTOR")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e lcvmix -v .true.
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 1
        ;;
    "KPPLT-ENTR")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e lcvmix -v .true.
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 2
        ;;
    "EPBL")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 100
	    ${cmd_nmlchange} -f MOMturb.nml -e Mode -v 0
        ;;
    *)
        echo "Turbulence method ${turbmethod} not supported. Stop."
        exit 1
esac

