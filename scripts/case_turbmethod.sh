# script for setting up namelist depending on turbulence method

# set turbulence model
case ${turbmethod} in
    "K-EPSILON-SG")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 2
        ${cmd_nmlchange} -f gotmturb.nml -e tke_method -v 2
        ${cmd_nmlchange} -f gotmturb.nml -e len_scale_method -v 8
        ${cmd_nmlchange} -f gotmturb.nml -e stab_method -v 3
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "SMC")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e tke_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e len_scale_method -v 9
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_method  -v 1
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_coeff -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e length_lim -v .true.
        ${cmd_nmlchange} -f gotmturb.nml -e compute_c3 -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "SMCLT")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e tke_method -v 4
        ${cmd_nmlchange} -f gotmturb.nml -e len_scale_method -v 11
        ${cmd_nmlchange} -f gotmturb.nml -e e3 -v 5.0
        ${cmd_nmlchange} -f gotmturb.nml -e e6 -v 6.0
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_method  -v 4
        ${cmd_nmlchange} -f gotmturb.nml -e scnd_coeff -v 3
        ${cmd_nmlchange} -f gotmturb.nml -e length_lim -v .false.
        ${cmd_nmlchange} -f gotmturb.nml -e compute_c3 -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .true.
        ;;
    "OSMOSIS")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 98
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .true.
        ;;
    "KPP-GOTM")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 0
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "KPP-ROMS")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 2
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "KPP-CVMix")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 1
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 0
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "KPPLT-EFACTOR")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 1
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 1
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "KPPLT-ENTR")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 1
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 2
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "KPPLT-RWHGK")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 99
        ${cmd_nmlchange} -f kpp.nml -e kpp_opt -v 1
        ${cmd_nmlchange} -f kpp.nml -e langmuir_method -v 3
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .true.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .true.
        ;;
    "EPBL")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 100
        ${cmd_nmlchange} -f MOMturb.nml -e mode -v 0
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "JHL")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 100
        ${cmd_nmlchange} -f MOMturb.nml -e mode -v 2
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    "EPBL-LT")
        ${cmd_nmlchange} -f gotmturb.nml -e turb_method -v 100
        ${cmd_nmlchange} -f MOMturb.nml -e mode -v 0
        ${cmd_nmlchange} -f epbl.nml -e lt_enhance_form -v 3
        ${cmd_nmlchange} -f gotmmean.nml -e lagrangian_mixing -v .false.
        ${cmd_nmlchange} -f gotmmean.nml -e stokes_coriolis -v .false.
        ;;
    *)
        echo "Turbulence method ${turbmethod} not supported. Stop."
        exit 1
esac

