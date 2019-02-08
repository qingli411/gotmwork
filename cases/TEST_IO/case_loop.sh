# This is the main loop over turbulent methods, relaxation time for
# velocity damping to suppress inertial oscillation.
# Used by: OCSPapa, OSMOSIS, COREII
#
# Qing Li, 20180718

#######################################################################
#                           Loop over cases                           #
#######################################################################

l_test="no"

if [ ${l_test} == "yes" ]; then
    # name of the turbulence model
    turblist=(SMCLT)
    # relaxation time
    vrelaxlist=(none)
else
    # name of the turbulence model
    turblist=(KPP-CVMix KPP-ROMS KPPLT-EFACTOR KPPLT-ENTR KPPLT-RWHGK OSMOSIS EPBL EPBL-LT SMC SMCLT K-EPSILON-SG)
    # relaxation time
    vrelaxlist=(1d 10d)
fi

# output file name
outname="gotm_out"

# 1 m vertical resolution
grid_method=0
ddu=0
ddl=0
let nlev=depth
dt=600

# loop over turbulent methods
for turbmethod in ${turblist[@]}; do

# loop over time step
for vrelax in ${vrelaxlist[@]}; do

    case ${vrelax} in
        "1d")
            trelax=86400
            ;;
        "10d")
            trelax=864000
            ;;
        *)
            echo "Relaxation time ${vrelax} not supported. Stop."
            exit 1
    esac

    # case name
    casename="TEST_IO/${title}/${turbmethod}_dampV_${vrelax}"
    echo ${casename}

    # set output frequency (3-hourly output)
    let nsave=10800/dt

    # create run directory
    rundir="${GOTMRUN_ROOT}/${casename}"
    mkdir -p ${rundir}
    cd ${rundir}

    # copy base case
    cp ${basecase}/*.dat ./
    cp ${basecase}/*.nml ./

    if [ -f ${basecase}/*.gz ]; then
        cp ${basecase}/*.gz ./
        # decompress input data
        for f in *.gz; do
            gunzip -f ${f}
        done
    fi

    # set run parameters
    if [ -n "${datestart}" ] && [ -n "${dateend}" ]; then
        start_time="${datestart:0:4}-${datestart:4:2}-${datestart:6:2} 00:00:00"
        stop_time="${dateend:0:4}-${dateend:4:2}-${dateend:6:2} 00:00:00"
        ${cmd_nmlchange} -f gotmrun.nml -e start -v "${start_time}"
        ${cmd_nmlchange} -f gotmrun.nml -e stop -v "${stop_time}"
    fi
    if [ -n "${sprof_file}" ]; then
        ${cmd_nmlchange} -f obs.nml -e s_prof_file -v ${sprof_file}
    fi
    if [ -n "${tprof_file}" ]; then
        ${cmd_nmlchange} -f obs.nml -e t_prof_file -v ${tprof_file}
    fi
    ${cmd_nmlchange} -f gotmrun.nml -e title -v ${title}
    ${cmd_nmlchange} -f gotmrun.nml -e out_fn -v ${outname}
    ${cmd_nmlchange} -f gotmrun.nml -e dt -v ${dt}
    ${cmd_nmlchange} -f gotmrun.nml -e nsave -v ${nsave}
    ${cmd_nmlchange} -f gotmrun.nml -e nlev -v ${nlev}
    ${cmd_nmlchange} -f gotmrun.nml -e eq_state_method -v 4
    ${cmd_nmlchange} -f gotmmean.nml -e grid_method -v ${grid_method}
    ${cmd_nmlchange} -f gotmmean.nml -e ddu -v ${ddu}
    ${cmd_nmlchange} -f gotmmean.nml -e ddl -v ${ddl}

    ${cmd_nmlchange} -f obs.nml -e vel_relax_tau -v ${trelax}

    # set turbulence method
    source ${scpt_case_turbmethod}

    # run
    ${cmd_gotm} 2> log.${outname}

    # plot some figures
    if [ ${l_test} == "yes" ]; then
        source ${curdir}/case_postproc.sh
    fi

    # clean up input data
    if [ $? == 0 ]; then
        rm -f *.dat
    fi

done
done
