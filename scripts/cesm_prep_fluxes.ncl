;
; This script gathers all the necessary variables to diagnose
; ustar and wstar and La
; on POP grid gx1v6
;
; Qing Li, 20170914
;          20180307, used here to prepare flux data for GOTM

begin

; define parameters
    casen = "b1850_f19_gx1_vr12-ma_v2"

; define file names
    dataroot = "/Volumes/Qing_work/data/"
    diri = dataroot+"/yellowstone/"+casen+"/"
    filocn = casen+".pop.h.nday1.0101.nc"
    filwav = casen+".ww3.daily.0101.gx16.nc"
    filo = casen+"_parSpace_daily_0101.nc"

; input file
    focn = addfile(diri+filocn,"r")
    fwav = addfile(diri+filwav,"r")
    lat2d = focn->TLAT
    lon2d = focn->TLONG
    regionmask = focn->REGION_MASK
    tarea = focn->TAREA
    tarea = tarea/1e4 ; cm^2 -> m^2
    tarea@units = "m^2"
    grav = focn->grav
    grav = grav/1e2 ; cm/s^2 -> m/s^2
    grav@units = "m/s^2"
    cp_sw = focn->cp_sw
    cp_sw = cp_sw/1e4 ; erg/g/K -> J/kg/K
    cp_sw@units = "J/kg/K"
    rho_sw = focn->rho_sw
    rho_sw = rho_sw*1e3 ; g/cm^3 -> kg/m^3
    rho_sw@units = "kg/m^3"
    rho_fw = focn->rho_fw
    rho_fw = rho_fw*1e3 ; g/cm^3 -> kg/m^3
    rho_fw@units = "kg/m^3"
    ice = focn->IFRAC
    SSS = focn->SALT(:,0,:,:)
    printVarSummary(SSS)
    RHO = focn->RHO(:,0,:,:)
    printVarSummary(RHO)

; create output data file
    system("/bin/rm -f "+diri+filo)
    setfileoption("nc","Format","LargeFile")
    out	= addfile(diri+filo,"c")

; missing value
    fmsg = default_fillvalue("float")

; add attributes
    fileAtt		= True
    fileAtt@title	= "Data for diagnosing parameter space"
    fileAtt@src_data	= diri+filocn
    fileAtt@src_data2	= diri+filwav
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "cesm_prep_fluxes.ncl"
    fileAtt@author	= "Qing Li"
    fileAtt@missing_value = fmsg
    ; fileAtt@_FillValue = fmsg
    fileattdef(out,fileAtt)

; save salinity and density
    RHO = where(ismissing(RHO), fmsg, RHO*1e3) ; g/cm^3 -> kg/m^3
    SSS = where(ismissing(SSS), fmsg, SSS)
    RHO = where(ice .gt. 0.01, fmsg, RHO)
    SSS = where(ice .gt. 0.01, fmsg, SSS)
    RHO@units = "kg/m^3"
    RHO@_FillValue = fmsg
    RHO@missing_value = fmsg
    SSS@_FillValue = fmsg
    SSS@missing_value = fmsg
    out->RHO = RHO
    out->SSS = SSS
; loop over variables
;    vars = getfilevarnames(f)
    time = ispan(0,364,1)
    varswav = (/"LaTurb","LaProj","LaSL","LaSLProj","efactor"/)
    varsocn = (/"SST","SHF","SHF_QSW","SHF_HEAT","TAUX","TAUY","EVAP_F","PREC_F","MELT_F","ROFF_F","HMXL_DR","HBLT"/)
    unitsocn = (/"degC","W/m^2","W/m^2","W/m^2","N/m^2","N/m^2","kg/m^2/s","kg/m^2/s","kg/m^2/s","kg/m^2/s","m","m"/)
    ufactor = (/1,1,1,1,0.1,0.1,1,1,1,1,0.01,0.01/)
    ;print(vars)
	nvaro = dimsizes(varsocn)
	nvarw = dimsizes(varswav)
	do iv = 0,nvaro-1
        print(varsocn(iv))
        if (varsocn(iv) .eq. "SHF_HEAT") then
            dat = focn->SHF
            dat = (/focn->SHF - focn->SHF_QSW/)
            dat@long_name = "Total Surface Heat Flux, excluding SW"
        else
	        ; read variable
	        var_in = varsocn(iv)
	        dat = focn->$var_in$
        end if
        dat = where(ismissing(dat), fmsg, dat*ufactor(iv)) ; convert untis
        tmpdat = where(ice .gt. 0.01, fmsg, dat)
        dat = (/tmpdat/)
        dat@units = unitsocn(iv)
        dat@_FillValue = fmsg
        dat@missing_value = fmsg
	    ; save variable
	    out->$varsocn(iv)$ = dat
        delete(dat)
	end do
	do iv = 0,nvarw-1
        print(varswav(iv))
	    ; read variable
	    var_in = varswav(iv)
	    dat = fwav->$var_in$
        dat = where(ismissing(dat), fmsg, dat)
        tmpdat = where(ice .gt. 0.01, fmsg, dat)
        dat = (/tmpdat/)
        dat@_FillValue = fmsg
        dat@missing_value = fmsg
	    ; save variable
	    out->$varswav(iv)$ = dat
        delete(dat)
	end do
    out->REGION_MASK = regionmask
    out->TAREA = tarea
    out->grav = grav
    out->cp_sw = cp_sw
    out->rho_sw = rho_sw
    out->rho_fw = rho_fw
    out->TLAT = lat2d
    out->TLONG = lon2d
    out->time = (/time/)

end
