;
; This script gathers all the necessary variables to diagnose
; ustar and wstar and La
; on POP grid gx1v6
;
; Li Qing, 20170914

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin

; define parameters
    casen = "b1850_f19_gx1_vr12-ma_v2"

; define file names
    diri = "/glade/u/home/qingli/pwork/data/"+casen+"/"
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
    grav = focn->grav
    cp_sw = focn->cp_sw
    rho_sw = focn->rho_sw
    rho_fw = focn->rho_fw
    ice = focn->IFRAC
    SSS = focn->SALT(:,0,:,:)
    printVarSummary(SSS)
    RHO = focn->RHO(:,0,:,:)
    printVarSummary(RHO)

; create output data file
    system("/bin/rm -f "+diri+filo)
    setfileoption("nc","Format","LargeFile")
    out	= addfile(diri+filo,"c")

; add attributes
    fileAtt		= True
    fileAtt@title	= "Data for diagnosing parameter space"
    fileAtt@src_data	= diri+filocn
    fileAtt@src_data2	= diri+filwav
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "gatherData.ncl"
    fileAtt@author	= "Qing Li"
    fileattdef(out,fileAtt)

; missing value
    fmsg = default_fillvalue("float")

; save salinity and density
    RHO = where(ice .gt. 0.05, fmsg, RHO)
    SSS = where(ice .gt. 0.05, fmsg, SSS)
    out->RHO = RHO
    out->SSS = SSS
; loop over variables
;    vars = getfilevarnames(f)
    time = ispan(1,365,1)
    varswav = (/"LaTurb","LaProj","LaSL","LaSLProj"/)
    varsocn = (/"SST","SHF","SHF_QSW","TAUX","TAUY","EVAP_F","PREC_F","MELT_F","ROFF_F","HMXL_DR","HBLT"/) 
    ;print(vars)
	nvaro = dimsizes(varsocn)
	nvarw = dimsizes(varswav)
	do iv = 0,nvaro-1
        print(varsocn(iv))
	    ; read variable
	    var_in = varsocn(iv)
	    dat = focn->$var_in$
        tmpdat = where(ice .gt. 0.05, fmsg, dat)
        dat = (/tmpdat/)
	    ; save variable
	    out->$varsocn(iv)$ = dat
        delete(dat)
	end do
	do iv = 0,nvarw-1
        print(varswav(iv))
	    ; read variable
	    var_in = varswav(iv)
	    dat = fwav->$var_in$
        tmpdat = where(ice .gt. 0.05, fmsg, dat)
        dat = (/tmpdat/)
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
    ;out->time = time
    out->TLAT = lat2d
    out->TLONG = lon2d

end 
