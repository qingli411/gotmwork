;
; This script interpolates the WW3 output data
; on POP grid gx1v6
;
; Li Qing, 20170914

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin

; define parameters
    casen = "b1850_f19_gx1_vr12-ma_v2"

; define file names
    dataroot = "/Volumes/Qing_work/data/"
    diri = dataroot+"yellowstone/"+casen+"/"
    fili = "b1850_f19_gx1_vr12-ma_v2.ww3.daily.0101.nc"
    filo = "b1850_f19_gx1_vr12-ma_v2.ww3.daily.0101.gx16.nc"
    wgtFileName = dataroot+"/NCL_DATA/map_ww3a_TO_gx1v6_splice_150428.nc"

; input file
    f = addfile(diri+fili,"r")

; read grid file
    gdir    = dataroot+"/yellowstone/grid/"
    gfile   = "gx16.grids.nc"
    gf      = addfile(gdir+gfile,"r")
    lat2d   = gf->TLAT
    lon2d   = gf->TLONG
    ndim    = dimsizes(lat2d)
    nlat    = ndim(0)
    nlon    = ndim(1)
    ;nlat    = ispan(1,ndim(0),1)
    ;nlon    = ispan(1,ndim(1),1)
    ;printVarSummary(lat2d)

; create output data file
    system("/bin/rm -f "+diri+filo)
    out	= addfile(diri+filo,"c")

; add attributes
    fileAtt		= True
    fileAtt@title	= "WW3 output data on POP gx1v6 grid"
    fileAtt@src_data	= diri+fili
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "cesm_ww3a_to_gx1v6.ncl"
    fileAtt@author	= "Qing Li"
    fileattdef(out,fileAtt)

; options
    opt                = True
    opt@PrintTimings   = True

; loop over variables
;    vars = getfilevarnames(f)
    time = ispan(1,365,1)
    vars = (/"LaTurb","LaProj","LaSL","LaSLProj","HS","TMean","efactor"/)
    ;print(vars)
	nvar = dimsizes(vars)
	do iv = 0,nvar-1
        print(vars(iv))
        if (vars(iv) .eq. "efactor") then
            la = f->LaSLProj
            alpha = f->AlphaSL
            dat = abs(cos(alpha))*sqrt(1.0+(1.5*la)^(-2)+(5.4*la)^(-4))
            dat = where(dat .gt. 5.0, 5.0, dat)
            dat@longname = "Enhancement factor"
            dat@units = "unitless"
            dat@_FillValue = la@_FillValue
        else
	        ; read variable
	        var_in = vars(iv)
	        dat = f->$var_in$
        end if
	    ; remapping
	    datRegrid = ESMF_regrid_with_weights(dat, wgtFileName, opt)
        datRegrid!0 = "time"
        datRegrid&time = time
        datRegrid!1 = "nlat"
        ;datRegrid&nlat = nlat
        datRegrid!2 = "nlon"
        ;datRegrid&nlon = nlon
        delete(datRegrid@lat2d)
        delete(datRegrid@lon2d)
	    printVarSummary(datRegrid)
	    ; save variable
	    out->$vars(iv)$ = datRegrid
        delete(datRegrid)
	end do
    out->time = time
    out->TLAT = lat2d
    out->TLONG = lon2d

end
