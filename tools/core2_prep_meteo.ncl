;
; This scripts merge meteorology data from CORE-II
;   Variable list:
;       u10 - Surface wind at 10 meters (m/s), x-component
;       v10 - Surface wind at 10 meters (m/s), y-component
;       slp - Sea level pressure (hPa)
;       sat - Surface air temperature (degK or degC)
;       hum - Specific huminidy (kg/kg) or relative huminity (%)
;       cld - Cloud cover
;
; Qing Li, 20180310

begin

; define parameters
    yyyy = "2008-2009"
    dlabel = "23OCT2012"

; define file names
    gotmdata_root = getenv("GOTMDATA_ROOT")
    diri = gotmdata_root+"/COREII_IAF/"
    filo = "meteo."+yyyy+"."+dlabel+".nc"
    filelist = (/"u_10", "v_10", "t_10", "q_10", "slp"/)
    varlist = (/"U_10_MOD", "V_10_MOD", "T_10_MOD", "Q_10_MOD", "SLP"/)

; missing value
    fmsg = default_fillvalue("float")
    dmsg = default_fillvalue("double")

; create output data file
    system("/bin/rm -f "+diri+filo)
    setfileoption("nc","Format","NetCDF4")
    out	= addfile(diri+filo,"c")

; input file
    nv = dimsizes(filelist)
    src_data = ""
    do i=0,nv-1
        fili = filelist(i)+"."+yyyy+"."+dlabel+".nc"
        src_data = fili+"\n"+src_data
        f = addfile(diri+fili, "r")
        if (i .eq. 0) then
            lat = f->LAT
            lon = f->LON
            time = f->TIME
            time_bnds = f->TIME_bnds
            out->LAT = lat
            out->LON = lon
            out->TIME = time
            out->TIME_bnds = time_bnds
        end if
        var = f->$varlist(i)$
        var = where(ismissing(var), fmsg, var)
        var@missing_value = fmsg
        var@_FillValue = fmsg
        if (varlist(i) .eq. "SLP" )
            var = (/var/)*0.01
            var@units = "hPa"
        end if
        out->$varlist(i)$ = var
    end do

    ; assume zero cloud cover
    cld = (/var/)*0.0
    cld!0 = "TIME"
    cld&TIME = time
    cld!1 = "LAT"
    cld&LAT = lat
    cld!2 = "LON"
    cld&LON = lon
    cld@long_name = "Cloud cover (all zero)"
    cld@units = "unitless"
    cld@missing_value = fmsg
    cld@_FillValue = fmsg
    out->CLD = cld

; add attributes
    fileAtt		= True
    fileAtt@title	= "Meteorology data from CORE-II"
    fileAtt@src_data	= src_data
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "core2_prep_meteo.ncl"
    fileAtt@author	= "Qing Li"
    fileAtt@missing_value = fmsg
    fileattdef(out,fileAtt)

end
