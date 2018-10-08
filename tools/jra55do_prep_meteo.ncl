;
; This scripts merge meteorology data from JRA55-do
;   Variable list:
;       u10 - Surface wind at 10 meters (m/s), x-component
;       v10 - Surface wind at 10 meters (m/s), y-component
;       slp - Sea level pressure (hPa)
;       sat - Surface air temperature (degK or degC)
;       hum - Specific huminidy (kg/kg) or relative huminity (%)
;       cld - Cloud cover
;
; Qing Li, 20180831

begin

; define parameters
    yyyy = "2008-2009"
    dlabel = "18Oct2017"

; define file names
    gotmdata_root = getenv("GOTMDATA_ROOT")
    diri = gotmdata_root+"/JRA55-do/"
    filo = "meteo."+yyyy+"."+dlabel+".nc"
; TODO: add sea level pressure <31-08-18, Qing Li> ;
    filelist = (/"u_10", "v_10", "t_10", "q_10", "slp"/)
    varlist = (/"uas_10m", "vas_10m", "tas_10m", "huss_10m", "psl"/)

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
            lat = f->latitude
            lon = f->longitude
            time = f->time
            time_bnds = f->time_bnds
            out->latitude = lat
            out->longitude = lon
            out->time = time
            out->time_bnds = time_bnds
        end if
        var = f->$varlist(i)$
        var = where(ismissing(var), fmsg, var)
        var@missing_value = fmsg
        var@_FillValue = fmsg
        if (varlist(i) .eq. "psl" )
            var = (/var/)*0.01
            var@units = "hPa"
        end if
        out->$varlist(i)$ = var
    end do

    ; assume zero cloud cover
    cld = (/var/)*0.0
    cld!0 = "time"
    cld&time = time
    cld!1 = "latitude"
    cld&latitude = lat
    cld!2 = "longitude"
    cld&longitude = lon
    cld@long_name = "Cloud cover (all zero)"
    cld@units = "unitless"
    cld@missing_value = fmsg
    cld@_FillValue = fmsg
    out->cld = cld

; add attributes
    fileAtt		= True
    fileAtt@title	= "Meteorology data from JRA55-do"
    fileAtt@src_data	= src_data
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "jra55do_prep_meteo.ncl"
    fileAtt@author	= "Qing Li"
    fileAtt@missing_value = fmsg
    fileattdef(out,fileAtt)

end
