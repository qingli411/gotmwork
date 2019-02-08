;
; This script extracts temperature and salinity profiles from the
; daily POP output
;
; Qing Li, 20180307

begin

; define parameters
    casen = "b1850_f19_gx1_vr12-ma_v2"

; define file names
    dataroot = "/Volumes/Qing_work/data/"
    diri = dataroot+"/yellowstone/"+casen+"/"
    filocn = casen+".pop.h.0101.nc"
    filo = casen+"_TSProfiles_monthly_0101.nc"

; input file
    focn = addfile(diri+filocn,"r")
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
    nlev = 39
    z_t = focn->z_t(0:nlev)*0.01 ; cm -> m
    salt = focn->SALT(:,0:nlev,:,:)
    printVarSummary(salt)
    temp = focn->TEMP(:,0:nlev,:,:)
    printVarSummary(temp)
    time_bound = focn->time_bound
    time = 0.5*(time_bound(:,0)+time_bound(:,1))-time_bound(0,0)-1

; create output data file
    system("/bin/rm -f "+diri+filo)
    setfileoption("nc","Format","LargeFile")
    out	= addfile(diri+filo,"c")

; missing value
    fmsg = default_fillvalue("float")
    dmsg = default_fillvalue("double")

; add attributes
    fileAtt		= True
    fileAtt@title	= "Monthly temperature and salinity profiles"
    fileAtt@src_data	= diri+filocn
    fileAtt@create_date	= systemfunc("date")
    fileAtt@src_code	= "cesm_prep_profiles.ncl"
    fileAtt@author	= "Qing Li"
    fileAtt@missing_value = fmsg
    fileattdef(out,fileAtt)

; save temperature and salinity
    temp = where(ismissing(temp), fmsg, temp)
    salt = where(ismissing(salt), fmsg, salt)
    temp@_FillValue = fmsg
    temp@missing_value = fmsg
    salt@_FillValue = fmsg
    salt@missing_value = fmsg
    out->TEMP = temp
    out->SALT = salt
    out->REGION_MASK = regionmask
    out->TAREA = tarea
    out->grav = grav
    out->cp_sw = cp_sw
    out->rho_sw = rho_sw
    out->rho_fw = rho_fw
    out->TLAT = lat2d
    out->TLONG = lon2d
    out->time = (/time/)
    out->time@units = "days since 0001-01-01-00000"
    out->z_t = (/z_t/)
    out->z_t@units = "m"
    out->z_t@valid_min = 5.0
    out->z_t@valid_max = 5375.0
end
