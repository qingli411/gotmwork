# script for post-processing

# extract variables to reduce the output file size
if [[ -f gotm_out.nc ]]; then
    # varlist="time,lat,lon,zi,z,u_taus,u10,v10,airt,airp,hum,precip,evap,I_0,qe,qh,qb,heat,tx,ty,sst,sss,u0_stokes,v0_stokes,delta,La_Turb,La_SL,La_SLP1,La_SLP2,theta_WW,theta_WL,Ekin,Epot,Eturb,temp,salt,rho,u,v,NN,avh,num,nuh,gamh,rad"
    # ncks -O -v $varlist gotm_out.nc gotm_out_s1.nc
    varlist="temp_obs,salt_obs,u_obs,v_obs,sst_obs,h,ga,bioshade,NNT,NNS,xP,drag,fric,gamu,gamv"
    ncks -O -x -v $varlist gotm_out.nc gotm_out_s1.nc
fi

if [[ $? == 0 ]]; then
    rm -f gotm_out.nc
fi
