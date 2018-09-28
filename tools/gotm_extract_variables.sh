#!/bin/bash

varlist="time,lat,lon,zi,z,u_taus,u10,v10,airt,airp,hum,precip,evap,I_0,qe,qh,qb,heat,tx,ty,sst,sss,u0_stokes,v0_stokes,delta,KPP_OSBL,La_Turb,La_SL,La_SLP1,La_SLP2,theta_WW,theta_WL,Ekin,Epot,Eturb,temp,salt,rho,u,v,NN,avh,num,nuh,gamh,rad"

ncks -O -v $varlist gotm_out.nc gotm_out_s1.nc
