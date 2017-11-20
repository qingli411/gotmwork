
      PROGRAM MYMAIN
C
C  This main program is used to drive a Second Moment Closure ocean model
C  through a diurnal cycle using an idealized, through prescribed
C  daily cycle of heating, or, an arbitrary (observed) time
C  series of surface flux data.
C
C  The program assumes the following logical units:
C   5 is the terminal for input,
C   6 is the terminal for output,
C   8 may be opened to read air/sea flux data (optional),
C   9 may be opened to read the initial T, S profile (optional),
C  10 may be opened to store matrix data (enabled),
C  11 is opened to store lists from each time step,
C  13 is opened to store profile data if requested.
C
C
C  The original code was written by Patrice Klein. The main program was
C  heavily revised by Jim Price in March, 1991 who claimed at that time to be  responsible for
C  any and all errors.  However, that code was appropriated and used by Ramsey Harcourt in 2010-12
C  as the basis for coding a modification to the Kantha & Clayson (2004)
C  second moment closure model for Langmuir turbulence. This predecessor of Jim's was chosen because of it's particularly
C  simple coding and attendant maleability for model development. (Thanks Jim!) The code was extensively altered, and Jim should
C  not be bothered about any surviving errors, though the ancient equation of state still needs an
C  update and the smoothing option seems to be non-conserving and hence nearly off.  The code version
c  transferred is the same as the one used to generate 'E6=5' SMC model results in
C  Harcourt, R.R. A second moment closure model of Langmuir turbulence, Journal of Physical Oceanography, 2013.
C  and further modified for near-surface pressure-strain & pressure-scalar correlations in Harcourt, R.R.
C  An Improved second moment closure model of Langmuir turbulence, Journal of Physical Oceanography, 2015.
C  This version has been passed initially to the Fox-Kemper group for use in collaborative work, and
C  if you did not receive this copy via that group with a report to R.R. Harcourt of the transfer, you would be wise
C  to check for any updates, as this model is a work in progress. Don't be shy.
C  Ramsey R Harcourt
C  harcourt@uw.edu
C  (+1) 206 221 4662
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KAPPA,KMS,KM,KH,KQ2,KQ2L,L
      REAL*8 KN,KBG,KHMIN,KMMIN,KBGCON
      CHARACTER*1 ANY
C     CHARACTER*16 IFILE
      CHARACTER*30 IFILE
C
      PARAMETER (NSOBS = 34999)
      PARAMETER (NWOBS = 34999)
      PARAMETER (NWFRQ = 199)
      PARAMETER (NZOBS = 399)
C
C  The arrays below are used to store the initial Z,T,S.
C
      DIMENSION ZO(NZOBS), TO(NZOBS), SO(NZOBS), DO(NZOBS) ,
     +          UO(NZOBS), VO(NZOBS), Q2O(NZOBS),Q2LO(NZOBS)
C
C  The arrays below store times series of surface flux data.
C
      DIMENSION DAYA(NSOBS), QIA(NSOBS), QLA(NSOBS), TXA(NSOBS),
     X TYA(NSOBS), EMPA(NSOBS)
C  The arrays below store times series of surface wave data.
C
      DIMENSION DAYW(NWOBS), FRQ(NWFRQ), H2SPC(NWOBS,NWFRQ),
     +          STKSSPC(NWOBS,NWFRQ), STKDSPC(NWOBS,NWFRQ),
     +          COSTKDR(NWOBS,NWFRQ), SNSTKDR(NWOBS,NWFRQ)
C
C
C  Where arrays have had to be dimensioned here or in subroutines,
C  they are set to 399 so that they can be easily found.
C
      PARAMETER (KB=250)
      DIMENSION ZZP(KB),DZP(KB),ZP(KB)
      DIMENSION ATTSCL(NWFRQ),ATT(KB,NWFRQ)
      COMMON/BLKT/TF(KB),T(KB),TB(KB),SF(KB),S(KB),SB(KB)
      COMMON/BLKU/UF(KB),U(KB),UB(KB),VF(KB),V(KB),VB(KB)
      COMMON/BLKS/STKX(KB),STKY(KB),STKX0,STKY0
      COMMON/BLKQ/Q2F(KB),Q2(KB),Q2B(KB),WW(KB)
      COMMON/BLKQL/Q2LF(KB),Q2L(KB),Q2LB(KB)
      COMMON/BLKR/RHO(KB),WT(KB)
      COMMON/BLKWU/WUSURF,WVSURF,USTR2,WUBOT,WVBOT
      COMMON/BLKWT/WTSURF,WSSURF,PROFIL(KB)
      COMMON/BLKZ/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),UMOL
      COMMON/BLKH/H,DT
      COMMON/BLKA/A(KB),C(KB),VH(KB),VHP(KB)
      COMMON/BLKK/KMS(KB),KM(KB),KH(KB),KQ2(KB),KQ2L(KB),L(KB),KBG(KB)
      COMMON/INDS/KBM1,KBM2
      COMMON/MODEL/GM(KB),GH(KB),SS(KB),SM(KB),SH(KB),KN(KB),FZS(KB),
     1   SPROD(KB),BPROD(KB),PROD(KB),DTEF(KB),VPROD(KB),GV(KB),GS(KB)
      COMMON/CSTUR/A1,B1,A2,B2,C1,C2,C3,E1,E2,VONKAR,SEF,QMIN,
     +          SQ2,SQ2L,COF1,COF2,COF3,COF4,COF5,ZRS,ZRSW,
     +          WAVEFAC,WAGE,ZRB,GHMIN,GHMAX,KHMIN,KMMIN,
     +          CKBG,KBGCON,E3,E4,E5,E6,C1S,C2S,
     +          DS0,DS1,DS2,DM0,DM1,DM2,DM3,DH0,DH1,DH2,DH3

      DIMENSION UU(KB),VV(KB)
C
c      PARAMETER (NT1=24)
c
c      DIMENSION UX(NT1,KB),VX(NT1,KB),TTX(NT1,KB)
c      DIMENSION DX(NT1,KB),EZ1(KB)
C
      DIMENSION TOX(35000),TOY(35000),RA(35000),HSURF(35000)
      DIMENSION TSEX(35000)
c
c
      DATA KAPPA/0.4D0/,EPS/1.E-10/,RPI/3.14159D0/
C
C  NZMAX should match the size of the initial profile arrays
C  dimensioned above, and NMET should match the flux arrays.
C
      NZMAX = 399
      NMET = 34999
      NWAV = 34999
C
      ISTART = 1
      IEND = 119999
C
C
        OPEN(10,FILE='MYMAT.DAT',STATUS='UNKNOWN',
C    +      FORM='UNFORMATTED')
     +      FORM='FORMATTED')
        OPEN(11,FILE='MYLIS.DAT',STATUS='UNKNOWN',
     +      FORM='FORMATTED')
        OPEN(13,FILE='MYPROF.DAT',STATUS='UNKNOWN',
     +      FORM='FORMATTED')

      BILEC1=0.D0
      BILEC2=0.D0
      BILECT=0.D0
      ECB=0.D0
C
c
      GRAV=9.806D0
C...not very conservative smoothing!      SMOTH=0.1D0
      SMOTH=0.01D0
C     UMOL=1.34E-5
      UMOL=1.83E-6
      UMOL=0.D0
C
      A1 = 0.92D0
      B1 = 16.6D0
      A2 = 0.74D0
      B2 = 10.1D0
      C1 = 0.08D0
C...  WALT version
C     C1S = C1
      C1S = 0.D0
C...  WALT version
      C2 = 0.7D0
C...  WALT version
C     C2S = C2

C..Test!!
      C2S = 0.D0
C     C2S = D0
C..Test!!
C...  WALT version
      C3 = 0.2D0
      E1 = 1.8D0
      E3 = 5.D0
      E2 = 1.D0
      E4 = 1.33D0
      E5 = 0.04D0
C.....Vortex production in q2l eqn is E6 l (Vprod)
C     E6 = 9.D0
C     E6 = 7.2D0
C     E6 = 5.D0
C     E6 = 4.D0
C     E6 = 7.D0
C...WALT6?
C     E6 = 3.6D0
C...WALT9?
      E6 = 6.D0
C...WALT6?

      VONKAR = 0.4D0
      SEF = 1.D0
      QMIN = 1.E-8
      GHMAX = 0.029D0
      GHMIN = -0.28D0
      SQ2 = 0.41D0
      SQ2L = 0.41D0
      ZRS = 0.1D0
C...for very young hurricane waves, else change to developed wage=1.2
C     WAGE = 0.6D0
      WAGE = 0.D0
C     WAVEFAC = 100.D0
      WAVEFAC = 0.D0
      ZRB = 0.1D0
C0103 KHMIN = 5.E-5
C0103 KMMIN = 1.E-4
C0103 KC04 uses smaller minima:
C123  KHMIN = 1.E-5
C123  KMMIN = 1.E-5
C0103 KC04 uses smaller minima:
C0104 This minimum is for background GM-driven mixing, not appropriate for comp w/LES
C123  KMMIN = 0.D0
      KHMIN = 0.D0
      KMMIN = 0.D0
C0104
      CKBG = 5.E-3
      KBGCON = 5.E-3

      COF1=A2*(1.D0-6.D0*A1/B1)
      COF2=3.D0*A2*B2*(1.D0-C3)+18.D0*A1*A2
      COF3=A1*(1.D0-3.D0*C1-6.D0*A1/B1)
      COF4=18.D0*A1*A1+9.D0*A1*A2*(1.D0-C2)
      COF5=9.D0*A1*A2

      DS0=A1*(1.D0-3.D0*C1S-6.D0*A1/B1)
      DM0=A1*(1.D0-3.D0*C1-6.D0*A1/B1)
      DH0=A2*(1.D0-6.D0*A1/B1)
C

      WRITE (6,800)
  800 FORMAT (
     X 1X, 'This program provides a means to run the Mellor',/,
     X 1X, 'turbulence closure model. To run the model you must',/,
     X 1X, 'prescribe air/sea fluxes, the initial',/,
     X 1X, 'temperature and salinity profile, and several',/,
     X 1X, 'paprameters defining the site (i.e., latitude).',/,
     X 1X, 'Some of these (Ug,Vg) are relevent to the bottom',/,
     X 1X, 'boundary layer which this model includes, and are',/,
     X 1X, 'different from the PWP model. So watch out!',/,
     X 1X, 'The site and model parameters are as follows',/,
     X 1X, '(default values and units given in parens.):',/)
      WRITE (6,802)
  802 FORMAT (
     X 1X, 'RLAT,  latitude (31.)',/,
     X 1X, 'BETA1, longwave extinction coefficient (0.6 m)',/,
     X 1X, 'BETA2, shortwave extinction coefficient (20. m)',/,
     X 1X, 'Ug, Vg, the imposed geostrophic current (0., 0. m/sec)',/,

     X 1X,/,
     X 1X, 'DT,  time step (400. sec)',/,
     X 1X, 'DAYS, the number of days to run (1. day)',/,
     X 1X, 'H, the water column depth (80. m)',/,
     X 1X,/)
C
C  Initialize RLAT, the latitude.
C
      RLAT = 31.D0
C
C  Initialize the solar radiation absorption;
C  BETA1 is the  non-penetrating part of the insolation,
C  BETA2 is the penetrating (shortwave) part.
C  Values below are appropriate for fairly clear waters.
C
      BETA1 = 0.6D0
      BETA2 = 20.D0
C
      UG = 0.D0
      VG = 0.D0
C
      WRITE (6,888)
  888 FORMAT (1X,'ENTER RLAT, BETA1, BETA2, Ug, Vg')
      READ (5,*) RLAT, BETA1, BETA2, UG, VG
C
C
      XM=.62D0
      ETA1=1.D0/BETA1
      ETA2=1.D0/BETA2
C
      COR = 2.D0*7.292E-5*SIN(RLAT*3.141D0/180.D0)
C
      PGX = COR*VG
      PGY = -COR*UG
C
C  Define the temporal and vertical resolution. In this model, the vertical
c  resolution (dz) is given by the water column depth, H, divided by the number
c  of grid cells, KB.
c
      DT1 = 400.D0
      DAYS = 1.D0
      H = 80.D0
C
      WRITE (6,706)
  706 FORMAT (1X,/,1X,'ENTER DT, DAYS, H')
      READ (5,*) DT1, DAYS, H
C
      DTD = DT1/8.64E4
C
C
      KBM1=KB-1
      KBM2=KB-2
      CALL DEPTH(Z,ZZ,DZ,DZZ,KB,KBM1)
      TIME=0.D0
      IDT1=INT(DT1)
      DT2=2.D0*DT1
C
C  Initialize the Z, T, S, U, V profiles.
C
C  Read in the data used to define the initial T,S profiles
C  from logical unit NUNIT (read in below). The data
C  is presumed to be free format, with each record being a
C  triple (Z, T, S).  The first record should be the surface
C  value, i.e., 0., 20., 36. is a reasonable surface value,
C  and the last true value should be deep enough to avoid
C  entering into the caluclation (unless intentional).
C  The end of file is defined either by the end
C  of the data (if read from a disk file), or by a z value less
C  than -10. These data are then linearly interpolated to give
C  profiles at DZ resolution needed by the model.
C
      WRITE (6,458)
  458 FORMAT (1X,/,
     X 1X,'The temperature and salinity  profile can',/,
     X 1X,'be defaulted to a profile that is hardwired into the',/,
     X 1X,'program (ENTER 1, the default), or enter Z, T, S data'/,
     X 1X,'from the terminal (ENTER 2), or read the data from a',/,
     X 1X,'disk file (ENTER 3).')
C
      MZTS = 1
      READ (5,*) MZTS
C
C
      IF(MZTS.EQ.1) THEN
      NTS = 3
      ZO(1) = 0.D0
      ZO(2) = 20.D0
      ZO(3) = 200.D0
      TO(1) = 20.D0
      TO(2) = 20.D0
      TO(3) = 5.D0
      SO(1) = 36.D0
      SO(2) = 36.D0
      SO(3) = 36.D0
      UO(1) = 0.D0
      UO(2) = 0.D0
      UO(3) = 0.D0
      VO(1) = 0.D0
      VO(2) = 0.D0
      VO(3) = 0.D0
      Q2O(1)=1.E-8
      Q2O(2)=1.E-8
      Q2O(3)=1.E-8
      Q2LO(1)=0.D0
      Q2LO(2)=1.E-8
      Q2LO(3)=1.E-8
      END IF
      IF(MZTS.EQ.1) GO TO 53
C
      IF(MZTS.EQ.3) THEN
      WRITE (6,3290)
 3290 FORMAT (1X,/,1X,'ENTER THE NAME OF THE Z,T,S DATA FILE')
      READ (5,116) IFILE
      OPEN (UNIT=9,FILE=IFILE)
      END IF
C
C
      NTS = 0
      DO 50 J=1,NZMAX
      IF(MZTS.EQ.2) THEN
      WRITE (6,876)
  876 FORMAT (1X,'ENTER A Z,T,S,U,V,Q2,Q2L (END = Z < -10)')
      READ (5,*) ZO(J), TO(J), SO(J), UO(J), VO(J), Q2O(J), Q2LO(J)
      END IF
C
      IF(MZTS.EQ.3) THEN
      READ(9,*,END=53) ZO(J),TO(J),SO(J),UO(J),VO(J),Q2O(J),Q2LO(J)
      END IF
C
      IF(ZO(J).LT.-9.) GO TO 53
      NTS = NTS + 1
C     IF(MZTS.EQ.3) WRITE (6,503) ZO(J),TO(J),SO(J),UO(J),VO(J),Q2O(J),Q2LO(J)
  503 FORMAT (1X,'INITIAL Z, T, S, U, V, Q2, Q2L ARE', 7F10.2)
   50 CONTINUE
   53 CONTINUE
C
C  Interpolate the initial data to get a profile at DZ resolution.
C  The arrays ZZP and DZP are first computed to find the grid levels.
C
      DO 52 J=1,KB
      ZZP(J) = -1.D0*H*ZZ(J)
      ZP(J) = -1.D0*H*Z(J)
C     DZP(J) = -DZZ(J)*H
C...Looks wrong. We want a dimensional dz that is centered at H*ZZ, so that is H*diff(Z), or +DZ(J)*H -RRH
C...DZP is used later below; ZP remains as above as depth in mid-layer for interpolation purposes.
      DZP(J) = DZ(J)*H
C...This changes the sign of heat and PE diagnostics below
      CALL INTX(ZO,TO,NZOBS,NTS,ZZP(J),TB(J))
      CALL INTX(ZO,SO,NZOBS,NTS,ZZP(J),SB(J))
      CALL INTX(ZO,UO,NZOBS,NTS,ZZP(J),UB(J))
      CALL INTX(ZO,VO,NZOBS,NTS,ZZP(J),VB(J))
      CALL INTX(ZO,Q2O,NZOBS,NTS,ZP(J),Q2B(J))
      CALL INTX(ZO,Q2LO,NZOBS,NTS,ZP(J),Q2LB(J))
c      WRITE (6,54) J,ZZP(J),ZZ(J),TB(J),SB(J),RHO(J)
      WRITE (6,54) J,ZZP(J),ZZ(J),TB(J),SB(J),Q2B(J),Q2LB(J)
   54 FORMAT (1X,'DENSITY PROF.',I4,2X,8F10.3)
   52 CONTINUE
C

      DO 5 K=1,KB
      UB(K)=UB(K)+UG
      U(K)=UB(K)
      UF(K)=U(K)
      VB(K)=VB(K)+VG
      V(K)=VB(K)
      VF(K)=V(K)
      WT(K)=0.D0
      T(K)=TB(K)
      TF(K)=T(K)
      S(K)=SB(K)
      SF(K)=S(K)
C
C....WARNING: various different initial conditions crash/don't or are relevant/not...
      Q2B(K)=1.E-8
C     Q2B(K)=0.D0
      Q2(K)=Q2B(K)
      Q2F(K)=Q2(K)
      Q2LB(K)=1.E-9
C     Q2LB(K)=0.D0
      Q2L(K)=Q2LB(K)
      Q2LF(K)=Q2L(K)
      L(K)=Q2LB(K)/(Q2B(K)+SMALL)
      DS0=A1*(1.D0-3.D0*C1S-6.D0*A1/B1)
      KMS(K)=DS0*L(K)*SQRT(Q2(K))
      DM0=A1*(1.D0-3.D0*C1-6.D0*A1/B1)
      KM(K)=DM0*L(K)*SQRT(Q2(K))
      DH0=A2*(1.D0-6.D0*A1/B1)
      KH(K)=DH0*L(K)*SQRT(Q2(K))
C     KMS(K)=0.D0
C     KM(K)=0.D0
C     KH(K)=0.D0
C     KMS(K)=KMMIN
C     KM(K)=KMMIN
C     KH(K)=KHMIN
      KQ2(K)=SQ2*KH(K)
      KQ2L(K)=SQ2L*KH(K)
C...WALT version
      FZS(K)=0.D0
      WW(K)=0.D0
C...WALT version
    5 CONTINUE
C
      CALL DENS(T,S,RHO)
C
C
C  Initialization of Z, T, S, U, V, Q2, Q2L is finished.
C
      DVT = 0.D0
      VVT = 10.D0
C
C  Set up to input air/sea flux data.
C
      WRITE (6,457)
  457 FORMAT (1X,//,
     X 1X,'You can prescribe the air/sea flux data by',/,
     X 1X,'specifying the parameters of an idealized diurnal',/,
     X 1X,'cycle (ENTER 1 (default)), or, read the air/sea',/,
     X 1X,'flux data from a disk file (ENTER 2).')
C
      METU = 1
      READ (5,*) METU
C
      IF(METU.EQ.1) THEN
C
      WRITE (6,2790)
 2790 FORMAT (1X,/,1X,'In order to prescribe an idealized',/,
     X 1X, 'diurnal cycle, the program will now ',/,
     X 1X, 'query for the values of the following variables, all ',/,
     X 1X, 'of which may be defaulted to values in parens.',/,
     X 1X, 'The default values are chosen to match day 131',/,
     X 1X, 'of the PWP JGR paper (91) C7, 8411-8427, 1986.',/)
      WRITE (6,801)
  801 FORMAT (
     X 1X, 'QIMAX, noon amplitude of insolation (978. W/m**2)',/,
     X 1X, 'QL, steady heat loss (-126. W/m**2)',/,
     X 1X, 'TX, steady east wind stress (0.07 Pa)',/,
     X 1X, 'EMP, evaporation minus precipitation (0. m/sec)',/,
     X 1X, 'PQFAC,  duration of insolation (10. hours)',/)
C
C  Initialize the amplitude of solar insolation (QI), the
C  heat loss (QL), and the wind stress (TX, uses east only).
C
      QIMAX = 978.D0
      QL = -126.D0
      TX = 0.07D0
      EMP = 0.D0
C
C  Initialize PQFAC, the duration of daylight in hours.
C
      PQFAC = 10.D0
C
C
      WRITE (6,777)
  777 FORMAT (1X,'ENTER QIMAX, QL, TX, EMP, PQFAC')
      READ (5,*) QIMAX, QL, TX, EMP, PQFAC
C
C  Convert PQFAC to a non-d form.
C
      IF(PQFAC.LT.0.001) PQFAC = 0.001D0
      PQFAC = PQFAC/12.D0
      END IF
C
C  Come here if you intend to read a time series of flux data from
C  a disk file.
C
      IF(METU.EQ.2) THEN
C
      WRITE (6,117)
  117 FORMAT (1X,/,1X,'ENTER THE NAME OF THE THE AIR/SEA FLUX FILE')
      READ (5,116) IFILE
  116 FORMAT (A30)
      OPEN (UNIT=8,FILE=IFILE)
C
C  In the file read below, the array EMPA is the integrated evaporation
C  minus precipitation (m), not the rate as above. This is done to insure
C  that the net E-P is correctly applied to the model regardless of the
C  time step.
C
      NMETO = 0
      DO 233 J=1,NMET
      READ (8,*,END=459) DAYA(J), QIA(J), QLA(J), TXA(J), TYA(J),
     X EMPA(J)
C
      IF(DAYA(J).LT.-9.) GO TO 459
C
      NMETO = NMETO + 1
C
      IF(NMETO.LE.5) THEN
      WRITE (6,8838) NMETO, DAYA(J), QIA(J), QLA(J), TXA(J), TYA(J),
     X EMPA(J)
 8838 FORMAT (1X,'FIRST 5 FLUX DATA', I4, F7.2,2F7.0,2F7.2,E10.3)
      END IF
C
C
  233 CONTINUE
  459 CONTINUE
C
C  Set NMET equal to the actual number of met observations input
C
      NMET = NMETO
C
      WRITE (6,8833)
 8833 FORMAT (1X,/,1X,'ENTER THE STARTING TIME (DAY)')
      READ (5,*) DAY0
C
      END IF

C............Wave Spectral Forcing Input Set-up
C.....Read Number of wave forcing frequencies
      READ (5,*) NFR
C
      IF(NFR.GT.0) THEN
      NFRP=NFR+1
C
C  Come here if you intend to read a time series of flux data from
C  a disk file.
C
      WRITE (6,119)
  119 FORMAT (1X,/,1X,'ENTER NAME OF SCALAR WAVE SPECTRUM x DF FILE')
      READ (5,118) IFILE
  118 FORMAT (A30)
      OPEN (UNIT=12,FILE=IFILE)

      WRITE (6,121)
  121 FORMAT (1X,/,1X,'ENTER NAME OF SCALAR STOKES SPECTRUM x DF FILE')
      READ (5,120) IFILE
  120 FORMAT (A30)
      OPEN (UNIT=14,FILE=IFILE)

      WRITE (6,123)
  123 FORMAT (1X,/,1X,'ENTER NAME OF STOKES DIRECTION SPECTRUM FILE')
      READ (5,122) IFILE
  122 FORMAT (A30)
      OPEN (UNIT=15,FILE=IFILE)

C
C  In the file read below, the first column is DAYW except for in the first
C  row containing the frequency vectors. Spectrum is already times df computed
C  trapezoidally with intervals cetered on the input frequecy
C
      READ (12,*) DUMMY, (FRQ(IFR), IFR=1,NFR)
      READ (14,*) DUMMY, (FRQ(IFR), IFR=1,NFR)
      READ (15,*) DUMMY, (FRQ(IFR), IFR=1,NFR)

      NWAVO = 0
      DO 234 J=1,NWAV
      READ (12,*,END=460) DAYW(J), (H2SPC(J,IFR), IFR=1,NFR)
      READ (14,*,END=460) DAYW(J), (STKSSPC(J,IFR), IFR=1,NFR)
      READ (15,*,END=460) DAYW(J), (STKDSPC(J,IFR), IFR=1,NFR)
C
      DO 255 IFR=1,NFR
      SNSTKDR(J,IFR)=SIN(STKDSPC(J,IFR))
      COSTKDR(J,IFR)=COS(STKDSPC(J,IFR))
  255 CONTINUE

      IF(DAYW(J).LT.-9.) GO TO 460
C
      NWAVO = NWAVO + 1
C
      IF(NWAVO.LE.5) THEN
      H2=0.D0
      SX=0.D0
      SY=0.D0
      DO 235 IFR=1,NFR
      H2=H2+H2SPC(J,IFR)
      SX=SX+STKSSPC(J,IFR)*SNSTKDR(J,IFR)
      SY=SY+STKSSPC(J,IFR)*COSTKDR(J,IFR)
  235 CONTINUE
      WRITE (6,8848) NMETO, DAYW(J), H2, SX, SY
 8848 FORMAT (1X,'FIRST 5 WAVE DATA', I4,1X,4F7.2)
      END IF
C
C
  234 CONTINUE
  460 CONTINUE
C
C  Set NWAV equal to the actual number of met observations input
C
      NWAV = NWAVO

      DO 225 IFR=1,NFR
      ATTSCL(IFR)=2.D0*((2.D0*RPI*FRQ(IFR))**2.D0)/GRAV
      DO 224 K=1,KBM1
C...See Harcourt & D'Asaro (2008) (Appendix?) for an explanation of the FILTFACT term
C...There are numerical issues for computing this when it gets very close to one, may be machine dependent, Making it 1 in this machine case.
      ATT(K,IFR)=EXP(-ATTSCL(IFR)*ZZP(K))
      IF((ATTSCL(IFR)*DZP(K)).LT.100.D0) THEN
      FILTFACT=SINH(ATTSCL(IFR)*DZP(K)/2.D0)/(ATTSCL(IFR)*DZP(K)/2.D0)
      ATT(K,IFR)=FILTFACT*EXP(-ATTSCL(IFR)*ZZP(K))
      ENDIF
c      WRITE (6,3377) IFR, K, ATT(K,IFR), FILTFACT, ATTSCL(IFR)
c 3377 FORMAT (1X,'STOKES DATA IFR=', 2I6, 5G10.5)
  224 CONTINUE
  225 CONTINUE

      END IF


C
C  Finished with specifying air/sea flux data.
C
C
      WSSURF=0.D0
      WUBOT=0.D0
      WVBOT=0.D0
C
      ROCP=4.2E+6
      DO 9 I=1,KBM1
      PROFIL(I)=(XM*(EXP(ETA1*Z(I)*H)-EXP(ETA1*Z(I+1)*H))
     1+(1.D0-XM)*(EXP(ETA2*Z(I)*H)-EXP(ETA2*Z(I+1)*H)))/DZ(I)/H/ROCP
  9   CONTINUE
C
C
      WRITE (6,8899)
 8899 FORMAT (1X,/,1X,'DO YOU WANT A SCREEN LISTING OF DATA ?',/,
     X 1X,'(0 = NO, N = EVERY Nth STEP (DEFAULT = 6)')
      IWRT = 6
      READ (5,*) IWRT
C
C
C  Section below stores some data onto a disk file for later
C  analysis or plotting.
C
      ISTOR = 1
      WRITE (6,3186)
 3186 FORMAT (1X,/,1X,'DO YOU WANT TO STORE THE DATA ON DISK ?',/,1X,
     X '(0 = NO, 1 = YES (DEFAULT)')
      READ (5,*) ISTOR
C
C
      ITFREQ = 20
      IZFREQ = 3
C
      IF(ISTOR.EQ.1) THEN
      WRITE (6,3187)
 3187 FORMAT (1X,/,1X,'ENTER THE DECIMATION RATE FOR TIME AND DEPTH',
     X /,1X,'(DEFAULT IS 20 AND 3)')
      READ (5,*) ITFREQ, IZFREQ
      END IF
C
C
C*********************************************************************
C      BEGIN TIME MARCH
C*********************************************************************
C
       DO 9000 ITER=ISTART,IEND
C
C
C
      TIME=TIME+DT1/3600.D0
C
C  Compute the solar insolation as a function of the time,
C  where TIMED is time since start in days.
C
C
      TIMED = TIME/24.D0
      DTD = DT1/8.64E4
      IF(TIMED.GT.DAYS) GO TO 499
C
      IF(METU.EQ.1) THEN
      P2 = 3.14159D0*2.D0
      TIMS = DMOD(TIMED,1.D0)
      TIMS = TIMS - 0.5D0
      QI = 0.D0
      IF(COS(P2*TIMS).GT.0.D0) THEN
      QI = QIMAX*COS(TIMS*P2/PQFAC)
      IF(QI.LT.0.D0) QI = 0.D0
      END IF
      EMP = 0.D0
      END IF
C
      IF(METU.EQ.2) THEN
C
      DAY = TIMED + DAY0
      CALL INTX (DAYA,QIA,NSOBS,NMET,DAY,QI)
      CALL INTX (DAYA,QLA,NSOBS,NMET,DAY,QL)
      CALL INTX (DAYA,TXA,NSOBS,NMET,DAY,TX)
      CALL INTX (DAYA,TYA,NSOBS,NMET,DAY,TY)
C
C  Compute the E-P rate from the integrated E-P in EMPA.
C
      DAYP = DAY + DTD
      CALL INTX (DAYA,EMPA,NSOBS,NMET,DAYP,EMP2)
      CALL INTX (DAYA,EMPA,NSOBS,NMET,DAY,EMP1)
      EMP = (EMP2 - EMP1)/(DTD*8.64E4)
C
      END IF
C
C
      IF(NFR.GT.0) THEN
C
      DAY = TIMED + DAY0

      H2=0.D0
      STKX0=0.D0
      STKY0=0.D0
      STKM0=0.D0
      DO 625 K=1,KB
      STKX(K)=0.D0
      STKY(K)=0.D0
  625 CONTINUE
      DO 645 IFR=1,NFR
      CALL INTX (DAYW,H2SPC(1,IFR),NWOBS,NWAV,DAY,DH2)
      CALL INTX (DAYW,STKSSPC(1,IFR),NWOBS,NWAV,DAY,DUS)
      CALL INTX (DAYW,SNSTKDR(1,IFR),NWOBS,NWAV,DAY,DSS)
      CALL INTX (DAYW,COSTKDR(1,IFR),NWOBS,NWAV,DAY,DCS)
c      WRITE (6,3366) IFR, DAY, DH2, DUS, DSS, DCS
c 3366 FORMAT (1X,'STOKES DATA IFR=', I6, 5F10.5)
      H2=H2+DH2
      STKX0=STKX0+DUS*DSS
      STKY0=STKY0+DUS*DCS
      STKM0=STKM0+DUS
      DO 635 K=1,KBM1
      STKX(K)=STKX(K)+DUS*DSS*ATT(K,IFR)
      STKY(K)=STKY(K)+DUS*DCS*ATT(K,IFR)
  635 CONTINUE
  645 CONTINUE
      STKX(KB)=STKX(KBM1)
      STKY(KB)=STKY(KBM1)
C.....Top production is half missing, so goose the uppermost one:
C.....This is assuming uniform DZ. Need to change by ration of DZZ/DZ
C??   STKX(1)=STKX(1)+0.5*STKX0
C??   STKY(1)=STKY(1)+0.5*STKY0
CXX   H2=4*H2
C... This is done by setting surface value of Q2

      STK0=SQRT(STKX0**2.D0+STKY0**2.D0)
C     STKE=STK0*EXP(-1.D0)
      STKE=STKM0*EXP(-1.D0)

      IF(STKE.GT.0.D0) THEN

      DSTK=0.D0
      KSTK=0
      STKL=STK0
      DO 655 K=1,KBM1
      STK=SQRT(STKX(K)**2.D0+STKY(K)**2.D0)
      IF(STK.GT.STKE) THEN
      DSTK=-ZZ(K)*H
      KSTK=K
      STKL=STK
      ENDIF

  655 CONTINUE

      STKN=SQRT(STKX(KSTK+1)**2.D0+STKY(KSTK+1)**2.D0)
      DZST=-H*ZZ(KSTK+1)-DSTK
      DSTK=DSTK+DZST*(STKE-STKL)/(STKN-STKL)


      ELSE
      DSTK=0.D0
      ENDIF


C.... specify 1.6 * significant wave height as wave breaking length scale
CXC
C??      ZRSW=DMAX1(ZRS,1.6*4*sqrt(H2))
         ZRSW=ZRS

C
      END IF
C
C
      TOX(ITER) = TX/1023.D0
      TOY(ITER) = TY/1023.D0
      RA(ITER) = QI
      HSURF(ITER) = QL
c      WRITE (6,3344) ITER, TIMS, RA(ITER), HSURF(ITER), TOX(ITER),
c     c   TOY(ITER)
 3344 FORMAT (1X,'FLUX DATA', I6, 5F10.5)
c
c
      WUSURF=TOX(ITER)
      WVSURF=TOY(ITER)
      WWSURF=(TOX(ITER)**2.D0+TOY(ITER)**2.D0)**.5D0
      WTSURF=HSURF(ITER)/ROCP
      WSSURF=EMP*S(1)
C
      CALL PROFQ(DT2)
C

      DO 325 K=1,KB
      Q2(K)=Q2(K)+.5D0*SMOTH*(Q2F(K)+Q2B(K)-2.D0*Q2(K))
      Q2B(K)=Q2(K)
 325  Q2(K)=Q2F(K)

      DO 335 K=1,KB
      Q2L(K)=Q2L(K)+.5D0*SMOTH*(Q2LF(K)+Q2LB(K)-2.D0*Q2L(K))
      Q2LB(K)=Q2L(K)
 335  Q2L(K)=Q2LF(K)
C
C
      CALL PROFT(TF,TB,WTSURF,DT2,RA(ITER))
      CALL PROFT(SF,SB,WSSURF,DT2,0.D0)
C
      DO 345 k=1,KB
      T(K)=T(K)+.5D0*SMOTH*(TF(K)+TB(K)-2.D0*T(K))
      TB(K)=T(K)
      T(K)=TF(K)
      S(K)=S(K)+.5D0*SMOTH*(SF(K)+SB(K)-2.D0*S(K))
      SB(K)=S(K)
      S(K)=SF(K)
  345 CONTINUE
C
C
      CALL DENS(T,S,RHO)
C
C

      DO 380 K=1,KBM1
C     UU(K) = UB(K) + DT2*(COR*V(K) - PGX)
C 380 VV(K) = VB(K) + DT2*(-1.D0*COR*U(K) - PGY)
C.... Adding Stokes drift coriolis force
      UU(K) = UB(K) + DT2*(COR*(V(K)+STKY(K)) - PGX)
  380 VV(K) = VB(K) + DT2*(-1.D0*COR*(U(K)+STKX(K)) - PGY)
C
      UU(KB)=UU(KBM1)
      VV(KB)=VV(KBM1)
C

C... Extra Stress from Stokes Shear * TKE is added explicitly as doesn't depend on u,v

      TAUXUP=0.D0
      TAUYUP=0.D0
      DO 390 K=1,KBM2

      TAUXDN=-KMS(K+1)*(STKX(K)-STKX(K+1))/(DZZ(K)*H)
      TAUYDN=-KMS(K+1)*(STKY(K)-STKY(K+1))/(DZZ(K)*H)

      UU(K) = UU(K) - DT2*(TAUXUP-TAUXDN)/(DZ(K)*H)
      VV(K) = VV(K) - DT2*(TAUYUP-TAUYDN)/(DZ(K)*H)

      TAUXUP=TAUXDN
      TAUYUP=TAUYDN

  390 CONTINUE

C
C
      CBC=DMAX1(.0025D0,.16D0/DLOG((ZZ(KBM1)-Z(KB))*H/.01D0)**2)
      CDB = CBC
      CBC=CBC*SQRT(UB(KBM1)**2+VB(KBM1)**2)
      CDB2 = CBC
C
C  TO DELETE BOTTOM B.L. SET CBC = 0. AS BELOW
C      CBC = 0.D0
C
      CALL PROFC(UU,UF,WUSURF,CBC,DT2)
      WUBOT1 = WUBOT
      WVBOT1 = WVBOT
      CALL PROFC(VV,VF,WVSURF,CBC,DT2)
      WUBOT2 = WUBOT
      WVBOT2 = WVBOT
CCC    print*,'UF(1),VF(1)', UF(1),VF(1)
C
C STORE BOTTOM STRESS
C
      WUBOT = WUBOT1
      WVBOT = WUBOT2
C
C
C      WRITE (11,7733) CDB,CDB2,WUBOT1,WVBOT1,WUBOT2,WVBOT2
C 7733 FORMAT (1X,'CDS, ETC', 6E10.3)
C
      DO 382 K=1,KB
      U(K)=U(K)+.5D0*SMOTH*(UF(K)+UB(K)-2.D0*U(K))
      V(K)=V(K)+.5D0*SMOTH*(VF(K)+VB(K)-2.D0*V(K))
      UB(K)=U(K)
      U(K)=UF(K)
      VB(K)=V(K)
  382 V(K)=VF(K)
C
C  AS A TEST, EXRAPOLATE ALL VARIABLES TO THE BOTTOM
C  (SEEMS FINE)
C
      T(KB) = T(KBM1)
      S(KB) = S(KBM1)
      Q2(KB) = 0.D0
      Q2L(KB) = 0.D0
C...B625      Q2(KB) = Q2(KBM1)
C...B625      Q2L(KB) = Q2L(KBM1)
      U(KB) = U(KBM1)
      V(KB) = V(KBM1)
C
C  CHECK MOMENTUM BUDGET FOR BOTTOM LAYER IN CASE WHERE F = 0.
C  (WORKS FINE)
C
C      KBB = KB - 15
C      UBSUM = 0.D0
C      DO 4445 J=KBB,KB
C      UBSUM = UBSUM + (UG - U(J))*DZP(J)
C 4445 CONTINUE
C
C      TXBSUM = TXBSUM + DT1*WUBOT
C
C  Evaluate the energy budget. Seems to make sense only with surface
C  stress (for now).
C
C  ecf is the kinetic energy
C  work is the work by the wind stress on the surface current
C  shrp is shear production
C  ecb is the kinetic energy at the previous time step
C
      IF(ITER.LE.2) GO TO 5555
      ECF=0.D0
      PGWORK = 0.D0
      DO 601 K=1,KB
      UE = 0.5D0*(UF(K) + UB(K))
      VE = 0.5D0*(VF(K) + VB(K))
C
      IF(K.EQ.(KB-2)) THEN
      UEB = UE
      VEB = VE
      END IF
C
      ECF = ECF + (UE**2 + VE**2)*DZ(K)*H/2.D0
      PGWORK = PGWORK - DT1*(UE*PGX + VE*PGY)*DZ(K)*H
601   CONTINUE
      WORK = DT1*(U(1)*TOX(ITER)+V(1)*TOY(ITER))
      WORK = WORK + DT1*(WUBOT*UEB + WVBOT*VEB)
C
      SHRP=0.D0
      BUOYP = 0.D0
C
      DO 602 K=2,KB
      AUD = (U(K-1) - U(K))/(H*DZ(K-1))
      AVD = (V(K-1) - V(K))/(H*DZ(K-1))
      RHODZ = (RHO(K-1) - RHO(K))/(H*DZ(K-1))
      BUOYP = BUOYP - 9.8D0*KH(K)*RHODZ*(H*DZ(K-1))
      SHRP = SHRP - KM(K)*(AUD**2 + AVD**2)*(H*DZ(K-1))
 602  CONTINUE
      SHRP = SHRP*DT1
      BUOYP = BUOYP*DT1
C
      BILEC1 = ECF - ECB
      BILEC2 = BILEC1 + (WORK + SHRP + PGWORK)
      BILECT = BILECT + DT1*BILEC2
      ECB = ECF
 5555 CONTINUE
C
 603  CONTINUE
C
C      WRITE (11,9344) BILEC1, WORK, PGWORK, SHRP, BUOYP, WUBOT, UE,
C     X  UBSUM, TXBSUM
C 9344 FORMAT (1X,/,1X,'E', 9E8.2)
C
C
C      NPLT = 9
C      ICON=1
C      MNPLT = ITER/NPLT
C      IF (ITER .EQ. MNPLT*NPLT) THEN
C      DO 401 K=1,KB
C        IF (ICON .EQ. 1) THEN
C         UX(MNPLT,K)=U(KB-K+1)
C         VX(MNPLT,K)=V(KB-K+1)
C         DX(MNPLT,K)=RHO(KB-K+1)
C         TTX(MNPLT,K)=T(KB-K+1)
C         EZ1(K)=FLOAT(K-1)*0.1D0
C        ELSE
C         UX(MNPLT,K)=U(K)*50.D0+FLOAT(MNPLT)*2.D0
C         VX(MNPLT,K)=V(K)*50.D0+FLOAT(MNPLT)*2.D0
C         DX(MNPLT,K)=RHO(k)*50.D0+FLOAT(MNPLT)*2.D0
C         TTX(MNPLT,K)=T(K)*50.D0+FLOAT(MNPLT)*2.D0
C         EZ1(K)=FLOAT(K-1)*.5D0
C        END IF
C 401  CONTINUE
C      END IF
C
C  Write out a few things on every IWRTth time step.
C
      IF(IWRT.LT.1) GO TO 1339
C
      IF(ITER.EQ.1) WRITE (6,444)
  444 FORMAT (1X,/,1X,'MODEL SOLUTION VARIABLES ARE:',/,1X,
     X '   M TIMED    QI     QL      TX   T(1)   S(1)   U(1) ',
     X '  V(1)  DML   TLD',/,1X)
C
      CALL MLDPTH(ZZ,RHO,Q2,KBM1,ZZMLD,TLD)
      ZZDMLD=-ZZMLD*H
      TLD = -TLD*H
C
      IF(MOD(ITER,IWRT).NE.1) GO TO 1339
      WRITE (6,133) ITER,TIMED,QI,QL,TX,T(1),S(1),
     C U(1),V(1),ZZDMLD,TLD
C 133 FORMAT (1X,I4,F6.2,2F7.0,F7.2,4F7.2,3F6.0)
C 133 FORMAT (1X,I4,F6.2,2F7.0,F7.2,4F7.2,3F7.2)
C 133 FORMAT (1X,I5,F6.2,2F7.0,F7.2,4F7.2,3F7.2)
  133 FORMAT (1X,I7,F9.2,2F7.0,F7.2,4F7.2,3F7.2)
C
 1339 CONTINUE
C
C  WRITE OUT DATA TO UNIT 11 EVERY TIME STEP
C
      TXB = 1023.D0*WUBOT
      TYB = 1023.D0*WVBOT
      TXS = 1023.D0*TOX(ITER)
      TYS = 1023.D0*TOY(ITER)
C
C
      IF(ITER.EQ.1) WRITE (11,305)
  305 FORMAT (1X,'ITER',2X,'TIMED',5X,'QI',5X,'QL',5X,'TX',5X,
     C 'TY',6X,'T',6X,'S',6X,'U',6X,'V',4X,'MLD',
     C 4X,'TXB',4X,'TYB',6X,'DEKE',4X,'BLWORK',4X,'PGWORK',4X,
     C 'SHRPRO')
      WRITE (11,500)   ITER,TIMED,RA(ITER), HSURF(ITER),
     C TXS,TYS,T(1),S(1),
     C U(1),V(1), ZZDMLD, TXB, TYB, BILEC1, WORK, PGWORK, SHRP
  500 FORMAT (1X,I4,F7.2, 2F7.0, 6F7.2, F7.0, 2F7.2, 4E10.3)
9001  CONTINUE
C
C  Store profile data on logical unit 13, if requested.
C
      IF(ISTOR.EQ.1) THEN
      IF(MOD(ITER,ITFREQ).EQ.1) THEN
      WRITE (13,3183)
 3183 FORMAT (1X)
      WRITE (13,3182)
 3182 FORMAT (1X,5X,'TIMED',9X,'Z',9X,'T',9X,'S',9X,'U',9X,
C    X 'V',8X,'Q2',8X,'Q2L',8X,'KH',8X,'KM',8X,'KQ2',8X,'KQ2L',
     X 'V',8X,'STKX',8X,'STKY',8X,'Q2',8X,'Q2L',8X,
C    X 'KH',8X,'KM',8X,'KQ2',8X,'KQ2L',8X,'SH',
     X 'KH',8X,'KM',8X,'KQ2',8X,'KQ2L',8X,'KMS',8X,'SH',
     X 8X,'SM',8X,'SS',5X,'BPROD',5X,'SPROD',5X,'VPROD',8X,'WW')
C    X 8X,'SM',5X,'BPROD',5X,'SPROD',5X,'VPROD',8X,'WW')
C    X 8X,'SH',8X,'SM',5X,'BPROD',5X,'SPROD',9X,'Z')
      DO 3189 K=1,KB,IZFREQ
      WRITE (13,3188) TIMED,ZZP(K),T(K),S(K),U(K),V(K),
     X STKX(K), STKY(K),
     X Q2(K), Q2L(K), KH(K), KM(K), KQ2(K), KQ2L(K), KMS(K),
C    X Q2(K), Q2L(K), KH(K), KM(K), KQ2(K), KQ2L(K),
     X SH(K), SM(K), SS(K), BPROD(K), SPROD(K), VPROD(K),WW(K)
C    X SH(K), SM(K), BPROD(K), SPROD(K), VPROD(K),WW(K)
C    X SH(K), SM(K), BPROD(K), SPROD(K), ZZP(K)

C.....write out same to unit 10 without headers
      WRITE (10,3188) TIMED,ZZP(K),T(K),S(K),U(K),V(K),
     X STKX(K), STKY(K),
     X Q2(K), Q2L(K), KH(K), KM(K), KQ2(K), KQ2L(K), KMS(K),
C    X Q2(K), Q2L(K), KH(K), KM(K), KQ2(K), KQ2L(K),
     X SH(K), SM(K), SS(K), BPROD(K), SPROD(K), VPROD(K),WW(K)
C    X SH(K), SM(K), BPROD(K), SPROD(K), VPROD(K),WW(K)

C3188 FORMAT (1X,6F10.2,8E10.3,F10.2)
C3188 FORMAT (8(1X,F10.4),12(1X,E11.3))
C3188 FORMAT (8(1X,F10.4),14(1X,E11.3))
C3188 FORMAT (4(1X,F12.6),18(1X,G12.5))
 3188 FORMAT (4(1X,F16.10),18(1X,G12.5))
 3189 CONTINUE
      END IF
      END IF
C
C
 9009 CONTINUE
C
C  Store the density profile if this is the first step.
C
      IF(ITER.EQ.1) THEN
      DO 2278 J=1,KB
      DO(J) = RHO(J)
 2278 CONTINUE
      END IF
C
C  This is the place to compute diagnostics.
C
C  Find max/min values on the first day.
C
      IF(ITER.EQ.1) DMLMIN = 100.D0
      IF(ITER.EQ.1) TMIN = T(1)
      IF(TIMED.LT.1.D0) THEN
      IF(T(1).LT.TMIN) TMIN = T(1)
      IF(T(1).GT.TMAX) TMAX = T(1)
      SPD = SQRT(U(1)**2 + V(1)**2)
      IF(SPD.GT.SPDMAX) SPDMAX = SPD
      END IF
C
C  Compute the net heat flux as a check on heat conservation.
C
      QNET = QNET + DT1*(RA(ITER) + HSURF(ITER))
C
C
 9000 CONTINUE
  499 CONTINUE
C
C  Compute heat and potential energy diagnostics.
C  DPE will be the change in potential energy over the
C  integration period, and HEAT is the change in heat
C  content. HEAT should nearly equal the time-integrated
C  surface heat flux, QNET, if heat has been conserved.
C
      DPE = 0.D0
      HEAT = 0.D0
C
      DO 8376 J=1,KBM1
      CALL INTX(ZO,TO,NZOBS,NTS,ZZP(J),TI)
      TANOM = T(J) - TI
      HEAT = HEAT + TANOM*DZP(J)
      IF(J.EQ.1) TANOMS = TANOM
C
      DANOM = RHO(J) - DO(J)
      DPE = DPE + 9.8D0*1000.D0*DANOM*ZZP(J)*DZP(J)
C
C      WRITE (13,4465) J,DO(J),RHO(J),DANOM,DPE
 4465 FORMAT (1X,'DO,D,DANOM,DPE',I4,4E12.4)
 8376 CONTINUE
C
      QNET = QNET/(1.023E3*4.183E3)
      WRITE (6,3877) HEAT, QNET, DPE
 3877 FORMAT (1X,/,1X,'HEAT, QNET, DPE ARE', 3F9.2)
C
C  End of diagnostics.
C
C
C          WRITE (10) UX
C          WRITE (10) VX
C          WRITE (10) DX
C          WRITE (10) TTX
C
C
   40 CONTINUE
C
C
      WRITE (6,870) TMIN,TMAX,DMLMIN,SPDMAX
  870 FORMAT (1X,/,1X,'TMIN, TMAX, DMLMIN, SPDMAX ARE',4F8.2)
C
C
      STOP
      END
C
C
      SUBROUTINE INTX(R,S,KB,N,RI,SI)
C
C  This subroutine interploates the array S to find the value
C  at the point RI, where R is the independent variable.
C  The interpolation is linear, and assumes that R increases
C  montonically, and that RI is within the range of R. If
C  RI is not within the range of R, SI is set equal to the
C  appropriate endpoint of S, and a warning is printed.

C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(KB), S(KB)
C
C  Check for out of range RI.
C
      IF(RI.LT.R(1)) THEN
      SI = S(1)
C     WRITE (6,20) RI, R(1)
  20  FORMAT (1X,/,1X,'WARNING, RI WAS BELOW R(1) IN CALL TO INTX',/,
     C 'RI AND R(1) WERE', 2E12.4)
      END IF
C
C
      IF(RI.GT.R(N)) THEN
      SI = S(N)
C     WRITE (6,21) RI, R(N)
  21  FORMAT (1X,/,1X,'WARNING, RI WAS ABOVE R(N) IN CALL TO INTX',/,
     C 'RI AND R(N) WERE', 2E12.4)
      END IF
C
C  Data are within bounds for an interpolation.
C
      NM = N - 1
      DO 2 K=1,NM
      IF(RI.GE.R(K).AND.RI.LE.R(K+1)) GO TO 1
    2 CONTINUE
C
      GO TO 9
C
    1 CONTINUE
      SI = (S(K)*(R(K+1) - RI) + S(K+1)*(RI - R(K)))/
     X (R(K+1) - R(K))
C
C
    9 CONTINUE
      RETURN
      END
C
C
c----------------------------------------------------------------------
      SUBROUTINE MLDPTH(ZZ,R,Q2,KBM1,ZZMLD,TLD)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(KB=250)
      DIMENSION ZZ(KB),R(KB),Q2(KB)
C
      REPS = 0.00002D0
      QEPS = 0.1D0
      KML = 1
      ZZMLD = 1.D0
C
      DO 100 K=1,KBM1
      IF(ABS(R(1) - R(K)).GT.REPS) GO TO 100
      ZZMLD = ZZ(K)
      KML = K
  100 CONTINUE
C
      QS = 0.D0
      DO 1 J=1,KML
      QS = QS + Q2(J)
    1 CONTINUE
      QA = QS/FLOAT(KML)
C
      TLD = -1.D0
      DO 2 J=2,KBM1
CRRH  IF(Q2(J).LT.(QEPS*QA)) GO TO 3
      IF((Q2(J).GT.(1.E-5)).AND.(Q2(J+1).LT.(1.E-5))) THEN
      TLD = ZZ(J)
      ENDIF
    2 CONTINUE
CRRH3 CONTINUE
C
      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE PROFQ(DT2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KMS,KM,KH,KQ2,KQ2L,L,KBG
      REAL*8 KAPPA,KN,KHMIN,KMMIN,KBGCON
      PARAMETER(KB=250)
      COMMON/BLKT/TF(KB),T(KB),TB(KB),SF(KB),S(KB),SB(KB)
      COMMON/BLKU/UF(KB),U(KB),UB(KB),VF(KB),V(KB),VB(KB)
      COMMON/BLKS/STKX(KB),STKY(KB),STKX0,STKY0
      COMMON/BLKQ/Q2F(KB),Q2(KB),Q2B(KB),WW(KB)
      COMMON/BLKQL/Q2LF(KB),Q2L(KB),Q2LB(KB)
      COMMON/BLKR/RHO(KB),WT(KB)
      COMMON/BLKWU/WUSURF,WVSURF,USTR2,WUBOT,WVBOT
      COMMON/BLKWT/WTSURF,WSSURF,PROFIL(KB)
      COMMON/BLKZ/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),UMOL
      COMMON/BLKH/H,DT
      COMMON/BLKA/A(KB),C(KB),VH(KB),VHP(KB)
      COMMON/BLKK/KMS(KB),KM(KB),KH(KB),KQ2(KB),KQ2L(KB),L(KB),KBG(KB)
      COMMON/INDS/KBM1,KBM2
      COMMON/MODEL/GM(KB),GH(KB),SS(KB),SM(KB),SH(KB),KN(KB),FZS(KB),
     1   SPROD(KB),BPROD(KB),PROD(KB),DTEF(KB),VPROD(KB),GV(KB),GS(KB)
      COMMON/CSTUR/A1,B1,A2,B2,C1,C2,C3,E1,E2,VONKAR,SEF,QMIN,
     +          SQ2,SQ2L,COF1,COF2,COF3,COF4,COF5,ZRS,ZRSW,
     +          WAVEFAC,WAGE,ZRB,GHMIN,GHMAX,KHMIN,KMMIN,
     +          CKBG,KBGCON,E3,E4,E5,E6,C1S,C2S,
     +          DS0,DS1,DS2,DM0,DM1,DM2,DM3,DH0,DH1,DH2,DH3

      DATA KAPPA/0.4D0/
      DATA GEE/9.806D0/
      DATA SMALL/1.E-10/,EPS/1.E-30/
C
      DO 8877 J=1,KB
      SM(J) = DM0
      SH(J) = DH0
      SS(J) = DS0
      GM(J) = 0.154D0
      GH(J) = 0.154D0
 8877 CONTINUE
C
      L(1)=Q2L(1)/DMAX1(Q2(1),SMALL)
      BPROD(1)=0.D0
      DO 10 K=2,KB
      L(K)=Q2L(K)/DMAX1(Q2(K),SMALL)

      BPROD(K)=GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H)

      DTEF(K)=SQRT(Q2(K))/(B1*L(K)+EPS)


   10 CONTINUE

      DO 100 K=2,KBM1
      A(K)=-DT2*(KQ2(K+1)+KQ2(K)+2.D0*UMOL)*.5D0/(DZZ(K-1)*DZ(K)*H*H)
      C(K)=-DT2*(KQ2(K-1)+KQ2(K)+2.D0*UMOL)*.5D0/(DZZ(K-1)*DZ(K-1)*H*H)
  100 CONTINUE
C*********************************************************************
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C        DT2*(KQ2*Q2')' - Q2*(2.D0*DT2*DTEF+1.D0) =Q2B
C
C********************************************************************
C...Standard surface value q2(z=0)=ustar**2*b1**(2/3), no breaking source
      CONST1=B1**(2.D0/3.D0)
      USTR2=SQRT(WUSURF**2+WVSURF**2)
      CLS2=USTR2**1.5D0+0.2D0*(STKX0*WUSURF+STKY0*WVSURF)
      CLS2=DMAX1(CLS2,0.3D0*USTR2**1.5D0)**(2.D0/3.D0)
C...Drichlet condition for wave breaking source
      CONST2=CONST1*(1.D0+3.D0*WAVEFAC*VONKAR/1.81D0)**(2.D0/3.D0)
      VH (1)=0.D0
C     VHP(1)=SQRT(WUSURF**2+WVSURF**2)*CONST2
C...WALT ??
C     VHP(1)=CLS2*CONST2
C     VHP(1)=Q2F(1)
CWALTB?
      VHP(1)=DMAX1(2*Q2F(2)-Q2F(3),USTR2*CONST1)
C     VHP(1)=DMAX1(Q2F(2),USTR2*CONST1)
CWALTB?
C...WALT ??
      Q2F(KB)=SQRT(WUBOT**2+WVBOT**2)*CONST1
      DO 101 K=2,KBM1
C     BPROD(K)=GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H)
C...  Grab gradient Ri
      KBG(K)=BPROD(K)
      BPROD(K)=-KH(K)*BPROD(K)
      SPROD(K)=((U(K)-U(K-1))**2+(V(K)-V(K-1))**2)/
     1(DZZ(K-1)*H)**2
C...ADDED Vortex Production VPROD
      VPROD(K)=((STKX(K)-STKX(K-1))*(U(K)-U(K-1))
     +         +(STKY(K)-STKY(K-1))*(V(K)-V(K-1)))/
     1(DZZ(K-1)*H)**2
C.....Tempted to modify the gradient Richardson number: Sprod > Sprod + Vprod
      KBG(K)=KBG(K)/(SPROD(K)+SMALL)
C     SPROD(K)=KM(K)*SPROD(K)
C...ADDED Stokes-gradient momentum flux shear Production SPROD
C...there are 2 cross terms, one contributing to VPROD
      SPROD(K)=KM(K)*SPROD(K)+KMS(K)*VPROD(K)
C...The extra cross product of momentum & stokes shear does NOT go into VPROD in WW source term, but goes into SPROD for production in down-windwave direction
      VPROD(K)=KM(K)*VPROD(K)
C...MORE PRODUCTION FROM STRESS COMPONENT PROPORTIONAL TO STOKES SHEAR
      VPROD(K)=VPROD(K)+KMS(K)*((STKX(K)-STKX(K-1))**2.D0
     +          +(STKY(K)-STKY(K-1))**2.D0)/(DZZ(K-1)*H)**2.D0

CWALT7...C... Diagnosis from ARSM
CWALT7...      WW(K)=Q2(K)*(1.D0-6.D0*A1/B1)/3.D0
CWALT7...      TSDIS=L(K)/SQRT(Q2(K)+EPS)
CWALT7...      WW(K)=WW(K)+6.D0*TSDIS*A1*(BPROD(K)+VPROD(K))
CWALT7...C...WALT version
CWALT7...      WW(K)=WW(K)-6.D0*TSDIS*A1*VPROD(K)*FZS(K)
CWALT7...C...WALT version
CWALT7...C...optional small correction for Coriolis, not includid in SMC ARSM
CWALT7...      WW(K)=WW(K)+6.D0*A1
CWALT7...     + *(KM(K)*(U(K)-U(K-1))+KMS(K)*(STKX(K)-STKX(K-1)))
CWALT7...     + *COR/(DZZ(K-1)*H)
CWALT7...C....

C 101 PROD(K)=SPROD(K)+BPROD(K)
  101 PROD(K)=SPROD(K)+BPROD(K)+VPROD(K)
      DO 102 K=2,KBM1
      VHP(K)=1.D0/(A(K)+C(K)*(1.D0-VH(K-1))-(2.D0*DT2*DTEF(K)+1.D0))
      VH(K)=A(K)*VHP(K)
      VHP(K)=(-2.D0*DT2*PROD(K)+C(K)*VHP(K-1)-Q2B(K))*VHP(K)
  102 CONTINUE
      DO 103 K=1,KBM1
      KI=KB-K
      Q2F(KI)=VH(KI)*Q2F(KI+1)+VHP(KI)
  103 CONTINUE
      DO 112 K=2,KBM1
      IF(Q2F(K).GT.SMALL) GO TO 112
      Q2F(K)=SMALL
  112 CONTINUE

C...Set Minimum Surface Value
C     Q2F(1)=DMAX1(SQRT(WUSURF**2+WVSURF**2)*CONST2,Q2F(1))
C1101 Q2F(1)=DMAX1(USTR2*CONST2,Q2F(1))
C...WALT ??
C     Q2F(1)=DMAX1(CLS2*CONST2,Q2F(1))
      Q2F(1)=DMAX1(DMAX1(2*Q2F(2)-Q2F(3),USTR2*CONST1),Q2F(1))
C     Q2F(1)=DMAX1(USTR2*CONST1,Q2F(1))
C...WALT ??

C...Ed. guess for wave-dependent surface roughness ZRSW, wave age WAGE:
C     ZRSW=DMAX1(ZRS,1E5*USTR2*(WAGE**1.5D0)/GEE)
C...about 1.8 times too high, it seems

      DO 1100 K=2,KBM1
      A(K)=-DT2*(KQ2L(K+1)+KQ2L(K)+2.D0*UMOL)*.5D0/(DZZ(K-1)*DZ(K)*H*H)
      C(K)=-DT2*(KQ2L(K-1)+KQ2L(K)+2.D0*UMOL)
     +         *.5D0/(DZZ(K-1)*DZ(K-1)*H*H)
 1100 CONTINUE

C*********************************************************************
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C        DT2*(KQ2L*Q2L')' - Q2L*(DT2*DTEF+1.D0) =Q2LB
C
C********************************************************************
      VH (1)=0.D0
C     VHP(1)=SQRT(WUSURF**2+WVSURF**2)*CONST2*ZRSW
C...WALT ??
C     VHP(1)=CLS2*CONST2*ZRSW
      VHP(1)=Q2F(1)*ZRSW
C...WALT ??
      Q2LF(KB)=SQRT(WUBOT**2+WVBOT**2)*CONST1*ZRB
      DO 1101 K=2,KBM1
      ZZD=-ZZ(K)*H
      ZD=-Z(K)*H

      STAB=GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H)
      VPSTAB=-((STKX(K)-STKX(K-1))*(U(K)-U(K-1))
     +        +(STKY(K)-STKY(K-1))*(V(K)-V(K-1)))/
     +        (DZZ(K-1)*H)**2

      STK0=SQRT(STKX0**2.D0+STKY0**2.D0)
CWALTB
CWALTB      VONKAR2=VONKAR*(1.D0+0.5D0*STK0/(sqrt(USTR2)+EPS))**(1./3.)
      FZSMOD=FZS(K)/(1.0+(1.0-6.0*A1/B1)/3.)
C     VONKAR2=VONKAR*(1.D0+0.5D0*STK0*FZSMOD/(sqrt(USTR2)+EPS))**(1./3.)
      VONKAR2=VONKAR
CWALTB
CWALT8
      DOCEAN=5.E3
CWALT8
      DTEF(K)=E2*DTEF(K)*(1.D0+E4*(L(K)*(DOCEAN+ZRSW+ZRB)/
     +         (VONKAR2*(ZD+ZRSW)*(DOCEAN-ZD+ZRB)))**2.D0)
CWALT8      DISL=DMIN1(L(K),VONKAR*(ZD+ZRSW))

C11   This one should really be done at higher level ...
CWALT8      DOCEAN=5.E3
CWALT8      DISL=DMIN1(DISL,VONKAR*(DOCEAN-ZD+ZRB))

CWALT8      IF(STAB.GT.0) THEN
CWALT8      DOZMI=SQRT(-GHMIN*Q2(K)/STAB)
CWALT8      ENDIF

C...B624 1101 PROD(K)=E1*SPROD(K) + E3*BPROD(K) + E6*VPROD(K)
CWALTB?
C      RSP=DMAX1(EPS,sqrt((GS(K)+2.*GV(K)+GM(K))*GM(K)))
C     RSP=DMIN1(1.D0,ABS(GV(K)+GM(K))/RSP)
C      RVP=DMAX1(EPS,sqrt((GS(K)+2.*GV(K)+GM(K))*GS(K)))
C     RVP=DMIN1(1.D0,ABS(GS(K)+GV(K))/RVP)
C 1101 PROD(K)=E1*RSP*RSP*SPROD(K)+E3*BPROD(K)+E6*RVP*RSP*VPROD(K)
C..WALT5?... 1101 PROD(K)=E1*RSP*RSP*SPROD(K)+E3*BPROD(K)+E6*RVP*RSP*VPROD(K)
C...WALT6?

CWALTB?
 1101 PROD(K)=E1*SPROD(K)+E3*BPROD(K)+E6*VPROD(K)
CWALTB?

CWALT9 1101 PROD(K)=E1*SPROD(K)+E3*BPROD(K)+E6*DMAX1(0.D0,VPROD(K))
CWALT9     +       +E1*DMIN1(0.D0,VPROD(K))
CWALT9?
C...WALT6? 1101 PROD(K)=DMAX1(E1*SPROD(K),DMAX1(0.D0,E6*VPROD(K)+E3*BPROD(K)))
C...WALT6?     +       +DMIN1(0.D0,E6*VPROD(K)+E3*BPROD(K))
C...WALT6?

      DO 1102 K=2,KBM1
      VHP(K)=1.D0/(A(K)+C(K)*(1.D0-VH(K-1))-(DT2*DTEF(K)+1.D0))
      VH(K)=A(K)*VHP(K)
      VHP(K)=(-DT2*L(K)*PROD(K)+C(K)*VHP(K-1)-Q2LB(K))*VHP(K)
 1102 CONTINUE
      DO 1103 K=1,KBM1
      KI=KB-K
      Q2LF(KI)=VH(KI)*Q2LF(KI+1)+VHP(KI)
 1103 CONTINUE
      DO 1112 K=2,KBM1
      IF(Q2LF(K).GT.SMALL*EPS) GO TO 1112
      Q2LF(K)=SMALL*EPS
 1112 CONTINUE

C...Set Minimum Surface Value
      Q2LF(1)=DMAX1(Q2F(1)*ZRSW,Q2LF(1))


C******************************************************************
C      THE FOLLOWING SECTION SOLVES FOR KM AND KH
C******************************************************************
      L(1)=Q2LF(1)/DMAX1(Q2F(1),SMALL)
      DBMAX=L(1)
CWALT8      KBMAX=0
      BFMAX=0.D0

CWALTB      DO 210 K=2,KB
      DO 210 K=2,KBM1
      ZD=-Z(K)*H
      L(K)=Q2LF(K)/DMAX1(Q2F(K),SMALL)
      STAB=(GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H))
      IF(STAB.GT.0.D0) THEN
      DOZMI=SQRT(-GHMIN*Q2F(K)/STAB)
      DBUOY=SQRT(Q2F(K)/STAB)
      IF(L(K).LT.DMAX1(DBUOY,ZD)) THEN
      DBMAX=DMAX1(DBMAX,L(K))
CWALT8      IF(L(K).EQ.DBMAX) THEN
CWALT8      KBMAX=K
CWALT8      ENDIF
      ENDIF
      IF(DBUOY.GE.DBMAX) THEN
      BFMAX=DMAX1(STAB,BFMAX)
      ENDIF
CWALT8      IF(L(K).GE.DBUOY) THEN
CWALT8      Q2LF(K)=DMIN1(Q2LF(K),Q2F(K)*DBMAX)
CWALT8      ENDIF
      DISL=DMIN1(L(K),(ZD+ZRSW))
      DOCEAN=5.E3
      DISL=DMIN1(DISL,(DOCEAN-ZD+ZRB))

CWALTB?      L(K)=DMIN1(L(K),DISL)
CWALTB?      L(K)=DMIN1(L(K),DOZMI)

      ENDIF

C..WALT7... diagnose WW with limited L's
C... Diagnosis from ARSM
      WW(K)=Q2(K)*(1.D0-6.D0*A1/B1)/3.D0
      TSDIS=L(K)/SQRT(Q2(K)+EPS)
      WW(K)=WW(K)+6.D0*TSDIS*A1*(BPROD(K)+VPROD(K))
      WW(K)=WW(K)-6.D0*TSDIS*A1*VPROD(K)*FZS(K)
CWALTB
C...WALT version
C...optional small correction for Coriolis, not includid in SMC ARSM
C(N/A)       WW(K)=WW(K)+6.D0*A1
C(N/A)     + *(KM(K)*(U(K)-U(K-1))+KMS(K)*(STKX(K)-STKX(K-1)))
C(N/A)     + *COR/(DZZ(K-1)*H)
C....
C...WALT7...


  210 CONTINUE

C...WALT version
      ZNRM=0.D0
      ZSCL=0.D0
      DO 207 K=2,KBM1
      ZNRM=ZNRM+DMAX1(0.D0,VPROD(K))
      ZD=-Z(K)*H
C 207 ZSCL=ZSCL+ZD*DMAX1(0.D0,VPROD(K))
  207 ZSCL=ZSCL+L(K)*DMAX1(0.D0,VPROD(K))

      ZSCL=ZSCL/DMAX1(SMALL,ZNRM)
C     print *, 'ZSCL=',ZSCL
      DO 209 K=2,KBM1

      ZD=-Z(K)*H
C...alternate fzS:
C     ZND=-VONKAR*ZD/ZSCL
C 209 FZS(K)=(1-ZND)/(1-ZND+ZND**2)
C...simplest:
C 209 FZS(K)=EXP(-0.4*ZD/ZSCL)
C... with echo-blocking: keeps this production from rotating into w2 via return to isotropy
CWALT9
CWALT9 209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*(1.D0+TANH(-0.4*ZD/ZSCL))
CWALT10
C  209 FZS(K)=EXP(-0.4*ZD/ZSCL)
C  209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*(1.D0+TANH(-0.4*ZD/ZSCL))
C...B2
C  209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*(1.D0+TANH(-0.25*ZD/ZSCL))
  209 FZS(K)=(1.D0+TANH(-0.25*ZD/ZSCL))
C...B2
C  209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*(1.D0+TANH(-0.25*ZD/ZSCL))
C  209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)
C     +      +(1.0+2.0*(1.0-6.0*A1/B1)/3.)*TANH(-0.3*ZD/ZSCL)
CWALT10
CWALT9
C 209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*EXP(-0.4*ZD/ZSCL)
C 209 FZS(K)=(1.0+(1.0-6.0*A1/B1)/3.)*EXP(-ZD/ZSCL)
C...WALT version

CWALTB!!
      DO 212 K=2,KBM1
CWALTB  212 GH(K)=L(K)**2/Q2F(K)*((U(K)-U(K-1))**2+(V(K)-V(K-1))**2)/
  212 GM(K)=L(K)**2/Q2F(K)*((U(K)-U(K-1))**2+(V(K)-V(K-1))**2)/
     1(-DZZ(K-1)*H)**2
CWALTB!!

      DO 2212 K=2,KBM1
 2212 GS(K)=L(K)**2/Q2F(K)*((STKX(K)-STKX(K-1))**2
     +       +(STKY(K)-STKY(K-1))**2)/(-DZZ(K-1)*H)**2

      DO 3212 K=2,KBM1
 3212 GV(K)=L(K)**2/Q2F(K)*((STKX(K)-STKX(K-1))*(U(K)-U(K-1))
     +   +(STKY(K)-STKY(K-1))*(V(K)-V(K-1)))/(-DZZ(K-1)*H)**2

      DO 214 K=2,KBM1
      STAB=(GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H))
      GH(K)=-(L(K)**2.D0)*STAB/DMAX1(Q2F(K),SMALL)

C...WALT6? version
C     DO 3214 K=2,KBM1
      GS(K)=GS(K)*(1.D0-FZS(K))**2
      GV(K)=GV(K)*(1.D0-FZS(K))
C3214 GV(K)=GV(K)*(1.D0-FZS(K))
C...WALT6? version

      GVMAX=0.024D0
CWALTB

      GH(K)=DMAX1(GH(K),GHMIN)
      GV(K)=DMAX1(GV(K),GHMIN)

      GH(K)=DMIN1(GH(K),GHMAX)
C     GV(K)=DMIN1(GV(K),(GHMAX-GH(K))*A2/A1)
C     GS(K)=DMIN1(GS(K),(GHMAX-GH(K))*A2/A1-GV(K))
      GV(K)=DMIN1(GV(K),GVMAX)
      GS(K)=DMIN1(GS(K),GVMAX)

CWALTB      IF (GV(K).GT.GVMAX) THEN
CWALTBC...reduce L^2 for purposes of comp
CWALTB      REDL2=GVMAX/GV(K)
CWALTB      GS(K)=REDL2*GS(K)
CWALTB      GM(K)=REDL2*GM(K)
CWALTB      GH(K)=REDL2*GH(K)
CWALTB      GV(K)=GVMAX
CWALTB      ENDIF

CWALTB      IF ((GH(K)+GV(K)).GT.GHMAX) THEN
CWALTBC...reduce L^2 for purposes of comp
CWALTB      REDL2=GHMAX/(GH(K)+GV(K))
CWALTB      GS(K)=REDL2*GS(K)
CWALTB      GM(K)=REDL2*GM(K)
CWALTB      GV(K)=REDL2*GV(K)
CWALTB      GH(K)=REDL2*GH(K)
CWALTB      ENDIF

CWALT6? C...WALT version
CWALT6? C     DO 3214 K=2,KBM1
CWALT6?       GS(K)=GS(K)*(1.D0-FZS(K))**2
CWALT6?       GV(K)=GV(K)*(1.D0-FZS(K))
CWALT6? C3214 GV(K)=GV(K)*(1.D0-FZS(K))
CWALT6? C...WALT version

C1208

      DS1=1.D0-9.D0*A1*(A2*GH(K)+A1*GV(K))
      DS2=-9.D0*A1*A2*C2S*GH(K)

      DM1=1.D0-9.D0*A1*(A2*GH(K)+4.D0*A1*GV(K))
      DM2=9.D0*A1*(2.D0*A1+A2*(1.D0-C2))*GH(K)
      DM3=27.D0*A1*A1*GS(K)

      DH1=1.D0-3.D0*A2*((6.D0*A1+B2*(1-C3))*GH(K)
     +    +3.D0*A2*(1.D0-C2)*GV(K)-3.D0*A2*C2S*GS(K) )
      DH2=9.D0*A2*GV(K)*(2.D0*A1+A2)
      DH3=9.D0*A2*GS(K)*(2.D0*A1+A2)

      SH(K)=(DH0*DM1*DS1+DH2*DM0*DS1+(DH3*DM1+DM3*DH2)*DS0)
     +      /(DH1*DM1*DS1-DH2*DM2*DS1-(DH3*DM1+DM3*DH2)*DS2)
      SH(K)=DMAX1(0.D0,SH(K))

      SS(K)=(DS0+DS2*SH(K))/DS1
      SS(K)=DMAX1(0.D0,SS(K))

      SM(K)=(DM0+DM2*SH(K)+DM3*SS(K))/DM1
      SM(K)=DMAX1(0.D0,SM(K))

  214 CONTINUE

C     TRIGGER_KBG=DMAX1(1.E-4,USTR2)*1.E-2
C     TRIGGER_KBG=1.E-6
      TRIGGER_KBG=0.2D0

      DO 216 K=1,KBM1
      L(K)=Q2LF(K)/DMAX1(Q2F(K),SMALL)
      KN(K)=L(K)*SQRT(Q2F(K))
CWALT10
CWALT10      KMS(K)=(KN(K)*SS(K)+KMS(K))*.5D0
      KMS(K)=(KN(K)*SS(K)*(1.D0-FZS(K))+KMS(K))*.5D0
CWALT10
      KM(K)=(KN(K)*SM(K)+KM(K))*.5D0
      KH(K)=(KN(K)*SH(K)+KH(K))*.5D0

      RI=KBG(K)
      IF(RI.GE.0.D0) THEN
      KBG(K)=1.D0-(KBG(K)/0.7D0)**2.D0
      KBG(K)=DMAX1(0.D0,KBG(K))
      KBG(K)=CKBG*KBG(K)**3.D0
      ELSE
      KBG(K)=KBGCON
      ENDIF

      IF(K.GE.2) THEN
      STAB=(GEE*(RHO(K)-RHO(K-1))/(DZZ(K-1)*H))
      IF(STAB.GT.BFMAX) THEN
      KM(K)=DMAX1(KM(K),KMMIN+KBG(K))
      KH(K)=DMAX1(KH(K),KHMIN+KBG(K))
      KMS(K)=DMAX1(KMS(K),KHMIN+KBG(K))
C...WALT4      XKBG=KM(K-1)*EXP(-3.D0*DZ(K)*H/DBMAX)
C...WALT4      KM(K)=DMAX1(KM(K),XKBG)
C...WALT4      XKBG=KH(K-1)*EXP(-3.D0*DZ(K)*H/DBMAX)
C...WALT4      KH(K)=DMAX1(KH(K),XKBG)
C...WALT4      XKBG=KMS(K-1)*EXP(-3.D0*DZ(K-1)*H/DBMAX)
C...WALT4      KMS(K)=DMAX1(KMS(K),XKBG)
      ENDIF
      ENDIF

C...WALT4 version
C     KQ2(K)=SQ2*KH(K)
C     KQ2L(K)=SQ2L*KH(K)
C...WALT5 version
C     KQ2(K)=0.2*KN(K)
C     KQ2L(K)=0.2*KN(K)
CWALT6?     KQ2(K)=DMAX1(0.2D0*KN(K),SQ2*KH(K))
      KQ2(K)=(KQ2(K)+SQRT((0.2D0)**2.+(SQ2*SH(K))**2.)*KN(K))*0.5
C     KQ2(K)=DMAX1(KQ2(K),0.5*KM(K))
CWALT6?      KQ2L(K)=DMAX1(0.2D0*KN(K),SQ2L*KH(K))
      KQ2L(K)=(KQ2L(K)+SQRT((0.2D0)**2.+(SQ2L*SH(K))**2.)*KN(K))*0.5
C     KQ2L(K)=DMAX1(KQ2L(K),0.5*KM(K))
C...WALT5 version
C...WALT4 version

C...WALT version
CWALT10    KMS(K)=KMS(K)*(1.D0-FZS(K))
C...WALT version

  216 CONTINUE

      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE PROFT(FF,FB,WFSURF,DT2,ra)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KMS,KM,KH,KQ2,KQ2L,L,KBG
      REAL*8 KAPPA
      PARAMETER(KB=250)
      COMMON/BLKT/TF(KB),T(KB),TB(KB),SF(KB),S(KB),SB(KB)
      COMMON/BLKWT/WTSURF,WSSURF,profil(kb)
      COMMON/BLKZ/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),UMOL
      COMMON/BLKH/H,DT
      COMMON/BLKA/A(KB),C(KB),VH(KB),VHP(KB)
      COMMON/BLKK/KMS(KB),KM(KB),KH(KB),KQ2(KB),KQ2L(KB),L(KB),KBG(KB)
      COMMON/INDS/KBM1,KBM2
      DIMENSION FF(KB),FB(KB)
      DATA KAPPA/0.4D0/
      UMOLPR=UMOL/13.4D0
C************************************************************************************
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C        DT2*(KH*FF')' -FF = FB
C*************************************************************************************
      do 96 k=1,kbm1
      fb(k)=fb(k)+profil(k)*ra*dt2
96    continue
      FF(KBM1)=FB(KBM1)
      DO 100 K=2,KBM1
      A(K-1)=-DT2*(KH(K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*H*H)
      C(K)=-DT2*(KH(K)+UMOLPR)/(DZ(K)*DZZ(K-1)*H*H)
  100 CONTINUE
C*************************************************************************************
C   INSTEAD OF HEAT FLUX ONE CAN IMPOSE A SURFACE TEMPERATURE
C   BY REPLACYNG THE NEXT TWO LINES BY VH(1)=0 and VHP(1)=F(1).
C   F(1) SHOULD BE SPECIFIED IN THE MAIN PROGAM
C***************************************************************************************
      VH(1)=A(1)/(A(1)-1.D0)
      VHP(1)=(-DT2*WFSURF/(DZ(1)*H)-FB(1))/(A(1)-1.D0)
   98 CONTINUE
      DO 101 K=2,KBM2
      VHP(K)=1.D0/(A(K)+C(K)*(1.D0-VH(K-1))-1.D0)
      VH(K)=A(K)*VHP(K)
      VHP(K)=(C(K)*VHP(k-1)-FB(K))*VHP(K)
  101 CONTINUE
C************************************************************************************
C   INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB),
C   ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM
C   COL. 1 OF THE NEXT TWO LINES.
C*************************************************************************************
C     FF(KBM1)=(C(KBM1)*VHP(KBM2)-FB(KBM1))
C    1   /(C(KBM1)*(1.D0-VH(KBM2))-1.D0)
   99 CONTINUE
      DO 102 K=2,KBM1
      KI=KB-K
      FF(KI)=VH(KI)*FF(KI+1)+VHP(KI)
  102 CONTINUE
      RETURN
      END
c----------------------------------------------------------------------
      SUBROUTINE PROFC(CB,CF,WSURF,CBC,DT2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KMS,KM,KH,KQ2,KQ2L,L,KBG
      PARAMETER(KB=250)
      COMMON/BLKU/UF(KB),U(KB),UB(KB),VF(KB),V(KB),VB(KB)
      COMMON/BLKS/STKX(KB),STKY(KB),STKX0,STKY0
      COMMON/BLKWU/WUSURF,WVSURF,USTR2,WUBOT,WVBOT
      COMMON/BLKZ/Z(KB),ZZ(KB),DZ(KB),DZZ(KB),UMOL
      COMMON/BLKH/H,DT
      COMMON/BLKA/A(KB),C(KB),VH(KB),VHP(KB)
      COMMON/BLKK/KMS(KB),KM(KB),KH(KB),KQ2(KB),KQ2L(KB),L(KB),KBG(KB)
      COMMON/INDS/KBM1,KBM2
      DIMENSION  CF(KB),CB(KB)
      DO 100 K=2,KBM1
      A(K-1)=-DT2*(KM(K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*H*H)
      C(K)=-DT2*(KM(K)+UMOL  )/(DZ(K)*DZZ(K-1)*H*H)
  100 CONTINUE
      VH(1)=A(1)/(A(1)-1.D0)
      VHP(1)=(-DT2*WSURF/(DZ(1)*H)-CB(1))/(A(1)-1.D0)
      DO 101 K=2,KBM2
      VHP(K)=1.D0/(A(K)+C(K)*(1.D0-VH(K-1))-1.D0)
      VH(K)=A(K)*VHP(K)
      VHP(K)=(C(K)*VHP(k-1)-CB(K))*VHP(K)
  101 CONTINUE
      CF(KBM1)=(C(KBM1)*VHP(KBM2)-CB(KBM1))/(CBC
     1     *DT2/(-DZ(KBM1)*H)-1.D0-(VH(KBM2)-1.D0)*C(KBM1))
      DO 103 K=2,KBM1
      KI=KB-K
      CF(KI)=VH(KI)*CF(KI+1)+VHP(KI)
  103 CONTINUE
   92 WUBOT=-CBC*CF(KBM1)
      RETURN
      END

c----------------------------------------------------------------------

      SUBROUTINE DENS(TI,SI,RHOO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(KB=250)
      DIMENSION TI(KB),SI(KB),RHOO(KB)
C
C   ......... THIS SUBROUTINE COMPUTES (DENSITY-1)*0.001 .............
C
      DO 61 K=1,KB
      RHOO(K)=SGT(TI(K),SI(K),SG)*1.E-3
 61   CONTINUE
      RETURN
      END
c----------------------------------------------------------------------

      SUBROUTINE DEPTH(Z,ZZ,DZ,DZZ,KB,KBM1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(KB),ZZ(KB),DZ(KB),DZZ(KB)
      KL1=3
C
CRRH...don't use this option with large KB
C     KL1=INT(0.3D0*KB)
      KL2=KB-2
C**********************************************************************
C   THIS SUBROUTINE ESTABLISHES THE VERTICAL RESOLUTION WITH LOG
C   DISTRIBUTIONS AT THE TOP AND BOTTOM AND A LINEAR DISTRIBUTION
C   BETWEEN. CHOOSE KL1 AND KL2. THE DEFAULT KL1 = .3*KB AND KL2 = KB-2
C   YIELDS A LOG DISTRIBUTION AT THE TOP AND NONE AT THE BOTTOM.
C   KL1=3 AND KL2=KB-2 YIELDS NONE LOG AT THE TOP AND NONE AT BOT.
C**********************************************************************
      BB=FLOAT(KL2-KL1)+4.D0
      CC=FLOAT(KL1)-2.D0
      DEL1=2.D0/BB/EXP(.693147D0*FLOAT(KL1-2))
      DEL2=2.D0/BB/EXP(.693147D0*FLOAT(KB-KL2-1))
      Z(1)=0.D0
      ZZ(1)=-DEL1/2.D0
      DO 3 K=2,KL1
      Z(K)=-DEL1*EXP(.693147D0*FLOAT(K-2))
    3 ZZ(K)=-DEL1*EXP(.693147D0*(FLOAT(K)-1.5D0))
      DO 4 K=KL1,KL2
      Z(K)=-(FLOAT(K)-CC)/BB
    4 ZZ(K)=-(FLOAT(K)-CC+0.5D0)/BB
      DO 5 K=KL2,KBM1
      Z(K)=(1.D0-DEL2*EXP(.693147D0*FLOAT(KB-K-1)))*(-1.D0)
    5 ZZ(K)=(1.D0-DEL2*EXP(.693147D0*(FLOAT(KB-K)-1.5D0)))*(-1.D0)
      Z(KB)=-1.D0
      ZZ(KBM1)=-1.D0*(1.D0-DEL2/2.D0)
      ZZ(KB)=-1.D0*(1.D0+DEL2/2.D0)
      DO 6 K=1,KBM1
      DZ(K)=Z(K)-Z(K+1)
    6 DZZ(K)=ZZ(K)-ZZ(K+1)
      RETURN
      END
c----------------------------------------------------------------------
C
C
      FUNCTION SG0(S)
C
C  A sigma-0 subroutine neede by the sigma-t subroutine
C  taken from SEAPROP.
C
C  SIGMA-0 KNUDSEN
      IMPLICIT REAL*8 (A-H,O-Z)
      SG0 = ((6.76786136E-6*S-4.8249614E-4)*S+0.814876577D0)*S
     X -0.0934458632D0
      RETURN
      END
C
C
      FUNCTION SGT(T,S,SG)
C
C  A sigma-t subroutine taken from SEAPROP.
C
C SIGMA-T KNUDSEN
      IMPLICIT REAL*8 (A-H,O-Z)
      SG = SG0(S)
   20 SGT = ((((-1.43803061E-7*T-1.98248399E-3)*T-0.545939111D0)*T
     X +4.53168426D0)*T)/(T+67.26D0)+((((1.667E-8*T-8.164E-7)*T
     X +1.803E-5)*T)*SG+((-1.0843E-6*T+9.8185E-5)*T-4.7867E-3)*T
     X +1.D0)*SG
      RETURN
      END
