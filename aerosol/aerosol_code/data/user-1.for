      SUBROUTINE BETA(npts,time,xvolum,yvolum,coef)
c---------------------------------------------------------------
c*    this subroutine calculates the coagulation coefficients
c*    beta, as a function of the volumes of the coagulating p
c*    particles
c*
c*                         Variables
c*
c*    akundx=knudsen number of first particle
c*    akundy=knudsen number of second particle
c*    coef=coagulation coefficient (output array cc/sec
c*    density=particle material density (grams/cc)
c*    difx=diffusivity of first particle (cm**2/sec)
c*    dify=diffusivity of second particle (cm)
c*    dx=diameter of first particle (cm)
c*    dy=diameter of second particle  (cm)
c*    freemn=mean free path of air (centimeters)
c*    npts=number of points (input)
c*    tempk=temperature (degrees kelvin)
c*    time=time (input in seconds)
c*    viscos=viscosity of air (poise)
c*    xvolum=volume of first particle (input in cc)
c*    yvolum=volume of second particle (input array in cc)
c*
c*    expression taken from fuchs, N.A.,1964, pp291-294.D0
c-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YVOLUM(1),COEF(1)
      COMMON/AERSL1/IPARM(10)
      DATA TEMPK,DENSTY,VISCOS,FREEMN/298.D0,1.D0,1.83D-4,.066D-4/
c
      DX=(1.90986*XVOLUM)**0.33333333333
      AKNUDX=FREEMN*2./DX
      TKV=1.46E-17*TEMPK/VISCOS
      DIFX=(1.+AKNUDX*(1.257+.4*EXP(-1.1/AKNUDX)))*TKV/DX
      SPEEDX=3.51E-16*TEMPK/(XVOLUM*DENSTY)
      PMEANX=2.546479*DIFX/SQRT(SPEEDX)
      GX=DX+PMEANX
      GX1=DX*DX+PMEANX*PMEANX
      GX=(GX*GS*GX-GX1*SQRT(GX1))/(3.*DX*PMEANX)-DX
      GX=GX*GX
C
      DO 10 I=1,NPTS
      DY=(1.90986*YVOLUM(I))**.33333333333
      AKNUDY=FREEMN*2./DY
      DIFY=(1+AKNUDY*(1.257+.4*EXP(-1.1/AKNUDY)))*TKV/DY
      SPEEDY=SQRT(3.51E-16*TEMPK/(YVOLUM(I)*DENSTY))
      SPEED=SQRT(SPEEDX+SPEEDY*SPEEDY)
      PMEANY=2.546479*DIFY/SPEEDY
      GY=DY+PMEANY
      GY1=DY*DY+PMEANY*PMEANY
      GY=(GY*GY*GY-GY1*SQRT(GY1))/(3.*DY*PMEANY)-DY
      G=SQRT(GX+GY*GY)
      COEF(I)=(DX+DY)*.5
      COEF(I)=COEF(I)/(COEF(I)+G)+4.*(DIFX+DIFY)/(SPEED*COEF(I))
   10 COEF(I)=6.2832*(DIFX+DIFY)*(DX+DY)/COEF(I)
      RETURN
C
    1 DO 11 I=1,NPTS
   11 COEF(I)=XVOLUM-YVOLUM(I)
      RETURN
C
    2 DO 12 I=1,NPTS
   12 COEF(I)=0.
      RETURN 
      END
C-----------------------------------------------------------------------
      FUNCTION BOUNDR(TIME)
C-----------------------------------------------------------------------
C*      THIS IS AN EXAMPLE OF A BOUNDARY SOURCE FUNCTION.              *
C*                                                                     *
C*                  VARIABLES                                          *
C*                                                                     *
C*    BOUNDR=SOURCE AT BOUNDARY (PARTICLE VOLUME=VA) (1/CC/SEC/CC/)    *
C*    TOUT(I)=I-TH OUTPUT TIME (SECONDS)                               *
C*    VA=LOWER LIMIT OF CONTINUOUS DISTRIBUTION IN PARTICLE VOLUME (CC)*
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/AERSL6/TOUT(16),TDUMY,VA,VB,VBOVA,ALGVBA
      BOUNDR=.3
      IF(TIME.GT.TOUT(2)) BOUNDR=BOUNDR*EXP(.0005*(TOUT(2)-TIME))
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GROWTH(NPTS,TIME,V,RATE,PARTAL)
C-----------------------------------------------------------------------
C*     THIS IS AN EXAMPLE OF A GROWTH SUBROUTINE USING THE CONTINUUM   *
C*     GROWTH LAW. (REF. GELBARD,F. AND SEINFELD,J.H.,JOURNAL OF       *
C*     colloid AND INTERFACE SCIENCE, VOL.68,PAGES      173-183 (1979))     *
C*                                                                     *
C*                 VARIABLES                                           *
C*                                                                     *
C*    D=PARTICLE DIAMETER  (CM)                                        *
C*    DIFUS=MOLECULAR DIFFUSIVITY  (CM**2/SECOND)                      *
C*    DEL=CONCENTRATION DIFFERENCE  (1/CC)                             *
C*    NPTS=NUMBER OF POINTS  (INPUT)                                   *
C*    PARTAL=PARTIAL OF RATE W.R.T. PARTICLE VOLUME (OUTPUT IN 1/SEC)  *
C*    RATE=GROWTH RATE OF A PARTICLE OF VOLUME V (OUTPUT IN CC/SEC)    *
C*    V=PARTICLE VOLUME (INPUT IN CC)                                  *
C*    VMOLEC=MOLECULAR VOLUME  (CC)                                    *
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(1),RATE(1),PARTAL(1)
      COMMON/AERSL1/IPARM(10)
      DATA PI,DIFUS,VMOLEC/3.14159,.1,6.54E-23/
      DEL=1.E8*EXP(-.001*TIME)
      DO I I=1,NPTS
      D=(1.9099*V(I))**.33333333333
      RATE(I)=2.*PI*D*DIFUS*DEL*VMOLEC
      PARTAL(I)=RATE(I)/(3.*V(I))
    1 IF(IPARAM(3).EQ.1) RATE(I)=-2.*RATE(I)
      RETURN
      END DO
C-----------------------------------------------------------------------
      SUBROUTINE XINTL(NPTS,V,DISTB)
C-----------------------------------------------------------------------
c*      this is an example of an initial distribution function.
c*      A sum of three log-normal functions or a logarithmic
c*      interpolation between data points can be used depending
c*      on iparm(2)
c*
c*                        variables
c*      distb=output array of initial distributions (no./cc/cc)
c*      npts=number of points (input)
c*      v=input array of particle volumes (cc)
c*-------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DP(19),F(19),V(1),DISTB(1)
      COMMON/AERSL1/IPARM(10)
      DATA DP/.0133,.0237,.0421,.075,.133,.237,.42,.513,.61,.704,.784,
     $ .847,.902,.953,1.01,1.07,1.16,1.28,1.41/
      DATA F/2.38E7,5.68E6,2.76E5,1.7E5,3.54E4,2.1E3,1.9E2,70.7,24.8,
     $ 13.9,10.5,7.78,5.9,4.57,3.22,2.36,1.58,.714,.705/
      DATA ICALL/0/

      ALNORM(TOT,GMEAN,GLSIG,VX)=.3989*TOT*EXP(-.5*(ALOG(VX/GMEAN)/
     $ GLSIG)**2)/(VX*GLSIG)

      IF(IPARM(2).EQ.1) GO TOUI 112
      DO 111 I=1,NPTS
 111  DISTB(I)=ALNORM(1.06E5,1.44E-18,1.76,V(I))+
     $ ALNORM(3.2E4,8.24E-17,2.31,V(I))+ALNORM(5.4,3.33E-13,2.38,V(I))
      RETURN

 112  IF(ICALL.EQ.1) GO TO 122
      ICALL=1
      DO 123 I=1,19
      DP(I)=ALOG(F(I))
 123  F(I)=ALOG(F(I))
 122  DO 3 NP=1,NPTS
      D=(6.E12*V(NP)/3.14159)**(.33333333333)
      DL=ALOG(D)
      DO 124 I=2,19
      IF(DL.LT.DP(I)) GO TO 125
 124  CONTINUE
 125  DL=E(I)+(DP(I)-DL)/(DP(I)-DP(I-1))*(F(I-1)-F(I))
   3  DISTB(NP)= EXP(DL)*D/(3.*V(NP))
      RETURN 
      END
