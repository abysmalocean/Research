        SUBROUTINE BETA(npts,time,xvolum,yvolum,coef)
c---------------------------------------------------------------
c*      this subroutine calculates the coagulation coefficients
c*      beta, as a function of the volumes of the coagulating p
c*      particles
c*
c*                         Variables
c*
c*      akundx=knudsen number of first particle
c*      akundy=knudsen number of second particle
c*      coef=coagulation coefficient (output array cc/sec
c*      density=particle material density (grams/cc)
c*      difx=diffusivity of first particle (cm**2/sec)
c*      dify=diffusivity of second particle (cm)
c*      dx=diameter of first particle (cm)
c*      dy=diameter of second particle  (cm)
c*      freemn=mean free path of air (centimeters)
c*      npts=number of points (input)
c*      tempk=temperature (degrees kelvin)
c*      time=time (input in seconds)
c*      viscos=viscosity of air (poise)
c*      xvolum=volume of first particle (input in cc)
c*      yvolum=volume of second particle (input array in cc)
c*
c*      expression taken from fuchs, N.A.,1964, pp291-294.D0
c-----------------------------------------------------------------------
        IMPLICIT REAL*8 (A-H,O-Z)
        dimension yvolum(1),coef(1)
        common/aersl1/iparm(10)
        data tempk,densty,viscos,freemn/298.D0,1.D0,1.83D-4,.066D-4/

        end
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
      common/aersl1/iparm(10)

      END
C-----------------------------------------------------------------------
      SUBROUTINE XINTL(NPTS,V,DISTB)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
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
      dimension X(30),f(30),v(1),distb(1)
      common/aersl1/iparm(10)

      END
