C     ******************************************************************
      PROGRAM MAIN
C     ******************************************************************
      DIMENSION Q(160),SRATE(160),DEPSIT(3,8),TOUT(10),DIAM(21)
c     ------------------------------------------------------------------
      DIMENSION DPI(20),V(21)
      CHARACTER*40 OUTFILE
C     Varaiables in the above two lines:
C        DPI - Array storing the geometric mean particle dia of section
C        V   - Array storing section limit sizes (i.e. particle volumes)
C        OUTFILE - variable storing output file name.
C     ------------------------------------------------------------------
      COMMON/PHYSPT/ACELOV,AFLROV,AWALOV,CHI,DELDIF,DELSAT,DENSTY,FSLIP,
     $FTHERM,GAMMA,PSAT,STICK,TGRADC,TGRADF,TGRADW,TKGOP,TURBDS,VOLUME,
     $WTCONM,WTMOL
C     ------------------------------------------------------------------
C     Specify default values for variables
C     ------------------------------------------------------------------
      DATA TGAS1,TGAS2,PGAS1,PGAS2/298.,450.,1.01E5,7.09E5/
      DATA M,KCOMP,IPRNT,IPLOT,ISCAL,NROW,NCOL,QMINPL,QMAXPL
     $/11,1,1,1,0,50,101,1.E-6,1./
      ACELOV=1.07
      AFLROV=1.07
      AWALOV=2.29
      CHI= 1.0
      DELDIF=1.E-5
      DELSAT= 1.E-11
      DENSTY= 1.E3
      FSLIP= 1.37
      FTHERM= 1.
      GAMMA=1.
      PSAT= 3.16E3
      STICK = 1.
      TGRADC=0.
      TGRADF=0.
      TGRADW= 0.
      TKGOP= .05
      TURBDS= .001
      VOLUME= 2.86
      WTCONM= 3.0E-26
      WTMOL= 28.8
C     ------------------------------------------------------------------
C     Setup Output file for storing output
      WRITE(*,30)
      READ(*,31)FILNAM1
   30 FORMAT(1X,'ENTER OUTPUT FILE NAME: ')
   31 FORMAT(A40)
C
      OPEN(UNIT=IPRNT,FILE=OUTFILE,STATUS='NEW')
C     ------------------------------------------------------------------
C     Set or reset the following variables to suit your problem.  You
C     can set them by assignment statements or by interactive input:
C  M - number of sections, should be less than 20
C       DENSTY - Particle density
C  KCOMP - number of components
C  DIAM(I) - Section boundary diameters       | Initial conditions
C  Q(I) - Mass concentration in each section  | specified here.
C     TOUT(I) - Array of output times, should be greater than 0
C     All other variables are already set
C     ------------------------------------------------------------------

C        SET SECTION BOUNDARIES



C        ------------------------------------------------------------------
C     SET TGAS=TGAS1 AND PGAS=PGAS1 TO USE SECTIONAL COEFFICIENTS
C     AT THESE TEMPERATURES AND PRESSURES IF NEEDED
C
      TGAS=TGAS1
      PGAS=PGAS1
      TIME=0.
      ROUND=6.0E-8
      IFLAG=1
      NEWCOF=2
C     ------------------------------------------------------------------
C     Call the Output subroutine and then MAEROS
C     ------------------------------------------------------------------
      CALL PRAERO(Q,DEPSIT,DIAM,TIME,VOLUME,M,KCOMP,IFLAG,IPRNT,IPLOT,
     $QMINPL,QMAXPL,ISCAL,NROW,NCOL,DPI,DENSITY)
      DO 3 ITIME=1,NTIME
      DELTIM=TOUT(ITIME)-TIME
      CALL MAEROS(M,KCOMP,TIME,DELTIM,Q,SRATE,TGAS,PGAS,DEPSIT,DIAM,
     $ TGAS1,TGAS2,PGAS1,PGAS2,ROUND,IPRNT,IFLAG,NEWCOF)
    3 CALL PRAERO(Q,DEPSIT,DIAM,TIME,VOLUME,M,KCOMP,IFLAG,IPRNT,IPLOT,
     $QMINPL,QMAXPL,ISCAL,NROW,NCOL,DPI,DENSITY)
C
      STOP
      END
