C     ******************************************************************
      PROGRAM MAIN
C     ******************************************************************
      DIMENSION Q(160),SRATE(160),DEPSIT(3,8),TOUT(50),DIAM(21)
      DIMENSION DPI(20),VI(21),SV(21)
      CHARACTER*40 OUTFILE
C     Varaiables in the above two lines:
C        DPI - Array storing the geometric mean particle dia of section
C        SV  - Array storing section limit sizes (i.e. particle volumes)
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
     +/10,1,10,10,-1,25,60,1.0E-6,1.0E-2/
      DATA NTIME /50/
      DO 50 I=1,NTIME
      TOUT(I)=2.0*I
  50  CONTINUE
      AFLROV=0.0
      AWALOV=0.0
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
      TURBDS=0.0
      VOLUME= 2.86
      WTCONM= 3.0E-26
      WTMOL= 28.8
C     ------------------------------------------------------------------
C     Setup Output file for storing output
      WRITE(6,30)
      READ(5,31) FILNAM1
   30 FORMAT(1X,'ENTER OUTPUT FILE NAME: ')
   31 FORMAT(A40)
C
      OPEN(UNIT=IPRNT,FILE='OUTFILE',STATUS='NEW')
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
      DIAM(1)=0.001E-6
      DIAM(M+1)=1.0E-6
      DO 10 I=2,M
      DIAM(I)=DIAM(1)*(DIAM(M+1)/DIAM(1))**(FLOAT(I-1)/FLOAT(M))
  10  CONTINUE
      DO 15 I=1,M+1
      SV(I)=(3.141592654*DIAM(I)**3)/6.0
  15  CONTINUE
      MKCOMP=M*KCOMP
      DO 20 I=1,MKCOMP
      Q(I)=0.0
      SRATE(I)=0.0
  20  CONTINUE
      CN0=5.0E14
      SV0=2.0E-28
      V0=CN0*SV0
      DO 40 I=1,M
      SVAVE=SQRT(SV(I)*SV(I+1))
      CON=CN0*SVAVE/V0
      IF(CON.GT.88) THEN
      Q(I)=0.0
      GO TO 40
      ELSE
      END  IF
      Q(I)=DENSTY*SVAVE*(CN0*CN0)*(EXP(-CON))*(SV(I+1)-SV(I))/V0
  40  CONTINUE
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
      CALL PRAERO(Q,DEPSIT,DIAM,TIME,VOLUME,M,KCOMP,IFLAG,IPRNT,IPLOT
     +,QMINPL,QMAXPL,ISCAL,NROW,NCOL,DPI,DENSTY)
      DO 100 ITIME=1,NTIME
      WRITE(*,1000) ITIME,IFLAG
 1000 FORMAT(2X,'ITIME = ',I4,5X,'IFLAG = ',I4)
      DELTIM=TOUT(ITIME)-TIME
      CALL MAEROS(M,KCOMP,TIME,DELTIM,Q,SRATE,TGAS,PGAS,DEPSIT,DIAM,
     + TGAS1,TGAS2,PGAS1,PGAS2,ROUND,IPRNT,IFLAG,NEWCOF)
      IFLAG=-1
      ISCAL=-1
      IF((ITIME.EQ.1).OR.(ITIME.EQ.10).OR.(ITIME.EQ.25).OR.
     +  (ITIME.EQ.50)) THEN
      CALL PRAERO(Q,DEPSIT,DIAM,TIME,VOLUME,M,KCOMP,IFLAG,IPRNT,IPLOT
     +,QMINPL,QMAXPL,ISCAL,NROW,NCOL,DPI,DENSTY)
      ELSE
      END IF
  100 CONTINUE
C----------------------------------------------------------------
      STOP
      END
