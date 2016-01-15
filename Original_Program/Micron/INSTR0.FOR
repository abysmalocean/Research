c
c       
c      THIS IS THE SUBROUTINE THAT SIMULATES THE CNC RESPONSE TO
c      PARTICLES FOR MICRON
c
c      NOTE: ONLY CHANNELS 375 - 499 ARE USED FOR THE INVERSION
c
c                                Spyros N. Pandis       
c                        Department of Chemical Engineering
c                           Carnegie Mellon University
c                             Pittsburgh, PA 15146
c
c      LAST MODIFIED : 5 JANUARY 1994 BY S.P.
c                      9 FEBRUARY 1994 BY S.P.
c                      3 MARCH 1994 by S.P.
c                     16 March 1994 by S.P. 
c                      7 July 1994 (new calibration data)              
C                     19 JUNE 1996 (STOLZENBURG ANALYTICAL EXPRESSIONS-BY RW 
c
c......VARIABLES
c      jcha : THE CHANNEL FOR WHICH THE RESPONSE IS DESIRED
c
c
       double precision function INST1(diam)  
c       
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       integer iinst, jcha, kcha
       integer minst,minstt,mpar,mpart,msg,nlog,mcha
       logical first
       double precision diam, cncmax, par,filt,FILTER,rmin
       COMMON/COEF/COF(4,3),PRM(3),DPP
       parameter (mpar=6)
       parameter (minst=10)
       PARAMETER(DPLGE=50.0)
       PARAMETER(CHLG=480.0)
       PARAMETER(ICHN0=375)
       PARAMETER(CNCMAX=16.0)
       PARAMETER(CNCMIN=2.5)
       common /cmkerf/ iinst, jcha, mcha
       save   /cmkerf/
       common /cmrmin/ rmin
       save   /cmrmin/
       common /cmpar/ par(minst,mpar)
       save   /cmpar/
       common /cmparm/ minstt, mpart
       save   /cmparm/ 
       save msg, nlog, first
       data msg,nlog,first / 0, 6, .true. /
C      
       if (first) then 
         first = .false.
         call PARCHK(minst, minstt,'minst','INST1')
         call PARCHK(mpar,mpart,'mpar','INST1')
         endif
c
c      WRITE A WARNING IF THE KERNEL FUNCTION IS COMPUTED WHEN
c      CALIBRATION DATA IS UNAVAILABLE
c       
c       if (diam .gt. cncmax) then
c         if (msg .eq. 0) then
c           write(nlog,1)cncmax 
c           msg=1
c         endif
c       endif
c       
1      format(/,' Warning:',/,
     & ' The instrument response for the ultrafine CNC',/,
     & ' is diameter independent above ',d10.3, ' nm')
c
c
C      
      JCH=ICHN0+JCHA-1
      P=PCHAN(CHLG,DPLGE)
      CHOF=CHLG-PRM(1)
      IF(ABS(PRM(3)).GE.1.D-8) CHOF=CHOF-(PRM(2)/PRM(3))
     1   *LOG(0.5D0*(1.D0+SQRT(1.D0+4.D0*PRM(3)**2)))
      COF(1,1)=COF(1,1)+CHOF
      DPP=0.D0
C
      CH=DBLE(JCH)
      INST1=PCHAN(CH,DIAM)
      RETURN
      END
C
C*******************************************************************      
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/COEF/COF(4,3),PRM(3),DPP
C
C values for coef2.prm
C
      DATA COF(1,1),COF(2,1),COF(3,1),COF(4,1)/482.0894549,
     & -513.4098062,-0.0230800,1.8770826/
      DATA COF(1,2),COF(2,2),COF(3,2),COF(4,2)/5.0569115,
     & 16.8464950,1.7770725,1.8688531/
      DATA COF(1,3),COF(2,3),COF(3,3),COF(4,3)/0.0142593,
     & 283.3491562,-0.6487761,4.8220668/
C
      END
C
C*************************************************************
      FUNCTION PCHAN(CHAN,DPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/COEF/COF(4,3),PRM(3),DPP
      DATA P/0.3275911D0/A1/0.254829592D0/A2/-0.284496736D0/
     1 A3/1.421413741D0/A4/-1.453152027D0/A5/1.061405429D0/
C
      IF(DPN.NE.DPP) THEN
        DO 10 K=1,3
          PRM(K)=COF(1,K)+COF(2,K)/(DPN-COF(3,K))**COF(4,K)
   10   CONTINUE
        DPP=DPN
      ENDIF
      PROB=0.D0
      U=(CHAN-PRM(1))/PRM(2)
      IF(ABS(PRM(3)).LT.1.D-8) THEN
        IF(ABS(U).LT.10.D0) PROB=EXP(-0.5D0*U**2)
      ELSE
        IF(ABS(PRM(3)*U).LT.50.D0) THEN
          U2=PRM(3)*U-0.5D0*((EXP(PRM(3)*U)-1.D0)/PRM(3))**2
          IF(U2.GT.-50.D0) THEN
            IF(ABS(PRM(3)).GT.0.1D0) THEN
              X=1.D0/(ABS(PRM(3))*SQRT(2.D0))
              T=1.D0/(1.D0+P*X)
              ERFC=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))*EXP(-X**2)
            ELSE
              ERFC=0.D0
            ENDIF
            PROB=EXP(U2)/(1.D0-0.5D0*ERFC)
          ENDIF
        ENDIF
      ENDIF
      PROB=PROB/(SQRT(4.D0*ASIN(1.D0))*PRM(2))
C
      PCHAN=PROB
      RETURN
      END
c      
c      
      double precision  function EZKERX(kiinst, kjcha,  kncha,  diam)
************************************************************************
c
c EZKERX provides the writer of a kernel function subroutine an easy
c method for calling other kernel function subroutines.
c
c INPUTS:
c kiinst - the instrument index of the instrument for which the response
c          is desired.
c kjcha  - the channel of the instrument for which the response is
c          desired.
c kncha  - the number of channels of the instrument for which the
c          response is desired.
c diam   - the diameter for which the instrument response is desired.
c
c LAST MODIFIED: 9 NOV 1987
c
************************************************************************
      integer           iinst,  jcha,   minst,  minstt, mpart
      integer           ncha
      integer           kiinst, kjcha,          kncha
      double precision  diam,   INST1,  INST2,  INST3,  INST4,  INST5
      double precision  INST6,  INST7,  INST8,  INST9,  INST10
      logical           first,  more
      parameter        (minst  = 10)
      common /cmkerf/   iinst,  jcha,   ncha
      save   /cmkerf/
      common /cmparm/   minstt, mpart
      save   /cmparm/
      save              first
      data first       /.true./
c
      if (first) then
        first = .false.
        call PARCHK(minst, minstt, 'minst', 'EZKERX')
      end if
      EZKERX = -1.d00
      if ((iinst .lt. 1)  .or.  (iinst .gt. minst)) return
      call ISWTCH(iinst, kiinst)
      call ISWTCH(jcha,  kjcha)
      call ISWTCH(ncha,  kncha)
      more = .true.
      if (more) then
        more = .false.
        if (iinst .eq. 1) then
          EZKERX = INST1(diam)
        else if (iinst .eq. 2) then
          EZKERX = INST2(diam)
        else if (iinst .eq. 3) then
          EZKERX = INST3(diam)
        else if (iinst .eq. 4) then
          EZKERX = INST4(diam)
        else if (iinst .eq. 5) then
          EZKERX = INST5(diam)
        else
          more = .true.
        end if
      end if
c
c
      if (more) then
        more = .false.
        if (iinst .eq. 6) then
          EZKERX = INST6(diam)
        else if (iinst .eq. 7) then
          EZKERX = INST7(diam)
        else if (iinst .eq. 8) then
          EZKERX = INST8(diam)
        else if (iinst .eq. 9) then
          EZKERX = INST9(diam)
        else if (iinst .eq. 10) then
          EZKERX = INST10(diam)
        end if
      end if
      call ISWTCH(iinst, kiinst)
      call ISWTCH(jcha,  kjcha)
      call ISWTCH(ncha,  kncha)
c
      return
c
*___end of EZKERX___*
c
      end
      
      
      double precision function SENS(d)
************************************************************************
      double precision d
c
      return
c
*___end of SENS___*
c
      end
      
      double precision function XK(x)
************************************************************************
c
      double precision x
       return
      end
      
      double precision function INST2(diam)
      double precision diam
      return
      end



      double precision  function INST3(diam)
************************************************************************
c
      double precision  diam
      return
c
*___end of INST3___*
c
      end
      
      double precision  function INST4(diam)
      return
c
*___end of INST4___*
c
      end


      double precision  function INST5(diam)
************************************************************************
c
      double precision   diam
      return
c
*___end of INST5___*
c
      end


      double precision  function INST6(diam)
      return
c
*___end of INST6___*
c
      end


      double precision  function INST7(diam)
      return
c
*___end of INST7___*
c
      end



      double precision  function INST8(diam)
      return
c
*___end of INST8___*
c
      end


      double precision function INST9(diam)
      return
c
*___end of INST9___*
c
      end

      double precision  function INST10(diam)
************************************************************************
c
      double precision  diam
      return
c
*___end of INST10___*
c
      end




      subroutine MATINT(jrowf, jrowx, mrowf, mrowx, nvec, arryf,
     &                 arryx, fun, x)
************************************************************************
c
c MATINT interpolates a row of a matrix to find the value of a function,
c f, represented as a vector, f(x1), f(x2), f(x3), ... .  The set of
c points x1, x2, x3, ... at which f has been evaluated need not be
c equally spaced.  The interpolation is logrithmic.
c Important: if x is less than x1 then fun = f(x1)
c            if x is greater than xnvec fun(xnvec)
c
c INPUTS:
c mrowf  - the declared row dimension of the array arrayf
c mrowx  - the declared row dimension of the array arrayx
c jrowf  - the row representing the function of interest
c nvec   - the number of values of f(x) which are represented
c          in the array;  these are assumed to be stored in
c          positions 1-nvec in the array arrayf.
c arrayf - the matrix which contains the function values
c arrayx - the matrix which contains the points at which the
c          function is evaluated
c x      - the argument at which the function is to be evaluated
c
c OUTPUTS:
c fun    - the function value at x
c
c LAST MODIFIED: 7 NOV 1987
c
************************************************************************
      integer          ibin,   ihigh,  ilow,   imid,   jrowf,  jrowx
      integer          mrowf,  mrowx,  nvec
      double precision arryf, arryx,  fun,    x,      XLNINT
      dimension        arryf(mrowf, nvec), arryx(mrowx, nvec)
c
c ##FIRST SEE IF WE NEED TO EXTRAPOLATE
c
      if (x .lt. arryx(jrowx, 1)) then
        fun = arryf(jrowf, 1)
      else if (x .gt. arryx(jrowx, nvec)) then
        fun = arryf(jrowf, nvec)
      else
c
c ##ELSE WE NEED TO FIND THE POINTS IN arrayx WHICH JUST BOUND x.
c ##A BINARY SEARCH ROUTINE IS USED
c
        ilow = 1
        ihigh = nvec
        do 300 ibin = 1, nvec + 1
          imid = int((ilow + ihigh)/2)
          if (x .lt. arryx(jrowx, imid)) then
            ihigh = imid
          else
            ilow = imid
          end if
          if ((ihigh - ilow) .eq. 1) go to 1000
300     continue
        write(*, *) ' Binary search routine has failed'
1000    fun = XLNINT(x, arryx(jrowx, ilow), arryx(jrowx, ihigh),
     &        arryf(jrowf, ilow), arryf(jrowf, ihigh))
      end if
c
      return
c
*___end of MATINT___*
c
      end



      double precision function XLNINT(d, dl, du, xl, xu)
************************************************************************
c XLNINT performs logrithmic interpolation.
c
c
c LAST MODIFIED: 7 NOV 1987
c
************************************************************************
      double precision d,   denlog, dl,     du,     rmin
      double precision xl,  xlogi,  xu
      common /cmrmin/  rmin
      save   /cmrmin/
c
      denlog = dlog(du/dl)
      if (dabs(denlog) .lt. rmin) then
        XLNINT = 0.5d00*(xl + xu)
      else
        xlogi = dlog(d/dl)/denlog
        XLNINT = xl + xlogi*(xu - xl)
      end if
c
      return
c
*___end of XLNINT___*
c
      end


      subroutine MVLI(jrow, marry, nfun, arry, fun, x, xl, xu)
************************************************************************
c
c MVLI interpolates a dicretized function whose values are stored in
c one row of a matrix.  The function is assumed to be evaluated at
c equally spaced points, and the interpolation is linear.
c Important: if x lies outside (xl, xu) then fun is set
c equal to the function value at the endpoint.
c
c INPUTS:
c marry  - the declared row dimension of arry
c arry   - the matrix which represent discretized functions
c jrow   - the row representing the function of interest
c nfun   - the number of discrete values given
c xl, xu - the endpoints of the interval over which the function
c          is defined
c x      - the argument at which we'd like the function evaluated
c OUTPUTS:
c fun    - function value at x
c
c LAST MODIFIED: 7 NOV 1987
c
************************************************************************
      integer          jrow,   marry,  mcha,   nfun,   nlow
      double precision arry,   delta,  fun,    rmin,   x,      xl
      double precision xlow,   xu
      common /cmrmin/  rmin
      save   /cmrmin/
      dimension        arry(marry, nfun)
c
c ##FIRST CHECK IF EXTRAPOLATION IS NEEDED
c
      if ((x .lt. xl) .or. (nfun .eq. 1)) then
        fun = arry(jrow, 1)
      else if (x .gt. xu) then
        fun = arry(jrow, nfun)
      else
c
c
        delta = (xu - xl)/dble(nfun - 1)
        nlow = int((x - xl)/(xu - xl)*dble(nfun - 1))
        xlow = xl + dble(nlow)*delta
        call LNINT(xlow, xlow + delta, arry(jrow, nlow + 1),
     &             arry(jrow, nlow + 2), x, fun, rmin)
      end if
c
      return
c
*___end of MVLI___*
c
      end



      double precision function FILTER(d, dup)
************************************************************************
c FILTER models a filter which passes all particles smaller than dup
c but passes no particles larger than 1.5*dup.
c
c LAST MODIFIED:10 NOV 1987
c
************************************************************************
      double precision d,      dup
c
      if (d .lt. dup) then
        FILTER = 1.d00
      else if (d .gt. (1.5d00*dup)) then
        FILTER = 0.d00
      else
        FILTER = 1.d00 - ((d - dup)/(0.5d00*dup))
      end if
c
      return
c
*___end of FILTER___*
c
      end


************************************************************************
*
* usrfmc.f contains the following subroutines
* (modify this list as needed)
* (**) are required by MICRON
*
**LOADR   loads the covariance function for the size distribution
**LSPC    loads user defined linear functionals for information
*           calculations
**INST1   impactor response (mass)
**INST2   impactor response (mass)
**INST3   optical particle counter response (number)
**INST4   cnc response (fitted eq. : number)
**INST5   diffusion battery response (number)
**INST6   perfect filter
**INST7   impactor response (mass)
**INST8   DMA response (number)
**INST9   DMA response (number)
**INST10  DMA response (number)
* MATINT  matrix interpolation
* MVLI    matrix interpolation (log spacing)
* EZKERX  easy INSTxx caller
* SENS    sensitivity function for EAA
* XK      calculates part of the Cunningham correction factor
* DCPART  particle diffusion coefficient
* XLNINT  log interpolation
* FILTER  filter model
* CC      Cunningham correction factor
* MFPAIR  mean free path  of air
* AIRVSC  air viscosity
* FDCHGI  part of DMA response function
* SRF     sensory response fumction for DMA
* DTXF    part of DMA response function
* CHDST   charge sitribution for DMA
* SWV     switching valve response
* DERF    error function
**OUTSPC  specialized output (see programmer manual)
**INPSPC  specialized input  (see programmer manual)
**TF1     lognormal test function
**TF2     power law test function
**TF3     data file test function
* INSTA   EAA response function
* INSTB   Low pressure impactor response function (number)
* INSTC   diffusion battery (collected fraction)
* INSTD   CNC (calibration data)
* INSTE   EAA without electrometer
*
************************************************************************


      subroutine LSPC  (nspc,   ndim,   f)
      integer           nspc,   ndim
      double precision  f(*)
************************************************************************
* LSPC loads user defined linear functionals for information evaluation
* nspc is indexed from 1 to ... by MICRON.  LSPC is not called as soon
* as there is no linear functioal corresponding to nspc.  At this point
* the subroutine should do
* nspc = -nspc
*
* Variables:
* nspc - input : requests the nspc'th linear functional
*      - output: unchanged if there is an nspc'th linear
*                functional; otherwise the nspc = -nspc
*
* LAST MODIFIED:30 OCT 1989
************************************************************************
      nspc = -nspc
      return
*
*___end of LSPC___*
*
      end




      double precision  function TF1(diam)
************************************************************************
*
* TF1 computes the function value for a log-normal distribution which
* represents d(Number)/d(log diam).
*
* VARIABLES:(type, input/output, array)
* isum      (I,I)   - the row in the array tfpar which holds the
*                     parameters for this function.
* msum      (I,I)   - the maximum number of test functions which
*                     can be summed to form the test distribution.
* msumt     (I,I)   - the value of the FORTRAN parameter msum
*                     declared in the main program.
* mtfpar    (I,I)   - the maximum number of function parameters.
* mtfprt    (I,I)   - the value of the FORTRAN parameter mtfpar
*                     declared in the main program.
* diam      (D,I)   - natural log of particle diameter in microns
* gmdiam    (D,I)   - geometric mean diameter, passed
*                     as the second element in the common
*                     variable tfpar
* gstdv     (D,I)   - geometric standard deviation, passed
*                     as the third element in the common
*                     variable tfpar
* tfpar     (D,I,A) - the array containing the parameters for
*                     the test functions
* ttlnum    (D,I)   - total particle concentration #/cc, passed
*                     as the first element in the common
*                     variable tfpar
*
* LAST MODIFIED:13 NOV 1987
************************************************************************
      integer           isum,   isumfn, msum,   msumt,  mtfpar, mtfprt
      double precision  diam,   gstdv,  gmdiam, rt2pi,  tfpar,  ttlnum
      logical           first
      parameter        (rt2pi  = 2.506628275d+0)
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
c
c ## Always save named commons.
c
      common /cmtfpa/   isum,   isumfn(msum), tfpar(msum, mtfpar)
      save   /cmtfpa/
      common /cmfun/    msumt,  mtfprt
      save   /cmfun/
      save              first
      data first       /.true./
*-----------------------------------------------------------------------
c
c ##FIRST CHECK IF THE NUMBER THE DIMENSION OF THE ARRAY tfpar
c ##IS CORRECTLY DECLARED
c
      if (first) then
        first = .false.
        call PARCHK(msum, msumt, 'msum', 'TF1')
        call PARCHK(mtfpar, mtfprt, 'mtfpar', 'TF1')
      end if
c
      ttlnum = tfpar(isum,1)
      gmdiam = tfpar(isum,2)
      gstdv  = tfpar(isum,3)
      TF1 =    ttlnum/rt2pi/dlog(gstdv) *
     &         dexp(-((diam-dlog(gmdiam))**2) /
     &         2.d00/(dlog(gstdv)**2))
      return
*
*___end of TF1___*
*
      end


      double precision  function TF2(diam)
************************************************************************
*
* TF2 computes one term of a power law size distribution function.
*
* INPUTS:
* diam   - natural log of particle diameter in microns.
*
* LAST MODIFIED:13 NOV 1987
************************************************************************
      integer           isum,   isumfn, msum,   msumt,  mtfpar, mtfprt
      double precision  diam,   tfpar
      logical           first
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
      common /cmtfpa/   isum,   isumfn(msum), tfpar(msum, mtfpar)
      save   /cmtfpa/
      common /cmfun/    msumt,  mtfprt
      save   /cmfun/
      save              first
      data first       /.true./
*-----------------------------------------------------------------------
      if (first) then
        first = .false.
        call PARCHK(msum, msumt, 'msum', 'TF2')
        call PARCHK(mtfpar, mtfprt, 'mtfpar', 'TF2')
      end if
c
      TF2 = tfpar(isum,1)*dexp(diam*tfpar(isum,2))
      return
*
*___end of TF2___*
*
      end


      double precision  function TF3(diam)
************************************************************************
*
* TF3 computes one term of a distribution function defined by
* micron.tf3.  micron.tf3 must have x,y format with 2 data per
* record.  A maximum of 1000 points is allowed.
*
* INPUTS:
* diam   - natural log of particle diameter in microns.
*
* LAST MODIFIED:17 JAN 1988
************************************************************************
      integer           ioerr,  iread,  isum,   isumfn, msum,   msumt,
     &                  mtfpar, mtfprt, mwork,  ntf3,   nvec
      double precision  diam,   ediam,  fun,    fx,     tfpar,  x
      logical           first
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
      parameter        (mwork  = 1000)
      dimension         fx(mwork), x(mwork)
      common /cmtfpa/   isum,   isumfn(msum), tfpar(msum, mtfpar)
      save   /cmtfpa/
      common /cmfun/    msumt,  mtfprt
      save   /cmfun/
      save              nvec,   fx,     x
      save              first
      data first       /.true./
*-----------------------------------------------------------------------
      if (first) then
        first = .false.
        call PARCHK(msum, msumt, 'msum', 'TF3')
        call PARCHK(mtfpar, mtfprt, 'mtfpar', 'TF3')
        call FUNIT(ntf3)
        open(ntf3, file='micron.tf3', status = 'old')
        nvec = 0
        do 100 iread = 1, mwork
          read(ntf3, *, err = 1000, end = 1000, iostat = ioerr)
     &        x(iread), fx(iread)
          nvec = nvec + 1
  100   continue
        if (ioerr .gt. 0) write(6,*) ' Reading failure: TF3'
      end if
c
 1000 ediam = dexp(diam)
      call MATINT(1, 1, 1, 1, nvec, fx, x, fun, ediam)
      TF3   = fun
      return
*
*___end of TF3___*
*
      end


      double precision  function TF4(diam)
************************************************************************
*
* TF4 computes one term of a root law size distribution function.
*
* INPUTS:
* diam   - natural log of particle diameter in microns.
*
* LAST MODIFIED:13 NOV 1987
************************************************************************
      integer           isum,   isumfn, msum,   msumt,  mtfpar, mtfprt
      double precision  diam,   tfpar
      double precision  dmax,   scale
      logical           first
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
      common /cmtfpa/   isum,   isumfn(msum), tfpar(msum, mtfpar)
      save   /cmtfpa/
      common /cmfun/    msumt,  mtfprt
      save   /cmfun/
      save              first
      data first       /.true./
*-----------------------------------------------------------------------
      if (first) then
        first = .false.
        call PARCHK(msum, msumt, 'msum', 'TF4')
        call PARCHK(mtfpar, mtfprt, 'mtfpar', 'TF4')
      end if
c
      scale = tfpar(isum,1)
      dmax  = dlog(tfpar(isum,2))
      if (diam .gt. dmax) then
        TF4 = 0.d+0
      else
        TF4 = scale * ( dmax - diam )**(1.d+0 / tfpar(isum,3))
      end if
      return
*
*___end of TF4___*
*
      end
