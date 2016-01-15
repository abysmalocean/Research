      subroutine FNDCNS(fail,   ichpnt, iconrc, id1,    id2,    ierrmd, 
     &                  inspnt, mclin,  mctotl, mdat,   mdim,   mfail,  
     &                  mrecl,  msglev, nclin,  ndim,   nlog,   npdat,  
     &                  opchr,  p,      bl,     bspc,   bu,     c,      
     &                  dspc,   u,      enorm,  fk,     wmach,  y,      
     &                  work,   cfiles)
      integer           fail,   ichpnt, iconrc, id1,    id2,    ierrmd, 
     &                  inspnt, mclin,  mctotl, mdat,   mdim,   mfail,  
     &                  mrecl,  msglev, nclin,  ndim,   nlog,   npdat,  
     &                  opchr,  p
      double precision  bl,     bspc,   bu,     c,      dspc,   u,      
     &                  enorm,  fk,     wmach,  y,      work
      character*(*)     cfiles
      dimension         fail(*), ichpnt(*), iconrc(mclin), inspnt(*), 
     &                  msglev(*), nlog(*), opchr(*), bl(*), bspc(*), 
     &                  bu(*), c(mclin,*), dspc(*), enorm(*), 
     &                  fk(mdat,mdim), u(mdat,*), wmach(*), work(*), 
     &                  y(*), cfiles(*)
************************************************************************
*
* FNDCNS puts the constraints in the lower and upper constraint vectors
* "bl" and "bu", and places the general constraints in the matrix "c".
*
* INPUTS:
* wmach(5) - the smallest positive number on the machine
* wmach(7) - the largest number on the machine
* delta    - (d2 - d1)/(ndim - 1)
*
* OUTPUTS:
* nclin  - the number of general linear constraints
*
* LAST MODIFIED:16 JAN 1988
*
************************************************************************
      integer           iconi,  icono,  idat,   idim,   izero
      double precision  DDOT,   fmax,   one,    rmin,   rmax,   ydot,   
     &                  ysum,   xbig,   zero
      logical           test
      parameter        (one  = 1.d+0)
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      fmax  = wmach(8) / 100.d+0
      rmin  = wmach(5)
      rmax  = wmach(7)
      xbig  = rmax/2.d+0
      nclin = 0
      test  = .true.
c
c ##FIRST INITIALIZE bl AND bu ELEMENTS ASSOCIATED WITH f(x)
c ##TO HELP WITH NUMERICAL STABILITY, ALL SOLUTIONS ARE BOUNDED ABOVE
c
      call DLOAD(ndim, zero, bl, 1)
      call DLOAD(ndim, fmax,  bu, 1)
      if ((bspc(1) .gt. -rmin) .and. (bspc(1) .le. xbig)) then
        bl(1)    = bspc(1)
        bu(1)    = bspc(1)
      end if
      if ((bspc(2) .gt. -rmin) .and. (bspc(2) .le. xbig)) then
        bl(ndim) = bspc(2)
        bu(ndim) = bspc(2)
      end if
c
c ##IF THE INVERSION INTERVAL CHOSEN LIES INSIDE THE KERNEL FUNCTION
c ##SUPPORT, THEN WE NEED TO CONSTRAIN THE ENDS OF THE SOLUTION VECTOR
c ##TO ZERO
c
      if (id1 .gt. 0) then
        bl(1)    = zero
        bu(1)    = zero
      end if
      if (id2 .gt. 0) then
        bl(ndim) = zero
        bu(ndim) = zero
      end if
c
c ##NOW FIND THE CONSTRAINTS ON THE INVERTED DATA
c
      call ILOAD(npdat, 0, iconrc, 1)
      do 100 idat = mdim + 1, mdim + npdat
        iconi = idat - mdim
        if ((bl(idat) .ge. -rmin) .or. (bu(idat) .ge. -rmin)) then
c
c ##  THEN WE HAVE A GENERAL CONSTRAINT
c
          nclin = nclin + 1
          if (nclin .le. mclin) then
            icono = ndim + nclin
            bl(icono) = bl(idat)
            bu(icono) = bu(idat)
            call DCOPY(ndim, fk(iconi, 1), mdat, c(nclin, 1), mclin)
c
c ##  KEEP TRACK OF WHICH DATA IS BOUNDED
c
            if (bl(icono) .ge. -rmin) then
              iconrc(iconi) = 1
            else
              bl(icono) = -rmax
            end if
            if (bu(icono) .ge. -rmin) then
              iconrc(iconi) = iconrc(iconi) + 2
            else
              bu(icono) = rmax
            end if
          end if
        end if
100   continue 
c
c ##ADD THE CONSTRAINTS IMPOSED BY THE ERROR ARRAY.  THESE ARE THE
c ##RESULT OF A ROW OF ZEROS IN THE ERROR ARRAY OR AN EIGENVALUE
c ##OF THE GENERALIZED INVERSE EQUAL TO ZERO.  THE VECTORS THAT
c ##REPRESENT THE CONSTRAINTS ARE STORED IN THE TRAILING ROWS OF 
c ## u 
c
      do 200 izero = 1, npdat - p
        nclin = nclin + 1
        if (nclin .le. mclin) then
          icono = nclin + ndim
          ydot = DDOT(npdat, y, 1, u(1,npdat+1-izero), 1)
          bl(icono) = ydot
          bu(icono) = ydot
          do 210 idim = 1, ndim
            c(nclin,idim) = 
     &        DDOT(npdat, u(1,npdat+1-izero), 1, fk(1,idim), 1)
210       continue
        end if
200   continue
c
c ##CHECK FOR VOLUME CONSTRAINTS
c
c
c ##CHECK FOR NUMBER CONSTRAINTS
c
      if ((bspc(3) .ge. -rmin) .or. (bspc(4) .ge. -rmin)) then
        nclin = nclin + 1
        if (nclin .le. mclin) then
          icono = nclin + ndim
          call DLOAD(ndim-2, dspc(7), c(nclin,2), mclin)
          c(nclin, 1)    = dspc(7) * 0.5d00
          c(nclin, ndim) = c(nclin, 1)
          bl(icono) = bspc(3)
          bu(icono) = bspc(4)
          if      (bspc(3) .lt. -rmin) then
            bl(icono) = -rmax
          else if (bspc(4) .le. -rmin) then
            bu(nclin) =  rmax
          end if
        end if
      end if
c
c ##BOUND THE TOTAL INVERTED
c ##ERROR.  RESCALING PREVENTS THE CONSTRAINT FROM BECOMING UNREALISTIC.
c
      if (opchr(2) .gt. 0) then
c
c ##FIRST FIND THE INTEGRAL OF THE KERNEL FUNCTIONS, STORED IN
c ##work(npdat+1), AND THE INVERSE STORED IN work.
c
        call RDMEAN(fail,   ichpnt, inspnt, mdat, mfail,  mrecl,  
     &              msglev(2), nlog(2), work(npdat+1), cfiles, 'INTG')
        call DLOAD(npdat, one, work, 1)
        call DDDIV(npdat, work(npdat+1), 1, work, 1)
        nclin = nclin + 1
        if (nclin .le. mclin) then
          icono = nclin + ndim
          ysum = DDOT(npdat, work, 1, y, 1)
          bl(icono) = ysum 
          bu(icono) = rmax 
          do 300 idim = 1, ndim
            c(nclin, idim) = DDOT(npdat, work, 1, fk(1,idim), 1) 
300       continue
        end if
      end if
c
1000  if (nclin .gt. mclin)
     &  call ERROR(93, nlog(2), msglev(2), nclin, 0.d+0, ' ', fail(1))
      return
c
*___end of FNDCNS___*
c 
      end
      subroutine SIGTXF(fail,   ierrmd, iwork,  kerr,   mdat,   merr,   
     &                  msglev, ndim,   npdat,  nlog,   opchr,  p,      
     &                  e,      einv,   enorm,  fk,     wmach,  work,   
     &                  y) 
      integer           fail,   ierrmd, iwork,  kerr,   mdat,   merr,   
     &                  msglev, ndim,   npdat,  nlog,   opchr,  p
      double precision  e,      einv,   enorm,  fk,     wmach,  work,  
     &                  y
      dimension         fail(*), iwork(*), msglev(*), nlog(*), opchr(*), 
     &                  e(mdat,*), einv(mdat,*), enorm(*), fk(mdat,*), 
     &                  wmach(*), work(*), y(*)
************************************************************************
*
* SIGTXF find the 'inverse' of the error array, "e".  The inversion for
* the diagonal case is treated separately.  This subroutine also finds
* the constaints that are imposed by the error matrix and these are
* stored in the array "e".
*
* VARIABLES:(type,input/output,array)
* work     :(D,S,A)  is a scratch array that provides workspace for
*                    calcualting the inverse in the non-diagonal case.
*                    work is divided up in the following manner:
*                    iw - marks the beginning of the array of length
*                         merr that is used as a scratch array for the
*                         subroutine LINSVD.  
*                    isr- the initial index of the array of length
*                         mdat that holds the numbers used to rescale
*                         each row of the error array
*                    is - the initial index of the array used to hold
*                         the eigenvalues of the error array
* 
* LAST MODIFIED:19 JAN 1989
************************************************************************
      integer           iadd,   IDAMAX, idiag,  IDRANK, ifk,    imax,   
     &                  inform, ipos,   is,     iscale, isort,  isr,    
     &                  iw,     idat,   izcol,  izero,  jdat,   jerr,   
     &                  msglog, nlog2,  nnonz,  nrepl,  nsort,  nzeig,  
     &                  nzero
      double precision  DNRM2,  eps,    esmall, ksmall, one,    rmin,   
     &                  scale,  ysmall, zero
      parameter        (zero = 0.d+0)
      parameter        (one  = 1.d+0)
*-----------------------------------------------------------------------
      if ((fail(1) .gt. 0) .or. (opchr(9) .eq. 0)) return
      rmin  = wmach(5)
      eps   = wmach(3)
c
c ##MARK THE WORKSPACE POINTERS
c
      iw  = 0
      isr = iw  + merr     
      is  = isr + mdat 
      msglog = msglev(2)
      nlog2  = nlog(2)
      nzero  = 0
      call ILOAD(npdat, 0, iwork, 1)
      if (ierrmd  .eq. 1) then
        p = npdat
c
c ##THEN THE ERROR ARRAY IS DIAGONAL, AND THE CALCULATION IS TRIVIAL
c
        do 100 idiag = 1, npdat
c
c ##FIRST CHECK IF THE ERROR IS NEGLIGIBLE.  THE ERROR IS NEGLIGIBLE
c ##IF IT IS NEARLY ZERO, SMALL COMPARED TO THE DATA OR SMALL COMPARED
c ##TO THE MAXIMUM OF THE KERNEL FUNCTION
c
          ifk    = IDAMAX(ndim, fk(idiag,1), mdat)
          ksmall = dabs(fk(idiag,ifk)) * wmach(6) * 100.d+0
          ysmall = dabs(y(idat))       * wmach(6) * 100.d+0
          esmall = dmax1(rmin, ysmall, ksmall)
          if (e(idiag,idiag) .le. esmall) then
            p                    = p  - 1
            nzero                = nzero + 1
            iwork(nzero)         = idiag
            einv(idiag,idiag)    = zero
          else
            einv(idiag,idiag)    = one / e(idiag,idiag)
          end if
100     continue
        if (nzero .ne. 0) call WARNG(21, nlog2, msglog, nzero,0.d+0,' ')
        call DCOPY(npdat, einv, mdat+1, enorm, 1)
c
c ##KEEP TRACK OF THE CONSTRAINTS
c
        do 150 izero = 1, nzero
          call DLOAD(npdat, zero, e(1,npdat-izero+1), 1)
          e(iwork(izero),npdat-izero+1) = one
150     continue
      else
c
c ##ELSE WE'LL NEED TO PERFORM AN SVD.
c
        call SCRNM(24, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
c ##FIRST RESCALE E SO SMALL EIGENVALUES HAVE MEANING.
c ##THE SCALING THAT MAKES THE LARGEST ELEMENT IN EACH ROW 
c ##EQUAL TO ONE IS A REASONABLE CHOICE. IF THE LARGEST ELEMENT
c ##IN THE ROW IS "ZERO" THEN THE DATUM BECOMES A CONSTRAINT,
c ##AND THIS ROW IS REMOVED FROM THE ERROR MATRIX.
c
        nrepl = 0
        do 200 idat = 1, npdat
          imax = IDAMAX( kerr, e(idat,1), mdat )
          if (dabs(e(idat,imax)) .lt. rmin) then
c
c ##REMOVE THE ROW FROM THE ERROR ARRAY
c
            nrepl          = nrepl + 1
            iwork(nrepl)   = idat
            enorm(idat)    = zero
            work(isr+idat) = zero
          else
c
c ##ELSE RESCALE THE ROW
c
            work(isr+idat) = one/e(idat,imax)
            call DSCAL (kerr, work(isr+idat), e(idat,1), mdat)
            if (nrepl .gt. 0)
     &        call DCOPY(kerr, e(idat,1), mdat, e(idat-nrepl,1),mdat)
            enorm(idat) = one/DNRM2(kerr, e(idat,1), mdat)
          end if
200     continue
        if (nrepl .gt. 0) call WARNG(21,nlog2,msglog, nrepl, 0.d+0, ' ')
        nnonz = npdat - nrepl
c
c ##REMOVE COLUMNS FROM e THAT CONTAIN ONLY ZERO'S
c
        izcol = 0
        do 300 izero = 1, kerr
          imax = IDAMAX(npdat, e(1,izero), 1)
          if (dabs(e(imax,izero)) .le. eps) then
            call  WARNG(23, nlog2, msglog, izero, 0.d+0, ' ')
            izcol = izcol + 1
          else if (izcol .ne. 0) then
            call DCOPY (npdat, e(1,izero), 1, e(1,izero-izcol), 1)
          end if
300     continue
        call EISSVD(mdat, nnonz, kerr-izcol, e, work(is+1), 
     &              .true., einv, .false., e, inform, work(iw+1))
        if (inform .ne. 0) then
          call ERROR(37, nlog2, msglog, inform, 0.d+0, ' ', fail(1))
c
c ##DECIDE WHICH OF THE EIGENVALUES ARE ZERO
c
        else if (dabs(work(is+1)) .le. rmin) then
          call ERROR(34, nlog2, msglog, 0, 0.d+0, ' ', fail(1))
        end if
        if (fail(1) .gt. 0) return
c
c ##REORDER THE EIGENVALUES AND SEE HOW MANY ARE NONZERO
c
        do 400 isort = 1, nnonz
          nsort = kerr - izcol + 1 - isort
          iadd  = IDAMAX(nsort, work(is+isort), 1) - 1
          if (iadd .ne. 0) then
            call DSWAP(1, work(isort+is), 1, work(isort+iadd+is), 1)
            call DSWAP(nnonz, einv(1,isort),1, einv(1,isort+iadd),1)
          end if
400     continue
        p     = IDRANK(nnonz, work(is+1), 1, eps)
        nzeig = nnonz  - p
        if (nzeig .ne. 0) 
     &    call WARNG(20, nlog2, msglog, nzeig, 0.d+0, ' ')
c
c ##RESCALE THE LEFT SINGULAR EIGENVECTORS BY THE 
c ##EIGENVALUES
c
        nzero  = nrepl + nzeig
        do 500 iscale = 1, p   
          scale  = one/work(is+iscale)
          call DSCAL(nnonz, scale, einv(1,iscale), 1)
500     continue
c
c ##TAKE THE TRANSPOSE
c
        call TXPOSE(p, mdat, einv)
c
c ##REPLACE THE ZERO COLUMNS IN THE INVERSE, AND ADD THE CORRESPONDING
c ##CONSTRAINTS TO THE MATRIX CONTAINING THE LEFT EIGENVECTORS.
c
        do 600 ipos = npdat, 1, -1
          if (nrepl .ne. 0) then
            if (iwork(nrepl) .eq. ipos) then
              call DLOAD(p, zero, e(1,nnonz+nrepl), 1)
              e(ipos,nnonz+nrepl) = one
              nrepl = nrepl - 1
            else
              call DCOPY(p, einv(1,ipos-nrepl), 1, einv(1,ipos),1)
            end if
          end if
600     continue
c
c ##UNDO THE SCALING WHICH HELPED US DECIDE WHICH EIGENVALUES ARE
c ##SMALL
c
        do 700 iscale = 1, npdat
          call DSCAL(p, work(isr+iscale), einv(1,iscale), 1)
700     continue
      end if
c
c ##OUTPUT INVERSE OF ERROR ARRAY
c
      if (msglog .ge. 70) then
        if (ierrmd .eq. 1) then
          write(nlog2,3)
          call NZPRN(npdat, rmin, einv, mdat+1, nlog2)
        else 
          write(nlog2, 1) p, npdat
          do 900 jdat = 1, p
            write(nlog2,2) jdat
            write(nlog2,4) (einv(jdat,jerr), jerr = 1, npdat)
900       continue
        end if 
      end if
      return
1     format(//,' Inverse of error array:',/,
     &      ' The inverse of the error array is ',i4,' by ',i4)
2     format(/,' Row ',i3, ' of inverse error array . . .')
3     format(//,' The diagonal elements of the inverse of the error',/,
     &' array:',/)
4     format(6(2x,d11.4))
*
*___end of SIGTXF___*
*
      end
      subroutine TGTNRM(fail,   mbrent, mintg,  msglev, nlog,   opchr,  
     &                  p,      obj,    sqemin, tnrmsq, tol,    wmach,  
     &                  yhatm)
      integer           fail,   mbrent, mintg,  msglev, nlog,   opchr,  
     &                  p
      double precision  obj,    sqemin, tnrmsq, tol,    wmach,  yhatm  
      dimension         msglev(*), nlog(*), opchr(*), tol(*), wmach(*), 
     &                  yhatm(*), fail(*)
************************************************************************
*
* TGTNRM computes "tnrmsq" which is the target value of 
* (1/p) || Yhat^M - Yhat^R ||^2 where Yhat^M denotes the measured data 
* and Yhat^R is the recovered data.  "tnrmsq" is found by computing
* (1/p) E(||Yhat^T - Yhat^M||^2), where E denotes expectation operator 
* and Yhat^T is the true data.  TGTNRM also determines the accuracy with 
* which MICRON solves (1/p) || Yhat^M - Yhat^R ||^2 = "tnrmsq". This is 
* found by computing the square root of the variance of
* "tnrmsq" - (1/p) || Yhat^T - Yhat^M ||^2.
*
* OUTPUTS:
* sqemin  - the minimum value of (1/p) ||Y^T - Y^M||^2.
* tnrmsq  - the target value of (1/p) ||Y^R - Y^M||^2.  
* tol(2)  - the accuracy with which we attempt to achieve tnrmsq
*
* LAST MODIFIED:21 JAN 1989
************************************************************************
      integer           ierbrn, iquit
      double precision  argln,  b,      DBRENT, denom,  DLNGAM, DNRM2,  
     &                  dp,     dphalf, dsctol, dsqe,   ERRM1,  garb,
     &                  one,    posml,  rmin,   tlu,    tmptar, tmptol, 
     &                  tolrt,  tuu,    two,    UPFUN,  upper,  VARCOR, 
     &                  ynrmsq, zcheck, zero
      character*4       test 
      parameter        (zero   = 0.d+0)
      parameter        (one    = 1.d+0)
      parameter        (two    = 2.d+0)
      parameter        (test   = 'TEST')
      parameter        (dsctol = 1.d+0)
      external          ERRM1, UPFUN, VARCOR
*-----------------------------------------------------------------------
      if ((fail(1) .gt. 0) .or. (opchr(4) .le. 0)) return
      rmin    = wmach(5)
      dp      = dble(p)
      posml   = dlog(rmin)
      dphalf  = dp/two
c
c ##FIRST COMPUTE THE MINIMUM ||Yhat^M - Yhat^T||^2
c
      ynrmsq = DNRM2(p, yhatm, 1)**2
      sqemin = two * obj 
      if (sqemin .lt. zero) then
        call WARNG(24, nlog(2), msglev(2), 0, 0.d+0, ' ')
        sqemin = zero
      else if (sqemin .gt. dp) then
        call WARNG(33, nlog(2), msglev(2), 0, 0.d+0, ' ')
      end if
c
c ##IF THE ERROR BALL LIES NEXT TO AN INFINITE PLANE
c
      if (test .eq. 'TEST') then
        tmptar = zero
        tnrmsq = dsctol*dp + sqemin
      else
c
c ##ELSE THE ERROR BALL LIES INSIDE CONCENTRIC SPHERES
c
      dsqe   = dsqrt(sqemin)
      if (ynrmsq .gt. rmin) then
        zcheck = sqemin/ynrmsq
      else
        zcheck = one
      end if
c
c ##IF sqemin IS ZERO A SIMPLE EXPRESSION EXIST FOR "tnrmsq" AND
c ##tol(2).  IF NOT WE MUST COMPUTE tmptar, THE CORRECTION TO "tnrmsq"
c
      if ((sqemin .le. rmin) .or. (zcheck .le. rmin)) then
        tmptar = zero
      else if (sqemin .le. one) then
        if (p .lt. 40) then
          call DINTG(ERRM1, iquit, mintg, denom, 0.d00, dsqe, tol(6), 
     &               garb, rmin)
          denom  = dexp(DLNGAM(dphalf)) * 
     &             (two**(dphalf - one)) - denom
          tmptar = dexp(dlog(dsqe) * dp + sqemin/two) / denom
        else
c
c ##FOR LARGE p we use a limiting formula
c
          argln  = dlog(sqemin)*dphalf - sqemin/two -
     &             DLNGAM(dphalf) - (dphalf - one)*dlog(two)
          tmptar = dexp(argln)
        end if
      else
c
c ##HERE WE COMPUTE TMPTAR FOR LARGER SQEMIN BY EVALUATING AN
c ##IMPROPER INTEGRAL.  UPPER, REPLACES INFINITY AS THE UPPER LIMIT
c ##OF THE INTEGRAL, AND IS FOUND BY DETERMINING FOR WHICH VALUE OF X
c ##VARCOR(X) = rmin.  THIS IS DETERMINED BY SOLVING A NONLINEAR EQUATION
c ##WITH DBRENT.
c
c ##FIRST BOUND THE ROOT WITH tlu AND tuu, AND DETERMINE A REASONABLE TOLERANCE 
c ##FOR FINDING THE ROOT
c
        tlu    = dsqrt(one - two*posml/sqemin)
        b      = (dp - one)/sqemin
        tuu    = b + dsqrt(b**2 + tlu**2)
        tolrt  = (tlu - one)/10.d+0
c
c
        upper  = DBRENT(ierbrn, UPFUN, tlu, tuu, tolrt, mbrent, wmach)
        upper  = upper + tolrt
        call DINTG(VARCOR, iquit, mintg, denom, one, upper, tol(6), 
     &             garb, rmin)
        tmptar = one/denom
      end if
      tnrmsq   = tmptar + dp*dsctol
c
c ##COMPUTE THE SQUARE ROOT OF THE VARIANCE OF TNRMSQ FOR THE
c ##METHOD OF DISCREPANCY
c
      end if
      if (iabs(opchr(1)) .eq. 4) then
        tmptol   = two*dp + (sqemin - dp + two - tmptar)*tmptar
        tmptol   = dsqrt(tmptol)
        if (opchr(3) .gt. 0) then
          tol(2) = tmptol
          if (msglev(2) .ge. 40) write(nlog(2), 3) tol(2)
        else if (msglev(2) .ge. 30) then
          write(nlog(2),4) tol(2), tmptol
        end if
      else if (tol(2) .lt. wmach(3)) then
        tol(2) = 0.1d+0
      end if
      if (msglev(2) .ge. 40) then
        write(nlog(2), 1) sqemin
        if (iabs(opchr(1)) .eq. 4) write(nlog(2), 2) tnrmsq
      end if
c
      return
1     format(/,' The minimum norm square of the transformed  error',/,
     &' vector is ',d11.4)
2     format(/,' The expected norm square of the transformed error',/,
     &' vector is ',d11.4)
3     format(/,' The allowed  difference between  the actual  norm',/,
     &         ' square  of  the inverted  error and the target is',/,
     &         ' (tol(2)) ',d11.4)
4     format(/,
     &' The value specified for tol(2), (', d10.3,'), will',/,
     &' be used instead of the computed value (',d10.3,').')
*
*___end of TGTNRM___*
*
      end
      double precision function UPFUN(y)
************************************************************************
*
* We find the zero of this function which is greater than 1 to determine
* a reasonble upper bound on the indefinite integral used to compute
* "tmptar" in TGTNRM.
*
* LAST MODIFIED:12 NOV 1987
************************************************************************
      integer          p
      double precision rmin,    sqemin, y
      common /cmrmin/  rmin
      save   /cmrmin/
      common /cmp/     p 
      save   /cmp/
      common /cmsqem/  sqemin
      save   /cmsqem/
c
      UPFUN = dlog(rmin) - dble(p - 1)*dlog(y) - 
     &        sqemin*(1.d00 - y**2)/2.d00
c
      return
*
*___end of UPFUN___*
*
      end
      double precision function ERRM1(y)
************************************************************************
*
* ERRM1 is integrated to determine a term in the expression for
* "tmptar" in TGTNRM for the case: sqemin .lt. 1, p less than 40
*
* LAST MODIFIED:12 NOV 1987
************************************************************************
      integer          p
      double precision y
      common /cmp/     p 
      save   /cmp/
c
      ERRM1 = dexp((-y**2)/2.d00)*(y**(p - 1))
c
      return
*
*___end of ERRM1___*
*
      end
      double precision function VARCOR(y)
************************************************************************
*
* VARCOR is integrated to determine the a term in the expression for 
* "tmptar" in TGTNRM for the case sqemin .gt. 1.
*
* LAST MODIFIED:12 NOV 1987
************************************************************************
      integer          p
      double precision y,      sqemin, exarg
      common /cmp/     p
      save   /cmp/
      common /cmsqem/  sqemin
      save   /cmsqem/
c
      exarg = sqemin * (1.d00 - y**2) / 2.d00 + dble(p - 1) *
     &        dlog(y)
c
      VARCOR = dexp(exarg)
      return
*
*___end of VARCOR___*
*
      end
