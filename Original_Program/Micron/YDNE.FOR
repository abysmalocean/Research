      subroutine NRMDIF(fail,   iter,   mdat,   mdim,   msglev, ndim,   
     &                  nlog,   opchr,  p,      bl,     drvhst, dspc,   
     &                  f,      fkhat,  obj,    objvec, solspc, smp,    
     &                  tnrmsq, wmach,  work,   yhatm,  yhatr,  smooth)
      integer           fail,   iter,   mdat,   mdim,   msglev, ndim,   
     &                  nlog,   opchr,  p
      double precision  bl,     drvhst, dspc,   f,      fkhat,  obj,    
     &                  objvec, solspc, smp,    tnrmsq, wmach,  work,   
     &                  yhatm,  yhatr
      character*(*)     smooth
      dimension         fail(*), msglev(*), nlog(*), opchr(*), bl(*),
     &                  drvhst(0:*), dspc(*), f(*), fkhat(mdat,*), 
     &                  objvec(0:*), solspc(*), wmach(*), work(*), 
     &                  yhatr(*), yhatm(*), smooth(*) 
************************************************************************
*
* NRMDIF computes dsd1 = tnrmsq  -  || Yhat^R - Yhat^M ||^2 
* where tnrmsq is the target of || Yhat^R - Yhat^M || computed in 
* TGTNRM and Yhat^M is the transformed measured data and Yhat^R is the 
* recovered transformed data.
*
* INPUTS:
* tnrmsq - 
*
* OUTPUTS:
* yhatr  -
* dsd1   - 
*
* LAST MODIFIED: 2 JAN 1988
************************************************************************
      integer           i
      double precision  DNRM2,   DNRMDF, DROUGH, yhmsq
      save              yhmsq
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (iter .eq. 1) yhmsq = (DNRM2(p, yhatm, 1))**2
c
c ##CLEAN UP SLOPPY REPORTING
c
      do 100 i = 1, ndim
        if (f(i) .lt. bl(i)) f(i) = bl(i)
100   continue
c
c
      if ((opchr(8) .gt. 0) .and. (iter .eq. 0)) then
        objvec(iter) = 2.d+0 * obj
      else
        objvec(iter) = dabs(2.d+0 * obj + yhmsq) / dble(p)
      end if
      call MMMLT(mdat, mdim, mdat, ndim, 1,  p, fkhat, f, yhatr)
      solspc(3) = DROUGH(ndim, dspc(7), f, work)
      if (iabs(opchr(1)) .eq. 4) drvhst(iter) = 
     &  (DNRMDF( p, yhatr, 1, yhatm, 1, work))**2 - tnrmsq
c
      if (msglev(5) .ge. 50) write(nlog(5),2) iter, smp, drvhst(iter)
      if (msglev(2) .ge. 40)
     &  write(nlog(2), 3) iter, smp, smooth(iabs(opchr(1))), 
     &                    drvhst(iter), solspc(3), objvec(iter)
c
2     format(
     &'     iter = ',i3,',',2x,' smp = ',d11.3,',',2x,' smcf = 'd11.3) 
3     format(//,' Smoothing parameter iteration: ',i3,/,
     &' smoothing parameter', t41, '= ',d14.7,/,
     &' Smoothing choice function (',a4,')', t41, '= ',d14.7,/,
     &' roughness of solution', t41, '= ',d14.7,/,
     &' objective function', t41, '= ',d14.7)
      return
*
*___end of NRMDIF___*
*
      end
      subroutine DONED (fail,   inform, iter,   itersv, itnsmp, mdim,   
     &                  miter,  msglev, ndim,   nlog,   drvhst, dsd,    
     &                  f,      f2,     sp,     spdsda, sphst,  spinfb, 
     &                  sqemin, tnrmsq, tol,    wmach,  done,   bnd)
      integer           fail,   inform, iter,   itersv, itnsmp, mdim,   
     &                  miter,  msglev, ndim,   nlog
      double precision  drvhst, dsd,    f,      f2,     sp,     spdsda, 
     &                  sphst,  spinfb, sqemin, tnrmsq, tol,    wmach
      logical           done,   bnd
      dimension         fail(*), msglev(*), nlog(*), drvhst(0:*), 
     &                  dsd(3), f(*), f2(*), sp(3), sphst(0:*), tol(*), 
     &                  wmach(*)
************************************************************************
*
* DONECH determines if the search for the proper smoothing parameter
* should be terminated.  The search is terminated if
* 1) any previous subroutine has failed
* 2) the requested number of smoothing iterations have been
*    performed
* 3) the discrepancy difference is not changing
* 4) the smoothing parameter is not changing
* 5) the number of smoothing parameters we've tried is greater than
*    miter
* 6) the discrepacy difference is within the specified
*    tolerance, tol(2).
* 7) unconstrained GCV or ridge regression is used.
*
* If a poor guess of the smoothing parameter is made, dsd(1) will
* be nearly equal to it's minimum value, and in this case it is 
* possible dsd(1) will be insensitive to the smoothing parameter.
* tol(7) is used to check if dsd(1) is too close to it's minimum
* (relative to tnrmsq) to terminate the search.
*
* LAST MODIFIED:25 JAN 1989
************************************************************************
      integer           msglog, nlog2
      double precision  diff,   dds1,   dds2,   ds1,    ds2,    dsdf,   
     &                  dsdmin, dtol,   one,    pct10,  zero
      parameter        (one   = 1.d+0)
      parameter        (zero  = 0.d+0)
      parameter        (pct10 = 0.10d+0)
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
c
c ##SOME DISCREPANCY INITIALIZATION
c
      msglog = msglev(2)
      nlog2  = nlog(2)
      done   = .false. 
      dsd(1) = drvhst(iter) 
      dtol   = wmach(5)
      itersv = iter
      ds1    = dsd(1)
      dds1   = dabs(ds1)
c
c ##FIRST DETERMINE IF THE PROPER SMOOTHING PARAMETER HAS BEEN BOUNDED
c
      if (iter .eq. 1) 
     &  call DCOPY(2, dsd(1), 0, dsd(2), 1)
      bnd = .not.((dsd(1)*dsd(2) .gt. zero) .and.
     &           (dsd(2)*dsd(3) .gt. zero))
      if (.not.((iter .eq. 1) .or. (dds1 .lt. tol(2)))) then
c
c ##NOW CHECK FOR COMPUTATIONAL DIFFICULTIES
c
        ds2  = drvhst(iter - 1)
        dds2 = dabs(ds2)
        if ((sp(1)*sphst(iter - 1)) .le. dtol) then
          diff = one
        else
          diff = dabs(dlog(sp(1)/sphst(iter - 1)))
        end if
        dsdf   = dabs(ds1 - ds2)/dmin1(dds1,dds2)
        dsdmin = (sqemin + dsd(1) + tnrmsq)/tnrmsq
        if (diff .lt. 0.5d+0*tol(3)) then
c
c ##  THEN THE DISCREPANCY DIFFERENCE IS TOO SENSITIVE TO THE SMOOTHING
c ##  PARAMETER
c
	  call ERROR(24, nlog2, msglog, 0, tol(3), ' ', fail(15))
        else if ((dsdf .lt. tol(1)) .and. (bnd .or. (iter .ge. 4)) 
     &            .and. (dsdmin .gt. tol(7))) then
c
c ##  THEN THE DISCREPANCY DIFFERENCE IS INSENSITIVE TO THE SMOOTHING
c ##  PARAMETER
c
	  call ERROR(25, nlog2, msglog, 0, 0.d00, ' ', fail(16))
        else if (iter .eq. miter) then
	  call ERROR(26, nlog2, msglog, miter, 0.d00, ' ', fail(17))
        end if
c
c ##CHECK IF THE DISCREPANCY DIFFERENCE IS INCREASING AS THE 
c ##SMOOTHING PARAMETER DECREASES.  IF SO WE NEED TO STOP AND
c ##MAKE SURE WE KEEP THE PREVIOUS f IN f2
c
        if ((sp(1) .lt. sphst(iter - 1)) .and. (ds1 .gt. ds2) .and.
     &     (dsdmin .gt. tol(7))) then
          itersv = itersv - 1
          call ERROR(31, nlog2, msglog, 0, 0.d00,  ' ', fail(18))
        end if
      end if
      if (fail(18) .le. 0) call DCOPY(ndim, f, 1, f2, 1)
c
c ##CHECK IF THE MAXIMUM NUMBER OF ITERATIONS HAS BEEN REACHED AND
c ##WRITE WARNING MESSAGES
c
      if ((dds1 .gt. tol(2)) .and. (itnsmp .eq. iter))
     &  call WARNG(9, nlog2, msglog, itnsmp, 0.d+0, ' ')
c
c 
C     if ((spdsda .gt. zero) .and.
C    &    (inform .eq. 0)    .and.  (dlog(sp(1)/spinfb) .le. pct10)) 
C    &      call ERROR(113, nlog2, msglog, 0, 0.d+0, ' ', fail(18))
c
c ##SEE IF THIS INVERSION IS DONE DUE TO COMPUTATIONAL DIFFICULTIES
c
      if ((fail(15) .gt. 0) .or. (fail(16) .gt. 0) .or. 
     &    (fail(17) .gt. 0) .or. (fail(18) .gt. 0))
     &    call UINDEX(fail(1))
C     done = done .or. 
C    &       ((dds1 .lt. tol(2)) .and. (inform .ne. 1))
      done = done .or. (dds1 .lt. tol(2))
      return
*
*___end of DONECH___*
*
      end
      subroutine FSPD  (fail,   inform, iter,   kfiles, mdim,   msglev, 
     &                  ndim,   nlog,   opchr,  p,      a,      dsd,    
     &                  dspc,   sp,     spdsda, spinfb, tnrmsq, tol,    
     &                  wmach,  bnd)
      integer           fail,   inform, iter,   kfiles, mdim,   msglev, 
     &                  ndim,   nlog,   opchr,  p
      double precision  a,      dsd,    dspc,   sp,     spdsda, 
     &                  spinfb, tnrmsq, tol,    wmach
      logical           bnd
      dimension         fail(*), msglev(*), nlog(*), opchr(*), 
     &                  a(mdim,*), dsd(3), dspc(*), sp(3), tol(*), 
     &                  wmach(*)
************************************************************************
*
* FSPD tries to find the smoothing parameter that makes the 
* discrepancy difference equal to zero.  
*
* INPUTS:
* sp(1)  - most recent smoothing parameter 
*
* OUTPUTS:
* sp(1)  - new smoothing parameter 
* 
* LAST MODIFIED:19 DEC 1987
* bug:This subroutine is too complicated.
************************************************************************
      integer           IDAMAX, IDAMIN, imax,   imin,   info
      double precision  d,      DDIV,   denom,  dsp1,   dsp2,   
     &                  dsp3,   dsda,   dsdb,   e,      epsinv, one,    
     &                  rmin,   scale,  sclrg,  
     &                  spdsdb, spinfa, spmax,  spmin,  temp,   
     &                  ten,    tol1,   zero
      logical           dfail,  initl
      parameter        (one   = 1.d+0)
      parameter        (ten   = 10.d+0)
      parameter        (zero  = 0.d+0)
      parameter        (sclrg = 3.d+0)
      save              e,      d,      dsda,   dsdb,   spdsdb, 
     &                  spinfa, spmax,  spmin,  initl
*-----------------------------------------------------------------------
      info = 0
      if (fail(1)  .gt. 0) return
      if (iabs(opchr(1)) .ne. 4) return
      epsinv = one / wmach(3)
      rmin   = wmach(5)
c
c ##ON THE FIRST ITERATION WE ONLY NEED TO DETERMINE A REASONABLE
c ##VALUE FOR THE SMOOTHING PARAMETER.  
c
      if ((iter .eq. 1)) then
        denom = dble(p) * 6.d+0 / dspc(7)**3
c
c ##COMPUTE THE MINIMUM AND MAXIMUM REASONABLE SMOOTHING PARAMETERS
c
        imax  = IDAMAX(ndim, a, mdim + 1)
        imin  = IDAMIN(ndim, a, mdim + 1)
        temp  = DDIV(a(imax,imax), wmach(3), dfail)
        spmax = DDIV(temp, denom,  dfail)
        temp  = DDIV(a(imin,imin), epsinv, dfail)
        spmin = DDIV(temp,         denom,  dfail)
        spinfa = -one
        spinfb = -one
        spdsda = -one
        spdsdb = -one
        initl  = .true.
        return
      end if
c
c ##RECORD THE MAXIMUM SMOOTHING PARAMETER THAT CAUSES A WEAK   MINIMUM
c ##RECORD THE MINIMUM SMOOTHING PARAMETER THAT CAUSES A STRONG MINIMUM
c
      if (info .eq. 0) then
        if ((sp(1) .lt. spinfa) .or. (spinfa .lt. zero)) spinfa = sp(1)
      else if (sp(1) .gt. spinfb) then
        spinfb = sp(1)
      end if
c
c ##RECORD  THE MAXIMUM  SMOOTHING PARAMETER THAT CAUSES A NEGATIVE 
c ##DISCREPANCY DIFFERENCE.  RECORD THE MINIMUM SMOOTHING PARAMETER 
c ##THAT CAUSES A POSITIVE DISCREPANCY DIFFERENCE.
c
      if ((dsd(1) .lt. zero) .and. (sp(1) .gt. spdsdb)) then
        spdsdb = sp(1)
        dsdb   = dsd(1)
      else if (dsd(1) .gt. zero) then
        if ((sp(1) .lt. spdsda) .or. (spdsda .lt. zero)) then
          spdsda = sp(1)
          dsda   = dsd(1)
        end if
      end if
c
c
      if ((.not. bnd) .or. (info .eq. 1)) then
        sp(2)  = sp(1)
        sp(3)  = sp(1)
        dsd(2) = dsd(1)
        dsd(3) = dsd(1)
      end if
c
c ##IF WE HAVE A WEAK SOLUTION, THEN WE'LL NEED TO INCREASE
c ##THE SMOOTHING PARAMETER
c
      if (info .eq. 1) then
        initl = .true.
        if (spinfa .lt. zero) then
          sp(1) = sp(1) * ten**sclrg
        else
          sp(1) = dsqrt(spinfa * spinfb)
        end if
c
c ##IF WE HAVEN'T BOUNDED THE OPTIMAL SMOOTHING PARAMETER 
c ##WE USE THE MONOTONICITY OF THE DISCREPANCY DIFFERENCE
c ##AS A  FUNCTION OF THE SMOOTHING PARAMETER TO INCREASE 
c ##OR DECREASE PREVIOUS SMOOTHING PARAMETER CHOICE.
c
      else if (.not. bnd) then
        if ((dsd(1) .ge. zero) .and. (spinfb .ge. zero)) then
          sp(1) = dsqrt(sp(1) * spinfb)
        else
          if (iter .eq. 1) scale = zero
          scale = dabs(scale) + one
          if (scale .gt. sclrg) scale = sclrg
          if (dsd(1) .gt. zero) scale = -scale
          sp(1)  = sp(1) * ten**(scale)
        end if
c
c
      else if (spdsda .le. spinfb) then
        initl  = .true.
        sp(2)  = sp(1)
        sp(3)  = sp(1)
        dsd(2) = dsd(1)
        dsd(3) = dsd(1)
        sp(1) = dsqrt(spinfa * spinfb)
c
c ##else we have the 'proper' smooothing parameter bounded and we can
c ##use the root finding algorithm DBRENT
c
      else
        if (initl) then
          initl = .false.
          dsd(1) = dsda
          dsd(2) = dsdb
          dsd(3) = dsd(1)
          sp(1)  = spdsda
          sp(2)  = spdsdb 
        end if
        dsp1   = dlog(sp(1))
        dsp2   = dlog(sp(2))
        dsp3   = dlog(sp(3))
        call DBADJ(dsd(1), dsd(2), dsd(3), dsp1, dsp2, dsp3, e, d)
        tol1   = 2.d00*rmin*dabs(dsp1) + tol(3)*0.5d00
        call DBRCOR(dsd(1), dsd(2), dsd(3), dsp1, dsp2, dsp3,
     &              e, d, tol1)
        sp(1)  = dexp(dsp1)
        sp(2)  = dexp(dsp2)
        sp(3)  = dexp(dsp3)
      end if
c
      if (sp(1) .lt. spmin) then
        call WARNG(18, nlog(2), msglev(2), 0, sp(1), 'small')
      else if (sp(1) .gt. spmax) then
        call WARNG(18, nlog(2), msglev(2), 0, sp(1), 'large')
      end if
      if ((spinfa .gt. zero) .and. (spinfa .lt. spinfb)) 
     &  call ERROR(112, nlog(2), msglev(2), 0, 0.d+0, ' ', fail(1))
      return
*
*___end of FSPDSC___*
*
      end
