      subroutine FRSTSP(fail,   mdat,   mdim,   msglev, ndim,   nlog,   
     &                  opchr,  p,      drvhst, dspc,   fkhat,  smp,    
     &                  tnrmsq, tol,    wmach,  work,   yhatm)
      integer           fail,   mdat,   mdim,   msglev, ndim,   nlog,   
     &                  opchr,  p
      double precision  drvhst, dspc,   fkhat,  smp,    tnrmsq, tol,    
     &                  wmach,  work,   yhatm
      dimension         fail(*), msglev(*), nlog(*), opchr(*), 
     &                  drvhst(0:*), dspc(*), fkhat(mdat,*), tol(*), 
     &                  wmach(*), work(mdim,*), yhatm(*)
************************************************************************
*
* FRSTSP chooses the smoothing parameter via 
* the minimization of Eq. (1.8), Craven & Wahba, Numer. Math. 31,
* 377-403.
*
* LAST MODIFIED: 12 FEB 1989
************************************************************************
      integer           i,      inform, ismth,  norm
      double precision  dlgcst, RLAM,   rvmin,  sigsq,  smptmp, VLAM,   
     &                  xmin,   zero
      logical           sigkn
      parameter        (zero   = 0.d+0)
*-----------------------------------------------------------------------
c ##COMPUTE THE INVERSE OF THE SMOOTHING MATRIX
c
      if (fail(1) .gt. 0) return
      call SCRNM(26, nlog(5), msglev(5), 0, 0.d+0, ' ')
      ismth = iabs(opchr(1))
      if ((ismth .eq. 4) .and. (smp .lt. -wmach(5))) then
        smp = -smp
        return
      end if
      norm = opchr(12)
c
c ##PREPARE TO LOAD THE ARRAYS (PRERV) NEEDED TO CALCULATE 
c ##RLAM OR VLAM
c
      do 100 i = 1, p
        call DCOPY(ndim, fkhat(i,1), mdat, work(i,1), mdim)
100   continue
      sigkn = (ismth .eq. 3) .or. (ismth .eq. 4) .or.
     &        (ismth .eq. 7)
      sigsq = tnrmsq / dble(p)
      call PRERV (fail,   mdat,   mdim,   msglev, ndim,  nlog,   
     &            norm,   p,      sigsq,  wmach,  work,  yhatm)
      if (fail(1) .gt. 0) return
c
c ##BOUND THE SMOOTHING PARAMETER 
c
      dlgcst = dlog( dble(p) / (dspc(7)**(norm + norm - 1)))
      call MINRV( inform, msglev, nlog,   p,      dlgcst, dspc,   
     &            rvmin,  tol,    wmach,  work,   xmin,   sigkn)
c
c ##IF THE INPUT SMOOTHING PARAMETER IS NEGATIVE, THEN MICRON
c ##IGNORES THE CALCULATED MINIMUM
c
      smptmp = dexp(xmin - dlgcst) 
      if (smp .ge. zero) then
        smp = smptmp
      else
        smp = dabs(smp)
        if ((ismth .eq. 5) .or. (ismth .eq. 7)) then 
          call WARNG(35, nlog(2), msglev(2), 0, smptmp, ' ')
          work(1,1) = dlog(smp) + dlgcst
          if (ismth .eq. 5) then
            rvmin = VLAM(work)
          else
            rvmin = RLAM(work)
          end if
        end if
      end if
      if ((ismth .eq. 5) .or. (ismth .eq. 7)) drvhst(1) = rvmin
      return
c
*___end of FRSTSP___*
c
      end
      subroutine FDINV (m,      n,      p,      k,     x)
      integer           k,      m,      n,      p
      double precision  x
      dimension         x(m,*)
************************************************************************
* FDINV transforms the matrix x by the inverse of the k'th order finite
* difference matrix
* INPUTS:
*      x: the matrix that will be transformed
*      n: the number of columns of x
*      m: the declared dimension of x in the calling program
*      k: the order of the finite difference matrix
*      p: the number of rows of x
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           i,      inorm,  j
      double precision  one,    zero
      parameter        (one  = 1.d+0)
      parameter        (zero = 0.d+0)
c-----------------------------------------------------------------------
      do 200 inorm = 1, k
c
c ##TRIANGULAR MATRIX IS UPPER TRIANGULAR
c
        if (mod(inorm,2) .eq. 0) then
          do 210 i = 1, n
            x(i,1) = -x(i,1)
            do 211 j = 2, n
              x(i,j) = x(i,j-1) - x(i,j) 
211         continue
210       continue
        else
c
c ##TRIANGULAR MATRIX IS LOWER TRIANGULAR
c
          do 220 i = 1, n
            do 221 j = n-1, 1, -1
              x(i,j) = x(i,j) + x(i,j+1)
221         continue
220       continue
        end if
200   continue
      return
c
*___end of FDINV___*
c
      end
      subroutine MINRV( inform, msglev, nlog,   p,      dlgcst, dspc,   
     &                  rvmin,  tol,    wmach,  work,   xmin,   sigkn)
      integer           inform, msglev, nlog,   p
      double precision  dlgcst, dspc,   rvmin,  tol,    wmach,  work,   
     &                  xmin
      logical           sigkn
      dimension         msglev(*), nlog(*), dspc(*), tol(*), wmach(*), 
     &                  work(*)
************************************************************************
* MINRV performs a global minimization to find the correct smoothing 
* parameter
*
* LAST MODIFIED: 9 FEB 1989
************************************************************************
      integer           infbr
      double precision  ax,     BRENT,  bx,     cx,     delta,  dlgmax, 
     &                  dlgmin, RLAM,   VLAM
      external          RLAM,   VLAM
*-----------------------------------------------------------------------
      inform = 0
      dlgmin = dlog(wmach(6)) + dlgcst
      dlgmax = dlog(wmach(8)) + dlgcst
      if (work(p+3) .gt. wmach(5))
     &  dlgmin = dmax1( dlgmin, dlog(work(p+3)) )
      dlgmax = dmin1( dlgmax, dlog(work(4)) )
      delta  = (dlgmax - dlgmin) / 50.d+0
      if (delta .gt. 1.d+0) delta = 1.d+0
c
c ##MINIMIZE RLAM OR VLAM.  FIRST PERFORM GLOBAL MINIMIZATION AND
c ##IDENTIFY THE INTERVAL CONTAINING THE MINIMUM
c
      if (sigkn) then
        call BRAKMN(RLAM,  delta,  work,   dlgmin, dlgmax, 
     &              ax,    bx,     cx,     inform)
        if (inform .eq. 0) then
          rvmin = BRENT(ax,     bx,    cx,     RLAM,   tol(2), xmin,   
     &                  work,   infbr)
          if (infbr .lt. 0) 
     &      call WARNG(31, nlog(2), msglev(2), infbr, 0.d+0, 'RLAM ')
        else
          call WARNG(32, nlog(2), msglev(2), inform, 0.d+0, 'RLAM ')
          xmin      = bx
          work(1)   = bx
          rvmin     = RLAM(work)
        end if
      else
c
c ##MINIMIZE V(lambda)
c
        call BRAKMN(VLAM,  delta,  work,   dlgmin, dlgmax, 
     &              ax,    bx,     cx,     inform)
        if (inform .eq. 0) then
          rvmin = BRENT(ax,     bx,    cx,     VLAM,   tol(2), xmin,   
     &                  work,   infbr)
          if (infbr .lt. 0) 
     &      call WARNG(31, nlog(2), msglev(2), infbr, 0.d+0, 'VLAM ')
        else
          call WARNG(32, nlog(2), msglev(2), inform, 0.d+0, 'VLAM ')
          xmin      = bx
          work(1)   = bx
          rvmin     = VLAM(work) 
        end if
      end if
*
*___End of MINRV___*
*
      end
