      subroutine CPHSWS(fail,   iter,   kx,     kxr,    mdim,   mdat,   
     &                  msglev, ndim,   nlog,   nranka, nrankr, opchr,  
     &                  p,      a,      ainp,   dspc,   r,      smp,    
     &                  wmach,  work)
      integer           fail,   iter,   kx,     kxr,    mdat,   mdim,   
     &                  msglev, ndim,   nlog,   nranka, nrankr, opchr,  
     &                  p
      double precision  a,      ainp,   dspc,   r,      smp,
     &                  wmach,  work
      dimension         fail(*), kx(*), kxr(*), msglev(*), nlog(*), 
     &                  opchr(*), a(mdim,*), ainp(mdim,*), dspc(*), 
     &                  r(mdat,*), wmach(*), work(mdim,*)
************************************************************************
*
* CPHSWS the upper triangular matrix 'r' that forms the Hessian H, where
* H = ainp' ainp = r' r + lambda*D
* where D is the matrix that represents the penalty term.
*
* Inputs    :kx, a, r, sp(1)
* Outputs   :ainp
*
* LAST MODIFIED:26 Jan 1989
************************************************************************
      integer           i,      idiag,  inform, iplus,  
     &                  j,      mdimnd, mdimp1, ndiag
      double precision  one,    spmin,  tolrnk, zero
      parameter        (one    = 1.d+0)
      parameter        (zero   = 0.d+0)
      save              mdimnd, mdimp1, ndiag,  spmin
c-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (iter .eq. 1) then
        spmin = wmach(5)
        mdimnd = mdim * ndim
        mdimp1 = mdim + 1
        if (opchr(12) .ge. 0) then
          ndiag = opchr(12) + 1
        else
          ndiag = 4
        end if
      end if
      call DLOAD(mdimnd, zero, work, 1)
c
c ##CALCULATE THE SMOOTHING TERMS THAT WILL BE ADDED 
c ##ADD THEM TO THE UPPER HALF OF THE UNSMOOTHED HESSIAN
c
      call PENLTY(mdim,   ndim,   opchr(12),   p,      
     &            dspc(7), smp, work)
c 
c ##DECOMPOSE ainp
c
      if (smp .gt. spmin) then
        call DCOPY(mdimnd, a, 1, ainp, 1)
c
c ##THE PENALTY MATRIX HAS ndiag NONZERO DIAGONALS
c
        iplus = ndim
        do 300 idiag = 1, ndiag
          call DXPY(iplus, work(1,idiag), mdimp1, 
     &              ainp(1,idiag), mdimp1)
        iplus = iplus - 1
300     continue
        call IDLOAD(ndim, 1, 1, kx, 1)
        tolrnk = wmach(4)
        call LSCHOL(mdim, ndim, nranka, tolrnk, kx, ainp, 
     &              inform, nlog(5))
        if (inform .ne. 0) then
          call ERROR(96, nlog(2), msglev(2), 0, 0.d+0, 
     &               'LSCHOL', fail(1))
        else if (nranka .ne. ndim) then
          call WARNG(25,  nlog(2),   msglev(2), nranka,  0.d+0,   
     &               'LSCHOL')
c
c ##STORE THE SMOOTHING PARAMETER THAT CAUSED THE LSCHOL ERROR
c
          spmin = smp * (one + 10.d+0 * wmach(3))
        end if
      end if
c
c ##A LOSS OF PRECISION WAS PROBABLY ENCOUNTERED
c ##TRY GENERALIZED CHOLESKY
c
      if (smp .le. spmin) then
        call DCOPY(mdimnd, work, 1, ainp, 1)
c
c ##REORDER THE SMOOTHING TERMS TO CORRESPOND TO r
c ##PERFORM SYMMETRIC INTERCHANGES
c ##INTERCHANGES ARE EASIER IF THE UPPER-TRIANGULAR PORTION OF
c ##ainp IS COPIED TO THE LOWER TRIANGULAR PORTION
c
        call DCOPY(ndim, ainp, mdimp1, work(1,mdimp1), 1)
        do 500 i = 1, ndim
          call DCOPY(ndim-i, ainp(i,i+1), mdim, ainp(i+1,i), 1)
500     continue
        do 600 i = 1, ndim
          do 610 j = i+1, ndim
            if (kxr(i) .gt. kxr(j)) then
              ainp(i,j) = ainp(kxr(i),kxr(j))
            else
              ainp(i,j) = ainp(kxr(j),kxr(i))
            end if
610       continue
          ainp(i,i) = work(kxr(i),mdimp1)
600     continue
        call GCHOL (ndim,   nranka, nrankr, r,      mdat,   ainp,   
     &              mdim,  wmach,  inform)
        if (inform .ne. 0) then
          call ERROR(96, nlog(2), msglev(2), 0, 0.d+0, 'GCHOL', fail(1))
        else if (nranka .ne. ndim) then
          call WARNG(25,  nlog(2),   msglev(2), nranka,  0.d+0, 'GCHOL')
        end if
c
c ##ADD THE R TO THE PETURBATION MATRIX
c
        iplus = ndim
        do 700 i = 1, nrankr
          call DXPY(iplus, r(i,i), mdat, ainp(i,i), mdim)
          iplus = iplus - 1
700     continue
        call ICOPY(ndim, kxr, 1, kx, 1)
      end if
      return
c
*___end of CPHSWS___*
c
1     format(/,5(2x,d10.3))
2     format(/,5(2x,i4))
      end
      subroutine PENLTY(mdim,   ndim,   norm,   p,      delta,  smp,    
     &                  work)
      integer           mdim,   ndim,   norm,   p
      double precision  delta,  smp,    work
      dimension         work(mdim,*)
************************************************************************
c
c PENLTY loads the diagonal terms of the penalty matrix
c
************************************************************************
      integer           mdimp1
      double precision  const,  cp2,    cm4,    cp6,    cm15,   cp20,   
     &                  cp10,   cm45,   cp120,  cm210,  cp252
c-----------------------------------------------------------------------
      mdimp1 = mdim + 1
      if (norm .ge. 0) then
        const = dble(p) * smp / delta**(norm + norm - 1)
**********************************************************************
c
c THIS SECTION IS DESIGNED TO REPLACE ALL THE DIFFERENT NORM
c SECTIONS.  WE STILL NEED TO CALCULATE THE CONSTANTS THAT
c MAKE UP scldia and sclnrm AND WE NEED TO DECIDE HOW TO STORE
c AND REFERENCE THEM.  WE'LL MAKE THE CALCULATIONS GENERAL FOR
c UP TO MNORM AND WILL CHECK TO MAKE SURE P AND NDIM ARE LARGE
c ENOUGH.
**********************************************************************
c       do 100 ibeg = 1, norm
c         do 110 jbeg = ibeg, norm
c         ainp(ibeg,jbeg) = work(ibeg,jbeg) + const * sclnrm(k)
c10       continue
c00     continue
c       do 200 idiag = norm + 1, 1, -1
c         call DAXPY(ndiag, scldia(idiag), const, 0, 
c    &               ainp(idiag,norm+1), mdimp1)
c00     continue
c       do 300 iend = ndim, ndim + 1 - norm, -1
c         do 310 jend = iend, ndim + 1 - norm, -1
c           ainp(jend,iend) = work(jend,iend) + const * sclnrm(k)
c10       continue
c00     continue
**********************************************************************
      end if
      if (norm .eq. 0) then
        call DLOAD(ndim, const, work, mdimp1)
      else if (norm .eq. 1) then
        work(1,1)        = const
        cp2  =  2.d+0 * const
        call DLOAD(ndim-norm, -const, work(1,2), mdimp1)
        call DLOAD(ndim-norm-1, cp2, work(2,2), mdimp1)
        work(ndim,ndim) = const
      else if (norm .eq. 2) then
        work(1, 1) = + const 
        work(1, 2) = - const * 2.d+0
        work(2, 2) = + const * 5.d+0
        cm4  = -4.d+0 * const
        cp6  =  6.d+0 * const
        call DLOAD(ndim-norm,   const, work(1,3), mdimp1)
        call DLOAD(ndim-norm-1, cm4,   work(2,3), mdimp1)
        call DLOAD(ndim-norm-2, cp6,   work(3,3), mdimp1)
        work(ndim,   ndim)   = + const 
        work(ndim-1, ndim)   = - const * 2.d+0
        work(ndim-1, ndim-1) = + const * 5.d+0
      else if (norm .eq. 3) then
        work(1, 1) = + const 
        work(1, 2) = - const *  3.d+0
        work(1, 3) = + const *  3.d+0
        work(2, 2) = + const * 10.d+0
        work(2, 3) = - const * 12.d+0
        work(3, 3) = + const * 19.d+0
        cp6  =   6.d+0 * const
        cm15 = -15.d+0 * const
        cp20 =  20.d+0 * const
        call DLOAD(ndim-norm,   -const, work(1,4), mdimp1)
        call DLOAD(ndim-norm-1, cp6,    work(2,4), mdimp1)
        call DLOAD(ndim-norm-2, cm15,   work(3,4), mdimp1)
        call DLOAD(ndim-norm-3, cp20,   work(4,4), mdimp1)
        work(ndim, ndim)     = + const 
        work(ndim-1, ndim)   = - const *  3.d+0
        work(ndim-2, ndim)   = + const *  3.d+0
        work(ndim-1, ndim-1) = + const * 10.d+0
        work(ndim-2, ndim-1) = - const * 12.d+0
        work(ndim-2, ndim-2) = + const * 19.d+0
      else if (norm .eq. 5) then
        work(1, 1) = + const 
        work(1, 2) = - const *   5.d+0
        work(1, 3) = + const *  10.d+0
        work(1, 4) = - const *  10.d+0
        work(1, 5) = + const *   5.d+0
        work(2, 2) = + const *  26.d+0
        work(2, 3) = - const *  55.d+0
        work(2, 4) = + const *  60.d+0
        work(2, 5) = - const *  35.d+0
        work(3, 3) = + const * 126.d+0
        work(3, 4) = - const * 155.d+0
        work(3, 5) = + const * 110.d+0
        work(4, 4) = + const * 226.d+0
        work(4, 5) = - const * 205.d+0
        work(5, 5) = + const * 251.d+0
        cp10  =   10.d+0 * const
        cm45  =  -45.d+0 * const
        cp120 =  120.d+0 * const
        cm210 = -210.d+0 * const
        cp252 =  252.d+0 * const
        call DLOAD(ndim-norm,   -const, work(1,6), mdimp1)
        call DLOAD(ndim-norm-1, cp10,  work(2,6), mdimp1)
        call DLOAD(ndim-norm-2, cm45,  work(3,6), mdimp1)
        call DLOAD(ndim-norm-3, cp120, work(4,6), mdimp1)
        call DLOAD(ndim-norm-4, cm210, work(5,6), mdimp1)
        call DLOAD(ndim-norm-5, cp252, work(6,6), mdimp1)
        work(ndim,  ndim)   = + const 
        work(ndim-1,ndim)   = - const *   5.d+0
        work(ndim-2,ndim)   = + const *  10.d+0
        work(ndim-3,ndim)   = - const *  10.d+0
        work(ndim-4,ndim)   = + const *   5.d+0
        work(ndim-1,ndim-1) = + const *  26.d+0
        work(ndim-2,ndim-1) = - const *  55.d+0
        work(ndim-3,ndim-1) = + const *  60.d+0
        work(ndim-4,ndim-1) = - const *  35.d+0
        work(ndim-2,ndim-2) = + const * 126.d+0
        work(ndim-3,ndim-2) = - const * 155.d+0
        work(ndim-4,ndim-2) = + const * 110.d+0
        work(ndim-3,ndim-3) = + const * 226.d+0
        work(ndim-4,ndim-3) = - const * 205.d+0
        work(ndim-4,ndim-4) = + const * 251.d+0
      end if
      return
c
*___end of PENLTY___*
c
      end
      subroutine GCHOL (n,      nranka, nrankr, r,      mrowr,  s,      
     &                  mrows,  wmach,  inform)
      integer           inform, mrowr,  mrows,  n,      nranka, nrankr
      double precision  r,      s,      wmach
      dimension         r(mrowr,*), s(mrows,*), wmach(*)
************************************************************************
* Given D symmetric and R upper triangular GCHOL finds S upper triangular 
* that satisfies
* (R + S)'(R + S) = R'R + D
* without forming R'R.
*
* INPUTS:
*          r      - the upper triangualr matrix R
*          s      - the upper triangular matrix D
* OUTPUTS:
*          s      - the upper triangular matrix S
*          inform = 0, the factorization was completed successfully
*        
* LAST MODIFIED: 6 Jan 1989
************************************************************************
      integer           k,      kp1,    nmk
      double precision  d,      one,    scale,  zero
      parameter        (one  = 1.d+0)
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      inform = 0
      do 100 k = 1, nrankr
         if (dabs(s(k,k)) .le. wmach(5)) then
           d = zero
         else
           d = s(k,k) / (dabs(r(k,k)) + dsqrt(s(k,k) + r(k,k)*r(k,k)))
           if (r(k,k) .lt. zero) d = -d 
         end if
         s(k,k) = d
         if (k .lt. n) then
           nmk = n - k
           kp1 = k + 1
           call DAXPY(nmk, -d, r(k,kp1), mrowr, s(k,kp1), mrows)
           scale = one / (r(k,k) + d)
           call DSCAL(nmk, scale, s(k,kp1), mrows)
           call DSYRM1('U', nmk, s(k,kp1), mrows,
     &                        s(kp1,kp1), mrows)
           call DGERM1(nmk, nmk, s(k,kp1), mrows, r(k,kp1), mrowr, 
     &               s(kp1,kp1), mrows)
           call DGERM1(nmk, nmk, r(k,kp1), mrowr, s(k,kp1), mrows, 
     &               s(kp1,kp1), mrows)
         end if
100   continue
c
c ##USE SIMPLIFIED CALCULATIONS WHEN THE ROWS OF r ARE ZERO
c
      do 200 k = nrankr + 1, n
        if (s(k,k) .le. wmach(5)) then
          nranka = k - 1
          return
        else
          d = dsqrt(s(k,k))
        end if
        s(k,k) = d
        if (k .lt. n) then
           kp1 = k + 1
           scale = one / d
           call DSCAL(n-k, scale, s(k,kp1), mrows)
           call DSYRM1('U', n-k, s(k,kp1), mrows,
     &                        s(kp1,kp1), mrows)
        end if
200   continue
      nranka = n
      return
c
*___end of GCHOL___*
c
      end
