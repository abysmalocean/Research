      function          RLAM(x)
      double precision  RLAM,  x
      dimension         x(*) 
************************************************************************
*
* RLAM computes R given in Eq. 1.8 in Craven and Wahba, Numer.Math., 31, 
* 377-403.  Special attention is given to avoid underflow.
*
* x:  1             log(p*lambda)
*     2             p
*     3             sigsq
*     4    -  p+3   sigma(ii)**2
*     p+4  - 2p+3   y(i)**2
*     2p+4 - 3p+3   workspace
*    
* LAST MODIFIED: 11 NOV 1988
*
************************************************************************
      integer           i,      iarr1,  iarr2,  iarr3,  p
      double precision  plam,   sigsq,  temps,  tempy,  zero
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      plam   = dexp(x(1))
      p      = nint(x(2))
      sigsq  = x(3)
      iarr1  = 3
      iarr2  = p + iarr1
      iarr3  = p + iarr2
      tempy  = zero
      temps  = zero
c
c ##FIRST FORM USEFUL TEMPORARY VARIABLES
c
      do 100 i = 1, p
        x(iarr3+i) = plam / (plam + x(iarr1+i))
  100 continue
c
c ##FORM THE TWO SUMS
c
      do 200 i = 1, p
        tempy = tempy + x(iarr3+i) * x(iarr3+i) * x(iarr2+i)
        temps = temps + x(iarr3+i)
  200 continue
      RLAM = sigsq + (tempy - 2.d+0 * sigsq * temps) / x(2)
      return
c
*___end of RLAM___*
c
      end
      function          VLAM(x)
      double precision  VLAM,  x
      dimension         x(*) 
************************************************************************
*
* VLAM computes V(lambda), where V is given in equation 2.3, Golub 
* et. al., Technometrics, 21, 215-223.  Here we assume the dimension 
* of the solution vector is greater than p.
*
* x:  1             log(p*lambda)
*     2             p
*     3             not used
*     4    -  p+3   sigma(ii)**2
*     p+4  - 2p+3   y(i)**2
*     2p+4 - 3p+3   workspace
*
* LAST MODIFIED:12 JAN 1989
************************************************************************
      integer           i,      iarr1,  iarr2,  iarr3,  p
      double precision  denom,  dnumer, plam,   zero
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      plam   = dexp(x(1))
      p      = nint(x(2))
      iarr1  = 3
      iarr2  = p + iarr1
      iarr3  = p + iarr2
      denom  = zero
      dnumer = zero
c
c ##PREPARE USEFUL TEMPORARY VARIABLES
c
      do 100 i = 1, p
        x(iarr3+i) = plam / (x(iarr1+i) + plam)
100   continue
c
c ##COMPUTE NUMERATOR AND DENOMINATOR
c
      do 200 i = 1, p
        dnumer = dnumer + x(iarr3+i) * (x(iarr3+i) * x(iarr2+i))
        denom  = denom  + x(iarr3+i)
200   continue
      VLAM  = x(2) * dnumer / denom / denom
      return
c
*___end of VLAM___*
c
      end
      subroutine BRAKMN(F,      delta,  x,      xlow,   xhigh,  ax,     
     &                  bx,     cx,     inform)
      integer           inform
      double precision  ax,     bx,     cx,     delta,  F,      x,      
     &                  xhigh,  xlow
      dimension         x(*)
************************************************************************
*
* BRAKMN performs a linear global search to bracket the minimum of the
* function F(x).
*
* LAST MODIFIED: 5 DEC 1988
************************************************************************
      double precision  temp,   xmin,   ymin
*-----------------------------------------------------------------------
      inform = 0
      x(1)   = xlow + delta
      xmin   = x(1)
      ymin   = F(x)
c
c ##DO UNTIL x(1) .GE. xhigh
c
100   x(1) = x(1) + delta
      if (x(1) .lt. xhigh) then
        temp = F(x)
        if (temp .lt. ymin) then
          ymin = temp
          xmin = x(1)
        end if
        go to 100
      end if
c
c ##SEE IF THE ENDPOINTS ARE THE MINIMUM
c
      x(1) = xlow
      temp = F(x)
      if (temp .lt. ymin) then
        ymin   = temp
        xmin   = xlow
        inform = 1
      end if
      x(1) = xhigh
      temp = F(x)
      if (temp .lt. ymin) then
        ymin   = temp
        xmin   = xhigh
        inform = 2
      end if
      bx = xmin
      ax = dmax1(bx-delta,xlow)
      cx = dmin1(bx+delta,xhigh)
      return
c
*___end of BRAKMN___*
c
      end
      subroutine PRERV (fail,   mdat,   mdim,   msglev, ndim,   nlog,   
     &                  norm,   p,      sigsq,  wmach,  work,   yhatm)
      integer           fail,   mdat,   mdim,   msglev, ndim,   nlog,   
     &                  norm,   nsort,  p
      double precision  sigsq,  wmach,  work,   yhatm
      dimension         fail(*), msglev(*), nlog(*), wmach(*), 
     &                  work(mdim,*), yhatm(*)
************************************************************************
* PRERV is called to load the vector required to compute RLAM (or VLAM).
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           iadd,   IDAMAX, inform, is,     isort,  isortm, 
     &                  iu,     iw,     iy,     pp4
*-----------------------------------------------------------------------
c
c ##INITIALIZE THE WORKSPACE ARRAYS
c
      iu    = mdim + 1
      is    = mdim + mdim + 1
      iw    = mdim + mdim + 3
      iy    = mdim + mdim + 4
      pp4   = p + 4
c
c ##COMPUTE THE LEFT HAND EIGENVECTORS AND EIGENVALUES OF dfkhat
c
      call FDINV (mdim, ndim, p, norm, work)
      call EISSVD(mdim,   p,      ndim,   work, work(1,is), 
     &            .true., work(1,iu),     .false.,      work,   
     &            inform,   work(1,iw))
      if (inform .ne. 0) then
        call ERROR(64, nlog(2), msglev(2), inform, 0.d+0, ' ', fail(1))
        return
      end if
c
c ##NEXT PUT THE ELEMENTS p OF u WITH THE LARGEST EIGENVALUES AT THE
c ##BEGINNING
c
      do 100 isort = 1, p 
        nsort = ndim + 1 - isort
        iadd  = IDAMAX(nsort, work(isort,is), 1) - 1
        if (iadd .ne. 0) then
          isortm = isort - 1
          call DSWAP(1, work(isort,is), 1, work(isort+iadd,is), 1)
          call DSWAP(p, work(1,iu+isortm),1, work(1,iu+isortm+iadd),1)
        end if
100   continue
      call TMMMLT(mdim,  mdat,  mdim,  p,  1,  p, 
     &            work(1,iu),  yhatm,  work(1,iy))
c
c ##LOAD THE ARRAYS NEEDED TO COMPUTE RLAM AND VLAM
c
      work(2,1) = dble(p)
      work(3,1) = sigsq
      call DCOPY(p, work(1,is), 1, work(4,1),  1)
      call DDSCL(p, work(1,is), 1, work(4,1),  1)
      call DCOPY(p, work(1,iy), 1, work(pp4,1), 1)
      call DDSCL(p, work(1,iy), 1, work(pp4,1), 1)
      return
*
*___End of PRERV___*
*
      end
