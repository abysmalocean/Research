      subroutine INVEK (fail,   ichpnt, ierrmd, mdat,   msglev, ndim,   
     &                  nlog,   npdat,  p,      einv,   fk,     fkhat,  
     &                  wmach,  work,   y,      yhatm)
      integer           fail,   ichpnt, ierrmd, mdat,   msglev, ndim,   
     &                  nlog,   npdat,  p
      double precision  einv,   fk,     fkhat,  y,      yhatm,  wmach,  
     &                  work
      dimension         fail(*), ichpnt(*), msglev(*), nlog(*), 
     &                  einv(mdat,*), fk(*), fkhat(*), wmach(*), 
     &                  work(*), y(*), yhatm(*)
************************************************************************
* INVEK scales the respose matrix, fk, by the inverse error matrix.
*
* LAST MODIFIED:19 OCT 1989
************************************************************************
      integer           idat,   ip,     isave
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      call SCRNM(29, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
c ##FIRST COPY THE USED ROWS OF fk TO work
c
      isave = 1
      idat  = 1
10    if ((isave .le. npdat) .or. (idat .le. mdat)) then
        if (ichpnt(idat) .gt. 0) then
          call DCOPY(ndim, fk(idat), mdat, work(isave), mdat)
          isave = isave + 1
        end if
        idat = idat + 1
        go to 10
      end if
      if (ierrmd .eq. 1) then
c
c ##  THE COMPUTATIONS ARE SIMPLE WHEN THE ERROR ARRAY IS DIAGONAL
c ##  NOTE WE NEED TO ACCOUNT FOR 0 STANDARD DEVIATIONS
c
        ip = 0
        do 200 idat = 1, npdat
          if (einv(idat,idat) .gt. wmach(5)) then
            ip = ip + 1
            call DCOPY(ndim, work(idat), mdat, fkhat(ip), mdat)
            call DSCAL(ndim, einv(idat,idat),  fkhat(ip), mdat)
            yhatm(ip) = y(idat)*einv(idat,idat)
          end if
200     continue
        if (ip .ne. p) write(6,*) ' SEVERE ERROR: IPNEP'
      else
c
c ##  ELSE WE NEED TO PERFORM SOME MATRIX MULTIPLICATIONS
c
        call MMMLT(mdat, mdat, mdat, npdat, ndim, p, einv, work, 
     &            fkhat)
        call MMMLT(mdat, mdat, mdat, npdat, 1,    p, einv, y, 
     &            yhatm)
      end if
*
*___end of INVEK___*
*
      end
      subroutine VLF   (mdat,   mdim,   ndim,   p,      f,      a,      
     &                  ainp,   cvec,   r,      vlfam,  vlfbm,  work,
     &                  yinf)
      integer           mdat,   mdim,   ndim,   p
      double precision  f,      a,      ainp,   cvec,   r,      vlfam,  
     &                  vlfbm,  work,   yinf
      dimension         f(*), a(mdim,*), ainp(mdim,*), cvec(*), 
     &                  r(mdat,*), work(*)
************************************************************************
* VLF evaluates the variance of the linear functional stored in f.
*
* LAST MODIFIED:17 OCT 1989
************************************************************************
      integer           i,      n
      double precision  DDOT,   temp
*-----------------------------------------------------------------------
      n = ndim
      do 100 i = 1, ndim
        temp = DDOT(n, f(i), 1, a(i,i), mdim)
        f(i) = temp
        n    = n - 1
100   continue
      vlfbm = dsqrt( DDOT(ndim, f, 1, f, 1) )
      call MMMLT (mdat, ndim, p, ndim, 1, p, r, f, work)
      call TMMMLT(mdim, p, p, p,  1,  p, ainp, work, f)
      call DDDIV ( p, cvec, 1, f, 1 )
      vlfam = dsqrt( vlfbm*vlfbm - DDOT(p, f, 1, f, 1) )
      yinf  = dlog( vlfbm / vlfam )
      return
*
*___end of VLF___*
*
      end
      subroutine UNORD (n,      kx,     x,      incx,   work)
      integer           n,      kx,     incx
      double precision  x,      work
      dimension         kx(*), x(*), work(*)
************************************************************************
* Given the order vector kx, UNORD computes
* x(kx(i)) = x(i)
*
* LAST MODIFIED: 1 NOV 1989
************************************************************************
      integer           i,      j
*-----------------------------------------------------------------------
      j = 1
      do 100 i = 1, n
        work(kx(i)) = x(j)
        j           = j + incx
100   continue
      call DCOPY(n, work, 1, x, incx)
      return
*
*___end of ORDER___*
*
      end
      subroutine ORDER (n,      kx,     x,      incx,   work)
      integer           n,      kx,     incx
      double precision  x,      work
      dimension         kx(*), x(*), work(*)
************************************************************************
* Given the order vector kx, ORDER computes
* x(i) = x(k(i))
************************************************************************
      integer           i,      j
*-----------------------------------------------------------------------
      do 100 i = 1, n
        j       = (kx(i) - 1) * incx + 1
        work(i) = x(j)
100   continue
      call DCOPY(n, work, 1, x, incx)
      return
*
*___end of ORDER___*
*
      end
      subroutine VMAT  (fail,   kx,     mdat,   mdim,   msglev, ndim,   
     &                  nlog,   p,      fkhat,  a,      ainp,   r,      
     &                  cvec,   wmach,  work)
      integer           fail,   kx,     mdat,   mdim,   msglev, ndim,   
     &                  nlog,   p
      double precision  fkhat,  a,      ainp,   r,      cvec,   work,   
     &                  wmach
      dimension         fail(*), kx(*), msglev(*), nlog(*), 
     &                  fkhat(mdat,*), a(mdim,*), ainp(mdim,*), 
     &                  r(mdat,*), cvec(*), work(*), wmach(*)
************************************************************************
* VMAT computes the matrices used to simplify the variance calculations
*
* LAST MODIFIED:17 OCT 1989
************************************************************************
      integer           i,      inform, j,      n,      nranka
      double precision  DDOT,   tolrnk
*-----------------------------------------------------------------------
c ##FIRST COMPUTE THE SQUARE ROOT OF THE COVARIANCE MATRIX
c
      if (fail(1) .gt. 0) return
      call IDLOAD(ndim, 1, 1, kx, 1)
      tolrnk = wmach(4)
      call LSCHOL(mdim, ndim, nranka, tolrnk, kx, a, 
     &            inform, nlog(5))
      if (inform .ne. 0) then
        call ERROR(96, nlog(2), msglev(2), 0, 0.d+0, 
     &             'LSCHOL', fail(1))
        return
      else if (nranka .ne. ndim) then
        call WARNG(25,  nlog(2),   msglev(2), nranka,  0.d+0,   
     &             'LSCHOL')
      end if
c
c ##REORDER fkhat
c
      do 100 i = 1, p
        call ORDER (ndim, kx, fkhat(i,1), mdat, work)
100   continue
c
c ##MULTIPLY fkhat AND a TRANSPOSE AND PLACE IN r
c
      do 300 i = 1, p
        n = ndim
        do 310 j = 1, ndim 
          r(i,j)    = DDOT(n, fkhat(i,j), mdat, a(j,j), mdim)
          ainp(i,j) = r(i,j)
          n         = n - 1
310     continue
300   continue
c
c ##PERFORM SVD ON r
c
      call EISSVD(mdim,   p,      ndim,   ainp, cvec,
     &            .true., ainp,   .false., work,   
     &            inform,   work(1))
      if (inform .ne. 0) then
        call ERROR(64, nlog(2), msglev(2), inform, 0.d+0, ' ', fail(1))
        return
      end if
      do 400 i = 1, p
        cvec(i) = dsqrt( 1.d+0 + cvec(i) * cvec(i) )
400   continue
      return
*
*___end of VMAT___*
*
      end
      subroutine CVVEC (fail,   kx,     mdat,   mdim,   ndim,   nlf,    
     &                  p,      a,      ainp,   cvec,   diam,   f,      
     &                  r,      vlfam,  vlfbm,  work,   yinf)
      integer           fail,   kx,     mdat,   mdim,   ndim,   nlf,    
     &                  p
      double precision  a,      ainp,   cvec,   diam,   f,      r,      
     &                  vlfam,  vlfbm,  work,   yinf
      dimension         fail(*), kx(*), a(mdim,*), ainp(mdim,*), 
     &                  cvec(*), diam(*), f(*), r(mdat,*), vlfam(*), 
     &                  vlfbm(*), work(*), yinf(*)
************************************************************************
* CVVEC creates and reports the varince and information vectors.
*
* LAST MODIFIED:17 OCT 1989
************************************************************************
      integer           i,      ibeg,   iend,   nspc
      double precision  delta,  one,    zero
      parameter        (one  = 1.d+0)
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
c ##FIRST CALCULATE THE VARIANCE OF EACH LINEAR SPLINE
c
      if (fail(1) .gt. 0) return
      delta = diam(2) - diam(1)
      do 100 i = 1, ndim
        call DLOAD(ndim, zero, f, 1)
        f(i) = one
        call ORDER (ndim,   kx,     f,      1,      work)
        call VLF   (mdat,   mdim,   ndim,   p,      f,      a,      
     &              ainp,   cvec,   r,      vlfam(i),  vlfbm(i),  
     &              work,   yinf(i))
100   continue
      nlf = ndim
c
c ##NEXT COMPUTE THE VARIANCE IN THE FIRST 3 MOMENTS
c
      do 200 i = nlf+1, nlf+3
        call DLOAD(ndim, one, f, 1)
        f(1)    = 0.5d+0
        f(ndim) = 0.5d+0
        call DSCAL(ndim, delta, f, 1)
        call DPSCAL(ndim, i-nlf-1, diam, 1, f, 1 )
        call ORDER (ndim,   kx,     f,      1,      work)
        call VLF   (mdat,   mdim,   ndim,   p,      f,      a,      
     &              ainp,   cvec,   r,      vlfam(i),  vlfbm(i),  
     &              work,   yinf(i))
200   continue
      nlf = nlf + 3
c
c ##NEXT COMPUTE THE VARIANCE OF THE INTEGRAL OVER
c ##THIRDS OF THE SOLUTION INTERVAL
c
      do 300 i = nlf+1, nlf+3
        call DLOAD(ndim, zero, f, 1)
        ibeg    = ((i - nlf - 1) * ndim) / 3 + 1
        iend    = ((i - nlf) * ndim) / 3
        call DLOAD(iend-ibeg+1, one, f(ibeg), 1)
        f(1)    = 0.5d+0
        f(iend) = 0.5d+0
        call DSCAL(ndim, delta, f, 1)
        call ORDER (ndim,   kx,     f,      1,      work)
        call VLF   (mdat,   mdim,   ndim,   p,      f,      a,      
     &              ainp,   cvec,   r,      vlfam(i),  vlfbm(i),  
     &              work,   yinf(i))
300   continue
      nlf = nlf + 3
c
c ##NEXT COMPUTE THE VARIANCE OF USER DEFINE FUNCTIONALS
c
      nspc = 1
400   if (nspc  .gt. 0) call LSPC(nspc, ndim, f)
      if (nspc .gt. 0) then
        i = nlf + nspc
        call ORDER (ndim,   kx,     f,      1,      work)
        call VLF   (mdat,   mdim,   ndim,   p,      f,      a,      
     &              ainp,   cvec,   r,      vlfam(i),  vlfbm(i),  
     &              work,   yinf(i))
        go to 400
      end if
      nlf  = nlf + iabs(nspc) -1
      return
*
*___end of CVVEC___*
*
      end
      subroutine VREP  (fail,   msglev, ndim,   nlf,    nlog,   vlfam,  
     &                  vlfbm,  yinf)
      integer           fail,   msglev, ndim,   nlf,    nlog
      double precision  vlfam,  vlfbm,  yinf
      dimension         fail(*), msglev(*), nlog(*), vlfam(*), 
     &                  vlfbm(*), yinf(*)
************************************************************************
* VREP reports the results of the information calculations..
*
* LAST MODIFIED:10 OCT 1989
************************************************************************
      integer           i,      nout
*-----------------------------------------------------------------------
      nout = nlog(3)
      call SCRNM(36, nlog(5), msglev(5), 0, 0.d+0, ' ')
      if (fail(1) .le. 0) then
        write(nout,3)
        write(nout,2)
        do 100 i = 1, ndim
          write(nout,1) vlfbm(i), vlfam(i), yinf(i)
100     continue
        write(nout,4)
        write(nout,2)
        do 200 i = ndim+1, ndim+3
          write(nout,1) vlfbm(i), vlfam(i), yinf(i)
200     continue
        write(nout,5)
        write(nout,2)
        do 300 i = ndim+4, ndim+6
          write(nout,1) vlfbm(i), vlfam(i), yinf(i)
300     continue
        write(nout,6)
        if (nlf .le. ndim+6) then
          write(nout,7)
        else
          write(nout,2)
          do 400 i = ndim+7, nlf
            write(nout,1) vlfbm(i), vlfam(i), yinf(i)
400       continue
        end if
      end if
c
c ##CLOSE FILES
c
      if ((fail(2) .le. 0) .and. (fail(3) .le. 0)) then
        if (msglev(2) .gt. 0) close(nlog(2))
        if (msglev(3) .gt. 0) close(nout)
      end if
      return
1     format(1x,d10.3,t15,d10.3,t30,d10.3)
2     format(1x,'Before data',t15,'After data',t30,'Information')
3     format(//1x,'Spline values:')
4     format(//1x,'Moments:')
5     format(//,1x,'1/3 interval integrals:')
6     format(//,1x,'Special functionals:')
7     format(1x,'None')
*
*___end of VREP___*
*
      end
      subroutine SAMPLE(fail,   kx,     mdim,   ndim,   nseed,  nsamp,  
     &                  a,      diam,   f,      work)
      integer           fail,   kx,     mdim,   ndim,   nseed,  nsamp
      double precision  a,      diam,   f,      work
      dimension         fail(*), kx(*), a(mdim,*), diam(*), f(*), 
     &                  work(mdim,*)
************************************************************************
* SAMPLE generates random sample described by the autocorrelation
* function, and a mean response described in an .sim file.
*
* LAST MODIFIED:19 OCT 1989
************************************************************************
      integer           i,      j,      nsmp
      double precision  DDOT,   temp,   TSTFN
      character*3       chrnum
      character*10      fname
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (nsamp   .le. 0) return
      do 100 i = 1, ndim
        f(i) = TSTFN(diam(i))
        work(i,2) = dexp(diam(i))
100   continue
      call ORDER (ndim, kx, f, 1, work)
      call FUNIT(nsmp)
      call DLOAD(ndim, 0.d+0, work, 1)
      do 200 i = 1, nsamp + 1
        call INTCHR(i, chrnum)
        fname = 'sample.'
        call MERGE(fname, chrnum)
        if (i .ne. 1) call DNLOAD(ndim, nseed,  work, 1)
        do 210 j = 1, ndim
          temp    = f(j) + DDOT(ndim+1-j, a(j,j), mdim, work(j,1), 1)
          work(j,1) = temp
210     continue
        open(nsmp, file = fname, status='unknown')
        call UNORD (ndim, kx, work, 1, work(1,3))
        do 300 j = 1, ndim
          write(nsmp,1) work(j,2), work(j,1)
300     continue
        close(nsmp)
200   continue
      return
1     format(2(2x,d10.3))
*
*___end of SAMPLE___*
*
      end
