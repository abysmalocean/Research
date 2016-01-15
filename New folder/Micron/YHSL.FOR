      subroutine HESLIN(fail,   ichpnt, ierrmd, iw,     kx,     kxr,    
     &                  mdat,   mdim,   msglev, ndim,   nlog,   npdat,  
     &                  nrankr, opchr,  p,      a,      ainp,   b,      
     &                  cvec,   einv,   fk,     fkhat,  r,      ssq1,   
     &                  wmach,  work,   y,      yhatm)
      integer           fail,   ichpnt, ierrmd, iw,     kx,     kxr,    
     &                  mdat,   mdim,   msglev, ndim,   nlog,   npdat,  
     &                  nrankr, opchr,  p
      double precision  a,      ainp,   b,      cvec,   einv,   fk,     
     &                  fkhat,  r,      ssq1,   y,      yhatm,  wmach,  
     &                  work
      dimension         fail(*), ichpnt(*), iw(*), kx(*), kxr(*), 
     &                  msglev(*), nlog(*), opchr(*), a(mdim,*), 
     &                  ainp(mdim,*), b(*), cvec(*), einv(mdat,*), 
     &                  fk(*), fkhat(*), r(mdat,mdim), wmach(*), 
     &                  work(*), y(*), yhatm(*)
************************************************************************
*
* HESLIN determines the Hessian, "a", and linear, "cvec", term in the 
* quadratic functional ||[einv] [[fk][f] - [y]]||^2.
*
* INPUTS:
* ierrmd = 1 if the error array is diagonal
* p      - the number of independent error sources
* einv   - the inverse of the error array
* y      - the data
*
* OUTPUTS:
* fkhat  - the kernel function matrix transformed by the inverse error
*          array
* yhatm  - the data vector transformed by the inverse error array
*
* LAST MODIFIED:12 JAN 1989
************************************************************************
      integer           i,      idat,   IDRANK, info,   ip,     isave,  
     &                  j,      jmax,   jsave,  mdtdm,  ntemp,  nxchng
      double precision  DNRM2,  one,    tolrnk, zero
      parameter        (one  = 1.d+0)
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      if ((fail(1) .gt. 0) .or. (opchr(6) .le. 0)) return
      call SCRNM(25, nlog(5), msglev(5), 0, 0.d+0, ' ')
      mdtdm = mdat * ndim
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
c
c ##COMPUTE THE QUADRATIC AND LINEAR TERMS
c
      call TMMMLT(mdat, mdat, mdim, ndim, ndim, p, fkhat, fkhat,
     &           a)
      call TMMMLT(mdat, mdat, mdim, ndim, 1,    p, fkhat, yhatm,
     &           cvec)
      call DSCAL(ndim, -1.d00, cvec, 1)
c
c ##PREPARE FOR LEAST SQUARES
c
      call DCOPY(p, yhatm, 1, b, 1)
c
c ##PERFORM QR ON fkhat
c
      call DCOPY(mdtdm, fkhat, 1, r, 1)
      call IDLOAD(ndim, 1, 1, kxr, 1)
      call DGEQRP( 'Column iterchanges', p, ndim, r, mdat,
     &              work, iw, work(ndim+1), info )
      nxchng = min0(p, ndim)
      do 230 j = 1, nxchng
         jmax = iw(j)
         if (jmax .gt. j) then
            jsave     = kxr(jmax)
            kxr(jmax) = kxr(j)
            kxr(j)    = jsave
         end if
  230 continue
c
c ##TRANSFORM THE DATA VECTOR
c
      ntemp = min(p,ndim)
      call DGEAPQ( 'transpose', 'separate', p, ntemp, r, mdat,
     &             work, 1, b, p, work(ndim+1), info )
      tolrnk = wmach(3) * 10.d+0
      nrankr = IDRANK( ntemp, r, mdat+1, tolrnk )
      if (p .gt. nrankr) then
        ssq1 = DNRM2 ( p-nrankr, b(nrankr+1), 1 )
      else
        ssq1 = zero
      end if
      call ICOPY(ndim, kxr, 1, kx, 1)
      do 300 i = 1, ndim
        call DCOPY(nrankr, r(1,i), 1, ainp(1,i), 1)
300   continue
      return
c
*___end of HESLIN___*
c
      end
      subroutine MMMLT (mrowa,  mrowb,  mrowc,  ncola,  ncolb,  nrowa, 
     &                  a,      b,      c)
************************************************************************
c
c MMMLT computes "c" where "c" = "a" x "b",  where "a" is an "nrowa" 
c by "ncola" matrix and "b" is an "ncola" by "ncolb" matrix.
c
c INPUTS:
c mrowa  - the declared dimension of the columns of a in the calling
c          program
c mrowb  - analagous to mrowa
c mrowc  - analagous to mrowa
c
c LAST MODIFIED:11 NOV 1987
c
************************************************************************
      integer           irowa,  jcolb,  mrowa,  mrowb,  mrowc,
     &                  ncola,  ncolb,  nrowa
      double precision  a,      b,      c,      DDOT
      dimension         a(mrowa,ncola), b(mrowb,ncolb), c(mrowc,ncolb)
c
      do 100 irowa = 1, nrowa
        do 110 jcolb = 1, ncolb
          c(irowa, jcolb) = DDOT(ncola,a(irowa,1),mrowa,b(1,jcolb),1)
110     continue
100   continue
      return
c
*___end of MMMLT___*
c
      end
      subroutine TMMMLT(mrowa,  mrowb,  mrowc,  ncola,  ncolb,  nrowa, 
     &                  a,      b,      c)
************************************************************************
c
c TMMMLT computes "c" where "c" = transpose("a")  x "b",  where "a" is 
c an "nrowa" by "ncola" matrix and "b" is an "ncola" by "ncolb" matrix.
c
c INPUTS:
c mrowa  - the declared dimension of the columns of a in the calling
c          program
c mrowb  - analagous to mrowa
c mrowc  - analagous to mrowa
c
c LAST MODIFIED:11 NOV 1987
c
************************************************************************
      integer           icola,  jcolb,  mrowa,  mrowb,  mrowc,  ncola,  
     &                  ncolb,  nrowa
      double precision  a,      b,      c,      DDOT
      dimension         a(mrowa,ncola), b(mrowb,ncolb), c(mrowc,ncolb)
c
      do 100 icola = 1, ncola
        do 110 jcolb = 1, ncolb
          c(icola,jcolb) = DDOT(nrowa,a(1,icola),1,b(1,jcolb),1)
110     continue
100   continue
      return
c
*___end of TMMMLT___*
c
      end
      double precision function DLNGAM (x)
************************************************************************
*
* DLNGAM evaluates the log of the gamma function for positive arguments
* A modified version of the IMSL routine.
*
* LAST MODIFIED: 11 NOV 1987
************************************************************************
      integer           j
      double precision  p1, q1, p2, q2, p3, q3,  p4,  pi,  x,    rmax
      double precision  y,  t,  a,  b,  top, den, rmin
      double precision  one
      parameter        (one = 1.d+0)
      common /cmrmin/   rmin
      save   /cmrmin/
      common /cmrmax/   rmax
      save   /cmrmax/
      dimension         p1(9), q1(8), p2(9), q2(8), p3(9), q3(8), p4(7)
      save              p1, q1, p2, q2, p3, q3, p4, pi
c
c ##coefficients for minimax approximation to ln(gamma(x)),
c ##0.5 .le. x .le. 1.5
c
      data               p1(1)/6.304933722864032d02/,
     &                   p1(2)/1.389482659233250d02/,
     &                   p1(3)/-2.331861065739548d03/,
     &                   p1(4)/-2.651470392943388d03/,
     &                   p1(5)/-8.953073589022869d02/,
     &                   p1(6)/-9.229503102917111d01/,
     &                   p1(7)/-1.940352203312667d00/,
     &                   p1(8)/4.368019694395194d00/,
     &                   p1(9)/1.279153645893113d02/
      data               q1(1)/6.689575153359349d02/,
     &                   q1(2)/2.419887329355996d03/,
     &                   q1(3)/3.354196974608081d03/,
     &                   q1(4)/1.860416170944268d03/,
     &                   q1(5)/3.944307810159532d02/,
     &                   q1(6)/2.682132440551618d01/,
     &                   q1(7)/3.440812622259858d-01/,
     &                   q1(8)/5.948212550303777d01/
c
c ##coefficients for minimax approximation to ln(gamma(x)),
c ##1.5 .le. x .le. 4.0
c
      data               p2(1)/1.071722590306920d04/,
     &                   p2(2)/6.527047912184606d04/,
     &                   p2(3)/1.176398389569621d05/,
     &                   p2(4)/-9.726314581896472d03/,
     &                   p2(5)/-1.266094622188023d05/,
     &                   p2(6)/-5.393468741199669d04/,
     &                   p2(7)/-3.895916159676326d03/,
     &                   p2(8)/5.397392180667399d00/,
     &                   p2(9)/5.334390026324024d02/
      data               q2(1)/5.314589562326176d03/,
     &                   q2(2)/5.493654949398033d04/,
     &                   q2(3)/2.205757574602192d05/,
     &                   q2(4)/3.602313576600391d05/,
     &                   q2(5)/2.273446951911101d05/,
     &                   q2(6)/4.560612434396495d04/,
     &                   q2(7)/1.702062439974796d03/,
     &                   q2(8)/1.670328399370593d02/
c
c ##coefficients for minimax approximation to ln(gamma(x)),
c ##4.0 .le. x .le. 12.0
c
      data               p3(1)/-6.114039864945718d07/,
     &                   p3(2)/-5.588132821261888d08/,
     &                   p3(3)/-9.078970022444525d08/,
     &                   p3(4)/3.662935130796460d09/,
     &                   p3(5)/4.658700336821218d09/,
     &                   p3(6)/-4.570725249206307d09/,
     &                   p3(7)/-2.220833171087439d09/,
     &                   p3(8)/-1.474322990113017d04/,
     &                   p3(9)/-1.959850795570400d06/
      data               q3(1)/-5.636057205056241d05/,
     &                   q3(2)/-2.691827587118628d07/,
     &                   q3(3)/-4.411606716771217d08/,
     &                   q3(4)/-2.774890551941383d09/,
     &                   q3(5)/-6.579874397740792d09/,
     &                   q3(6)/-4.980644951174248d09/,
     &                   q3(7)/-6.677373781427094d08/,
     &                   q3(8)/-2.722530175870899d03/
c
c ##coefficients for minimax approximation to ln(gamma(x)),
c ##12.0 .le. x
c
      data               p4(1)/8.40596949829d-04/,
     &                   p4(2)/-5.9523334141881d-04/,
     &                   p4(3)/7.9365078409227d-04/,
     &                   p4(4)/-2.777777777769526d-03/,
     &                   p4(5)/8.333333333333333d-02/,
     &                   p4(6)/9.189385332046727d-01/,
     &                   p4(7)/-1.7816839846d-03/
      data               pi/3.141592653589793d0/
*-----------------------------------------------------------------------
      t   = x
c
c ##in this version of DLNGAM we don't want to mess with negative numbers
c
      if      (x .lt. rmin   ) then
        y = -rmax
      else if (x .lt. 1.5d00) then
        if (t .lt. 0.5d0 ) then
          b = -dlog(t)
          a = t
          t = t + one
        else
          b   = 0.0d0
          a   = t - one
        end if
        top = p1(8)*t + p1(9)
        den = t       + q1(8)
        do 100 j = 1, 7
          top = top*t + p1(j)
          den = den*t + q1(j)
100     continue
        y = (top/den)*a + b
c
      else if (x .le. 4.0d00) then
        b   = t - one
        top = p2(8)*t + p2(9)
        den = t + q2(8)
        a   = b - one
        do 200 j = 1, 7
          top = top*t+p2(j)
          den = den*t+q2(j)
200     continue
        y = (top/den)*a
c
      else if (x .le. 12.d00) then
        top = p3(8)*t + p3(9)
        den = t + q3(8)
        do 300 j = 1, 7
          top = top*t + p3(j)
          den = den*t + q3(j)
300      continue
         y = top/den
c
c ##else x .gt. 12, use asymptotic expansion
c
      else
        top = dlog(t)
        top = t*(top - one) - .5d0*top
        t   = one/t
        y   = top
        if (t .ge. rmin) then 
          b = t*t
          a = p4(7)
          do 400 j = 1, 5
            a = a*b + p4(j)
400       continue
          y = a*t+p4(6)+top
        end if
      end if
      DLNGAM = y
c
      return
c
*___end of DLNGAM___*
c
      end
      subroutine MINERR(fail,   inform, istate, iterqp, iwork,  kx,   
     &                  lenw,   mclin,  mdim,   msglev, nclin,  ndim,   
     &                  nlog,   nrankr, opchr,  ainp,   b,      bl,     
     &                  bu,     c,      clamda, cvec,   f,      obj,    
     &                  ssq1,   work)
      integer           fail,   inform, istate, iterqp, iwork,  kx,    
     &                  lenw,   mclin,  mdim,   msglev, nclin,  ndim,   
     &                  nlog,   nrankr, opchr
      double precision  ainp,   b,      bl,     bu,     c,      clamda, 
     &                  cvec,   f,      obj,    ssq1,   work
      dimension         fail(*), istate(*), iwork(*), kx(*), msglev(*), 
     &                  nlog(*), opchr(*), ainp(*), b(*), bl(*), bu(*), 
     &                  c(*), clamda(*), cvec(*), f(*), work(*)
************************************************************************
*
* MINERR computes the least squares solution when opchr(8) is positive
* and assumes the minimum error is zero otherwise.
*
* LAST MODIFIED: 4 FEB 1988
************************************************************************
      double precision  zero
      parameter        (zero = 0.d+0)
*-----------------------------------------------------------------------
      if ((opchr(8) .gt. 0) .and. (fail(1) .le. 0)) then
        call SCRNM(28, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
c ##GUESS INITIAL f
c
        call DLOAD(ndim, 1.d+0, f, 1)
c
c ##SOLVE THE LEAST SQUARES PROBLEM
c
        call LSOPTN('Cold Start')
        call LSOPTN('Problem type = ls3')
        call LSSOL (fail,   nrankr, ndim,   nclin,  mclin, mdim,   
     &              c,      bl,     bu,     cvec,   istate, kx,     
     &              f,      ainp,   b,      inform, iterqp, obj,    
     &              clamda, iwork,  mdim,   work,   lenw,   ssq1 )
        if ((inform .eq. 2) .or. (inform .eq. 3) .or. 
     &     (inform .eq. 6)) then
          call ERROR(95, nlog(2), msglev(2), inform, 0.d+0, ' ', 
     &               fail(1))
        else if ((inform .eq. 4) .or. (inform .eq. 5)) then
          call WARNG(27, nlog(2), msglev(2), inform, 0.d+0, ' ')
        end if
      else if (fail(1) .le. 0) then
        call SCRNM(31, nlog(5), msglev(5), 0, 0.d+0, ' ')
        call SCRNM(31, nlog(2), msglev(2), 0, 0.d+0, ' ')
        obj =  zero
      end if
      ssq1  = zero
      return
c
*___end of MINERR___*
c
      end
      subroutine QPREP (fail,   inform, istate, iterqp, kx,     msglev, 
     &                  nclin,  ndim,   nlog,   clamda, smp,    weak)
      integer           fail,   inform, istate, iterqp, kx,     msglev, 
     &                  nclin,  ndim,   nlog
      double precision  clamda, smp
      logical           weak
      dimension         fail(*), istate(*), kx(*), msglev(*), nlog(*), 
     &                  clamda(*)
************************************************************************
c
c QPREP reports the results of the QP solver
c
************************************************************************
      if (fail(1) .gt. 0) return
      if (inform  .eq. 0) then
        weak = .false. 
      else if (inform .eq. 1) then
        call WARNG(22, nlog(2), msglev(2), inform, smp, ' ')
        weak = .true.
      else if ((inform .eq. 4) .or.  (inform .eq. 5)) then
        call WARNG(34, nlog(2), msglev(2), inform, 0.d+0, ' ')
      else 
        call ERROR(95, nlog(2), msglev(2), inform, 0.d+0, ' ', fail(1))
      end if
      return
c
*___end of QPREP___*
c
      end 
