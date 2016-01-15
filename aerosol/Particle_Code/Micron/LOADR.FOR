      subroutine LOADR (flag,   mrho,   n,      rho,    diam)
      integer           flag,   mrho,   n
      double precision  diam,   rho
      dimension         diam(*), rho(mrho,*)
************************************************************************
* LOADR loads the covariance matrix and must be maintained by each
* user.
*
* LAST MODIFIED:17 OCT 1989
************************************************************************
      integer           i,      ilen,   j
      double precision  delta,  sdev,   TSTFN
*-----------------------------------------------------------------------
c
c ##SIZE DISTRIBUTION IS UNCORRELATED IF THE SIZE RATIO IS
c ##LESS THAN \SQRT{10}
c
      delta = diam(2) - diam(1)
      ilen  = nint ( dlog(10.d+0) / delta / 2.d+0 )
      if (flag .le. 0) then
        do 100 i = 1, n
          sdev = TSTFN(diam(i))
          rho(i,i) = sdev * sdev 
100     continue
        do 200 i = 1, n
          do 110 j = 1, i-1
            if ((i-j) .ge. ilen) then
              rho(i,j) = 0.d+0
            else
              rho(i,j) = dsqrt( rho(i,i)*rho(j,j) ) *
     &                  (1.d+0 - dble(i-j) / dble(ilen))
            end if
            rho(j,i) = rho(i,j)
110       continue
200     continue
      end if
      return
*
*___end of LOADR___*
*
      end
