      subroutine DINTG (fun,    iquit,  limit,  ans,    aa,     bb, 
     &                  epss,   xmax,   xmin)
      integer           iquit,  limit
      double precision  aa,     ans,    bb,     epss,   fun,    xmax,   
     &                  xmin
************************************************************************
*
* DINTG estimates the integral of fun between aa and bb using the adaptive
* gaussian method described by Ian Robinson in the Australian Computer
* Journal, Vol. 8, Nov. 1976.
*
* VARIABLES(type, input/output, array)
* ans   :estimate of the integral
* xmax  :the maximum value of the function
* epss  :the target relative accuracy of ans
* aa, bb:the lower and upper interval of integration
* limit :maximum number of function evaluations
* xmin  :the smallest positive machine number
* iquit := 0  if convergence is achieved  using less than mdiv
*                levels of subdivision and less than 'limit'
*                function evaluations
*        = 1  if mdiv levels  of  subdivisions  are required at 
*                some stage of the routine (in this case, in-
*                vestigation of the particular subinterval is
*                terminated and the procedure continues on to
*                the next appropriate subinterval)
*        = 2  if more than limit function evaluations are re-
*                quired ( in  this case  the best estimate is
*                returned)
*
* LAST MODIFIED:10 NOV 1987
************************************************************************
      integer           i,      j,      jj,     k,      kk,     l,      
     &                  m,      mdiv,   nf,     npoint
      double precision  a,      arg,    at,     best,   bt,     dif,    
     &                  eps,    est,    est0,   estt1,  estt2,  
     &                  estt3,  f,      ft1,    ft2,    ft3,    ft4,    
     &                  ft5,    ft6,    p,      q,      r,      r1,     
     &                  r2,     r3,     s,      sum,    sum0,   t,      
     &                  two,    u3,     w1
      parameter        (mdiv   = 40)
      parameter        (npoint = 3)
      parameter        (two    = 2.d+0)
      dimension         j(mdiv), kk(mdiv), a(4,mdiv), at(npoint,mdiv), 
     &                  bt(npoint,mdiv), dif(npoint,mdiv), eps(mdiv), 
     &                  est(npoint,mdiv), estt1(npoint,mdiv), 
     &                  estt2(npoint,mdiv), estt3(npoint,mdiv), 
     &                  f(9,mdiv), ft1(npoint,mdiv), ft2(npoint,mdiv), 
     &                  ft3(npoint,mdiv), ft4(npoint,mdiv), 
     &                  ft5(npoint,mdiv), ft6(npoint,mdiv), 
     &                  sum(npoint,mdiv)
      external          fun
      data              p /4.44444444444444444444D-1/
      data              q /2.77777777777777777777D-1/
      data              r /2.25403330758516622964D-1/
      data              s /1.51278584715470076567D-2/
      data              t /9.84872141528452992343D-1/
*-----------------------------------------------------------------------
      xmax   = xmin
      iquit  = 0
      i      = 1
      ans    = xmin
      kk(1)  = 1
      j(1)   = 5
c
c ##INITIAL 3-POINT ESTIMATE
c
      a(1,1) = aa
      a(4,1) = bb
      r1     = bb - aa
      r2     = r*r1/two
      arg    = aa + r2
      f(2,1) = fun(arg)
      f(5,1) = fun((aa + bb)/two)
      arg    = bb - r2
      f(8,1) = fun(arg)
      est0   = r1*p*(625.d-03*(f(2,1) + f(8,1)) + f(5,1))
      eps(1) = epss
c
c ##PERFORM FIRST SUBDIVISION AND FORM 3 ESTIMATES AT LEVEL 1
c
      r2     = 2.d-01*r1
      r1     = r*r1
      a(2,1) = a(1,1) + r1
      a(3,1) = a(4,1) - r1
      w1     = a(1,1) + r2
      u3     = a(4,1) - r2
      arg    = a(2,1) - r2
      f(1,1) = fun(arg)
      f(3,1) = fun(w1)
      arg    = two*w1 - (a(1,1) + a(2,1))/two
      f(4,1) = fun(arg)
      arg    = two*u3 - (a(3,1) + a(4,1))/two
      f(6,1) = fun(arg)
      f(7,1) = fun(u3)
      arg    = a(3,1) + r2
      f(9,1) = fun(arg)
      nf     = 9
      xmax   = dmax1(f(1,1), f(2,1), f(3,1), f(4,1), f(5,1), f(6,1), 
     &               f(7,1), f(8,1), f(9,1), xmax)
c
      est(1,1) = r1*(q*(f(1,1) + f(3,1)) + p*f(2,1))
      est(2,1) = (a(3,1) - a(2,1))*(q*(f(4,1) + f(6,1)) + p*f(5,1))
      est(3,1) = r1*(q*(f(7,1) + f(9,1)) + p*f(8,1))
      sum0     = (est(1,1) + est(2,1) + est(3,1) - s*est0)/t
      best     = dabs(sum0)
      dif(1,1) = dabs(sum0 - est0)
      if (dif(1,1) .le. eps(1)*best) then
        ans = sum0
        return
      end if
c
c ##SUBDIVIDE AT LEVEL I AND FORM NEW GAUSSIAN SUMS ON EACH SUBINTERVAL
c
 1    i = i + 1
      if (nf .gt. limit) go to 15
      eps(i) = eps(i - 1)/1.2d00
      best   =best - dabs(sum0)
      l  = 1
      r2 = a(2,i - 1) - a(1,i - 1)
      r1 = r*r2
      r2 = 2.d - 1*r2
 2    at(l,i)  = a(l,i - 1)   +r1
      bt(l,i)  = a(l+1,i - 1) -r1
      r3  = bt(l,i) - at(l,i)
      w1  = a(l,i - 1)   + r2
      u3  = a(l+1,i - 1) - r2
      arg = at(l,i)    - r2
      ft1(l,i) = fun(arg)
      ft2(l,i) = fun(w1)
      arg      = two*w1 - (a(l,i - 1) + at(l,i))/two
      ft3(l,i) = fun(arg)
      arg      = two*u3 - (a(l+1,i-1) + bt(l,i))/two
      ft4(l,i) = fun(arg)
      ft5(l,i) = fun(u3)
      arg      = bt(l,i) + r2
      ft6(l,i) = fun(arg)
      xmax     = dmax1(ft1(l,i), ft2(l,i), ft3(l,i), ft4(l,i), ft5(l,i),
     &                ft6(l,i), xmax)
c
      estt1(l,i) = r1*(q*(ft1(l,i) + ft2(l,i)) + p*f(3*l-2,i-1))
      estt2(l,i) = r3*(q*(ft3(l,i) + ft4(l,i)) + p*f(3*l-1,i-1))
      estt3(l,i) = r1*(q*(ft5(l,i) + ft6(l,i)) + p*f(3*l,i-1))
      sum(l,i)   = (estt1(l,i) + estt2(l,i) + estt3(l,i)-s*est(l,i-1))/t
      dif(l,i)   = abs(sum(l,i) - est(l,i-1))
      best       = best + abs(sum(l,i))
      l = l + 2
      if (l - 4) 2, 4, 3
 3    l = 2  
      r2 = a(3,i-1) - a(2,i-1)
      r1 = r*r2
      r2 = 2.d-01*r2
      go to 2
 4    nf = nf+18
c
c ##DETERMINE ORDER OF PROCESSING ACCORDING TO MAGNITUDE OF ERROR
c ##ON EACH SUBINTERVAL
c
      k  = 1
      if (dif(2,i) .gt. dif(1,i)) k = 2
      if (dif(3,i) .gt. dif(k,i)) k = 3
      jj = k
      if (jj - 2) 6, 7, 8
 6    l = 2
      m = 3
      go to 9
 7    l = 1
      m = 3
      go to 9
 8    l     = 2
      m     = 1
 9    kk(i) = jj
      j(i)  = l
      if (dif(m,i) - dif(l,i)) 12, 10, 11
10    if (dif(m,i) .lt. -xmin) j(i) = 4
      go to 12
11    j(i) = m
c
c ##TEST FOR CONVERGENCE
c
12    if (dif(jj,i) .le. eps(i)*best) go to 14
      if (i .lt. mdiv) go to 13
      iquit = 1
      go to 14
c
c ##IF NO CONVERGENCE, DEFINE POINTS FOR USE IN NEXT SUBDIVISION
c
13    a(1,i) = a(jj,i-1)
      a(2,i) = at(jj,i)
      a(3,i) = bt(jj,i)
      a(4,i) = a(jj+1,i-1)
      f(1,i) = ft1(jj,i)
      f(2,i) = f(3*jj-2,i-1)
      f(3,i) = ft2(jj,i)
      f(4,i) = ft3(jj,i)
      f(5,i) = f(3*jj-1,i-1)
      f(6,i) = ft4(jj,i)
      f(7,i) = ft5(jj,i)
      f(8,i) = f(3*jj,i-1)
      f(9,i) = ft6(jj,i)
      est(1,i) = estt1(jj,i)
      est(2,i) = estt2(jj,i)
      est(3,i) = estt3(jj,i)
      sum0     = sum(jj,i)
      go to 1
c
c ##WRAP UP INTERVAL AFTER ALL 3 SUBINTERVALS AT LEVEL I HAVE CONVERGED
c
 14   ans = ans + sum(1,i) + sum(2,i) + sum(3,i)
      i = i - 1
      k = kk(i)
      sum(k,i) = xmin 
c
c ##PROCEED TO NEXT APPROPRIATE LEVEL
c
      jj = j(i)
      k  = kk(i)
      dif(k,i) = -1.0
      go to (6, 7, 8, 14, 16), jj
c
c ##IF FAIL EXIT
c
15    iquit = 2
c
c ##NORMAL EXIT
c
16    return
c
*___end of DINTG___*
c
      end
      function DBRENT  (ier,    fun,    xlow,   xhigh,  tolrt,  maxitr, 
     &                  wmach)
************************************************************************
c
c DBRENT is a root finding algorithm described in detail in
c Numerical Recipies by Press, Flannery, Teukolsy, and Vetterling.
c
c
************************************************************************
      integer           iter,   maxitr, ier
      double precision  DBRENT, xlow,   xhigh,  tolrt,  wmach,  fa,     
     &                  fb,     fc,     a,      tol1,   rmin,  rmax,  
     &                  fun,    b,      e,      d,      c,     two,    
     &                  zero
      parameter        (two  = 2.d+0)
      parameter        (zero = 0.d+0)
      dimension         wmach(*)
c-----------------------------------------------------------------------
      rmax = wmach(7)
      rmin = wmach(5)
      fa   = fun(xlow)
      fb   = fun(xhigh)
c
c ##first make sure the root is bounded
c
      if (fb*fa .gt. zero) then
        DBRENT = -rmax
        ier    = 1
        return
      end if
      a  = xlow
      b  = xhigh
      fc = fb
      do 100 iter = 1, maxitr
c
c ##rename variables and adjust interval containing root
c
        call DBADJ (fb, fa, fc, b, a, c, e, d)
c
c ##check for convergence
c
        tol1 = two*rmin*dabs(b) + tolrt/two
        if ((dabs((c-b)/two) .le. tol1) .or. (fb .eq. zero)) then
          DBRENT = b
          return
        end if
c
c ##find best estimate of root
c
        call DBRCOR(fb, fa, fc, b, a, c, e, d, tol1)
        fb = fun(b)
100   continue
c
c ##maximum number of iterations exceeded
c
        ier = 2
      DBRENT = b
      return
c
*___end of DBRENT___*
c
      end
      subroutine DBADJ (fb, fa, fc, b, a, c, e, d)
************************************************************************
c
c This subroutine renames varibles and adjusts the interval for the 
c root finding algorithm.
c
************************************************************************
      double precision fb, fa, fc, b, a, c, e, d
c
      if (fb*fc .gt. 0.d00) then
        c  = a
        fc = fa
        d  = b-a
        e  = d
      end if
      if (dabs(fc) .lt. dabs(fb)) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      end if
      return
c
*___end of DBADJ___*
c
      end
      subroutine DBRCOR(fb, fa, fc, b, a, c, e, d, tol1)
************************************************************************
c
c DBRCOR is the core of DBRENT, the root finding algorithm described 
c in detail in Numerical Recipies by Press, Flannery, Teukolsy, 
c and Vetterling.
c
************************************************************************
      double precision fb, fa, fc, b, a, c, e, d, tol1
      double precision one,    two,   xm, p,  q,  r, s
      parameter        (one = 1.d+0)
      parameter        (two = 2.d+0)
c
      xm = (c - b)/two
      if((dabs(e) .ge. tol1) .and. (dabs(fa) .gt. dabs(fb))) then
c
c ##attempt inverse quadratic interpolation
c
        s = fb/fa
        if (a .eq. c) then
          p = two * xm * s
          q = one - s
        else
          q = fa/fc
          r = fb/fc
          p = s * (two * xm * q * (q - r) - (b - a) * (r - one))
          q = (q - one) * (r - one) * (s - one)
        end if
        if (p .gt. 0.d00) q = -q
        p = dabs(p)
c
c ##check if inverse quadratic interpolation falls within interval
c
        if (two*p .lt. dmin1(3.d00*xm*q - dabs(tol1*q), dabs(e*q))) then
          e = d
          d = p/q
        else
c
c ##interpolation failed, use bisection
c
          d = xm
          e = d
        end if
      else
c
c ##bounds decreasing too slowly, use bisection
c
        d = xm
        e = d
      end if
c
c ##move last best guess to a
c
      a  = b
      fa = fb
      if(dabs(d) .gt. tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      end if
c
      return
c
*___End of DBRCOR___*
c
      end
      function    BRENT(ax,     bx,    cx,     F,      tol,    xmin,   
     &                  z,      inform)
      integer           inform
      double precision  ax,     BRENT, bx,     cx,     F,     tol,    
     &                  xmin,   z
      dimension         z(*)
************************************************************************
*
* BRENT minimizes a function (F(z)) of 1 variable.  Taken from 
* NUMERICAL RECIPIES, Press et. al., 1986.
*
* inputs:
* ax, bx, cx: abscissas, bx is between ax and cx, and f(bx) is less than
*             both f(ax) and f(cx).
* z         : a vector that is passed the the function being minimized.
*
* outputs:
* inform    = -1 if the maximum number of evaluations were used.
* 
* LAST MODIFIED: 5 DEC 1988
************************************************************************
      integer           iter,   itmax
      double precision  two,    x,      xm,     y
      parameter        (itmax = 100)
      parameter        (two   = 2.d+0)
      dimension         x(8), y(4)
*
* Refer to numerical recipies:
*         x(1:8) = u, v, w, x, a, b, d, e
*         y(1:4) = fu, fv, fw, fx
c-----------------------------------------------------------------------
      inform = 1
      x(5)   = dmin1(ax,cx)
      x(6)   = dmax1(ax,cx)
      x(2)   = bx
      x(3)   = x(2)
      x(4)   = x(2)
      x(8)   = 0.d+0
      z(1)   = x(4)
      y(4)   = F(z)
      y(2)   = y(4)
      y(3)   = y(4)
      do 100 iter = 1, itmax
        xm   = (x(5) + x(6)) / two
c
c ##CHECK AND SEE IF WE ARE DONE
c
        if (dabs(x(4) - xm) .le. 
     &     (two*tol - (x(6) - x(5))/two)) goto 1000
c
c ##FIND THE NEXT ITERATE
c
        call DBMCOR(x,   y,   tol)
        z(1) = x(1)
        y(1) = F(z)
c
c ##DECIDE WHAT TO DO WITH THE FUNCTION EVALUATION, BOOK KEEPING
c ##FOLLOWS
c
        call DBMORD(x, y)
        inform = inform + 1
100   continue
      inform = -1
1000  xmin   = x(4)
      BRENT  = y(4)
      return
c
*___end of BRENT___*
c
      end
      subroutine DBMORD(x,      y)
      double precision  x,      y
      dimension         x(*), y(*)
************************************************************************
* DBMORD performs housekeeping for BRENT.
*
* LAST MODIFIED:12 DEC 1988
************************************************************************
      if (y(1) .le. y(4)) then
        if (x(1) .ge. x(4)) then
          x(5) = x(4)
        else
          x(6) = x(4)
        end if
        x(2) = x(3)
        y(2) = y(3)
        x(3) = x(4)
        y(3) = y(4)
        x(4) = x(1)
        y(4) = y(1)
      else
        if (x(1) .lt. x(4)) then
          x(5) = x(1)
        else
          x(6) = x(1)
        end if
        if ((y(1) .le. y(3)) .or. (x(3) .eq. x(4))) then
          x(2) = x(3)
          y(2) = y(3)
          x(3) = x(1)
          y(3) = y(1)
        else if ((y(1) .le. y(2)) .or. (x(2) .eq. x(4)) .or. 
     &           (x(2) .eq. x(3))) then
          x(2) = x(1)
          y(2) = y(1)
        end if
      end if
      return
c
*___End of DBMORD___*
c
      end
      subroutine DBMCOR(x,   y,   tol)
      double precision  x,   y,   tol
      dimension         x(*), y(*)
************************************************************************
* DBMCOR determines the next iterate to use in the minimization routine
* BRENT. 
*
* LAST MODIFIED:12 DEC 1988
************************************************************************
      double precision  cgold,  etemp,  p,      q,      r,      two,    
     &                  xm
      logical           golden
      parameter        (cgold = 0.3819660d+0)
      parameter        (two   = 2.d+0)
c-----------------------------------------------------------------------
      xm = (x(5) + x(6)) / two
      if (dabs(x(8)) .gt. tol) then
        golden = .false.
c
c ##CONSTRUCT A TRIAL PARABOLIC FIT
c
        r = (x(4) - x(3)) * (y(4) - y(2))
        q = (x(4) - x(2)) * (y(4) - y(3))
        p = (x(4) - x(2)) * q - (x(4) - x(3)) * r
        q = two * (q - r)
        if (q .gt. 0.d+0) p = -p
        q     = dabs(q)
        etemp = x(8)
        x(8)     = x(7)
        if ( (dabs(p) .ge. dabs(q * etemp / two)) .or. 
     &       (p .le. q * (x(5) - x(4))) .or. 
     &       (p .ge. q * (x(6) - x(4)))) then
          golden = .true.
        else
c
c ##TAKE THE PARABOLIC STEP
c
          x(7) = p / q
          x(1) = x(4) + x(7)
          if ((x(1) - x(5) .lt. two*tol) .or. 
     &        (x(6) - x(1) .lt. two*tol)) 
     &         x(7) = dsign(tol,xm - x(4))
        end if
      else
        golden = .true.
      end if
c
c ##GOLDEN SECTION
c
      if (golden) then
        if (x(4) .ge. xm) then
          x(8) = x(5) - x(4)
        else
          x(8) = x(6) - x(4)
        end if
        x(7) = cgold * x(8)
      end if
      if (dabs(x(7)) .ge. tol) then
        x(1) = x(4) + x(7)
      else
        x(1) = x(4) + dsign(tol,x(7))
      end if
      return
c
*___end of DBMCOR___*
c
      end
