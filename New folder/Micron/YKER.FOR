      subroutine KERINF(fail,   maxfev, maxtok, mfail,  minfev, minst,  
     &                  mintg,  mpar,   mrecl,  msgscr, nch,    npar,   
     &                  nscr,   decdiv, dspc,   par,    tol,    wmach,  
     &                  xkrvec, cfiles, inme,   parnm,  tokens)
      integer           fail,   maxfev, maxtok, mfail,  minfev, minst,  
     &                  mintg,  mpar,   mrecl,  msgscr, nch,    npar,   
     &                  nscr
      double precision  decdiv, dspc,   par,   tol,    wmach,  xkrvec
      character*(*)     cfiles, inme,   parnm,  tokens
      dimension         fail(*), nch(*), npar(*), dspc(*), par(minst,*), 
     &                  tol(*), wmach(*), xkrvec(*), cfiles(*), 
     &                  inme(*), parnm(*), tokens(*)
************************************************************************
*
* KERINF creates a data file (micron.kri) which contains the following 
* information about the kernel functions:
* 1  the number of unconnected intervals which make up the support,
*    their location, boundaries, integral, and maximum height.
* 2  the log-mean diameter
* 3  integral
* 4  maximum height
* 5  variance
*
* INPUTS:
*
* LAST MODIFIED:30 JAN 1989
************************************************************************
      integer           ibnd,   ibndrc, icha,   ichrec, iinst,  ioerr,  
     &                  ipar,   ISUM,   itmpr,  ncha,   ncht,   nfev,   
     &                  ninf
      double precision  dadd,   dmin,   dmax
      character*8       name
      common /cmkerf/   iinst,  icha,   ncha
      save   /cmkerf/
c-----------------------------------------------------------------------
      call SCRNM(15, nscr, msgscr, 0, 0.d+0, ' ')
      ninf = 0
      dmin = dspc(3)
      dmax = dspc(4)
c
c ##initialize the pointers used in the micron.kri file
c
      call OPNKRI(ioerr, 3, mrecl, ninf, cfiles)
      if (ioerr .ne. 0) go to 1000
      ncht   = ISUM(minst, nch, 1)
      ichrec = minst*3
      ibndrc = ichrec + 2*ncht
c
c ##and determine a reasonable number of function evaluations
c 
      call FNFEN (maxfev, minfev, nfev, dadd, decdiv, dmax, dmin)
c
c ##begin looking at all the instruments
c
      do 200 iinst = 1, minst
      if (fail(mfail) .le. 0) then
c
c ##here we make sure an instrument with no channels points at rec=0
c 
        if (nch(iinst) .gt. 0) then
          itmpr = ichrec + 1
        else
          itmpr = 0
        end if
        if (iinst .eq. 1) then
          write(ninf, rec = iinst)
     &      itmpr, nch(iinst), inme(iinst), npar(iinst), minst, mrecl
        else
          write(ninf, rec = iinst) 
     &      itmpr, nch(iinst), inme(iinst), npar(iinst)
        end if
        if (npar(iinst) .gt. 0) then
          write(ninf, rec = iinst + minst) 
     &    (par(iinst,ipar), ipar = 1, npar(iinst))
          write(ninf, rec = iinst + minst + minst) 
     &      (parnm(ipar), ipar = iinst, (npar(iinst) - 1)*minst + iinst, 
     &      minst)
        end if
c
c ##begin looking at each channel in the instrument
c
        ncha = nch(iinst)
        name = inme(iinst)
        do 210 icha = 1, ncha
          ibnd = 0
          call CHNINF(fail,   ibnd,   ibndrc, ichrec, icha,   maxfev, 
     &                mfail,  mintg,  msgscr, nfev,   ninf,   nscr, 
     &                dadd,   dmin,   tol,    wmach,  xkrvec, name)
          ibndrc = ibndrc + ibnd
210     continue
      end if
200   continue
c
1000  if (ioerr .ne. 0) then
        call ERROR(42, nscr, msgscr, ioerr, 0.d+0, 'open', 
     &             fail(mfail))
      else if (ninf .ne. 0) then
        if (fail(mfail) .gt. 0) 
     &    write(ninf, rec = 1) fail(mfail)
        close(ninf)
      end if
      return
*
*___end of KERINF___*
*
      end
      subroutine FNFEN (maxfev,  minfev, nfev,   dadd,   decdiv, dmax, 
     &                  dmin)
************************************************************************
*
* FNFEN is used by KERINF and UPDATE to find a reasonable number of 
* divisions for the interval (dmin,dmax).  The output number, nfev,
* is used to estimate the maximum of the kernel function.
* Two items to watch out for when choosing nfev is 
* 1) too large an nfev wastes storage and computation time
* 2) too small an nfev can cause a interval of kernel function support
*    to be entirely missed.
*
* INPUTS:
* maxfev - the maximum number of divisions for (dmin,dmax)
* minfev - the minimum number of divisions for (dmin,dmax)
* decdiv - the desired number of divisions per decade
*
* OUTPUTS:
* nfev   - the number of divisions for (dmin,dmax)
*
* LAST MODIFIED: 29 OCT 1987
************************************************************************
      integer           maxfev,  minfev, nfev
      double precision  bse10,   dadd,    decdiv, dmax,   dmin,   xchop
      parameter        (bse10  = 2.303d00)
c
      xchop = (dmax - dmin)*decdiv*bse10
      nfev  = dint(xchop)
      if (nfev .lt. minfev) then
        nfev = minfev
      else if (nfev .gt. maxfev) then
        nfev = maxfev
      end if
      dadd  = (dmax - dmin)/dble(nfev - 1)
      return
c
*___end of FNFEN___*
c
      end
      subroutine CHNINF(fail,   ibnd,   ibndrc, ichrec, icha,   maxfev, 
     &                  mfail,  mintg,  msgscr, nfev,   ninf,   nscr, 
     &                  dadd,   dmin,   tol,    wmach,  xkrvec, name)
************************************************************************
*
* CHNINF determines the intervals of a kernel function, and 
* finds the integral, variance, maximum, and mean diameter of the
* kernel function over each
* interval.  Additionally these properties are found for the kernel
* function over the interval (dmin, dmax).  The results of the 
* calculation are written to the data file micron.kri.
*
* LAST MODIFIED:18 NOV 1987
* bug: 'name' dimension should be adjustable
************************************************************************
      integer           fail,   ibnd,   ibndrc, icha,   ichrec, IDAMAX, 
     &                  ierr,   ierx,   iexx,   iker,   imax,   jker,   
     &                  maxfev, mfail,  mintg,  msgscr, nfev,   ninf,   
     &                  nscr 
      double precision  dadd,   diam,   djerr,  dlow,   dmin,   dup,    
     &                  ediam,  garb,   KERX,   parea,  pmax,   pmean,  
     &                  pvar,   rmin,   tarea,  tmean,  tol,    tvar,   
     &                  wmach,  XKERX,  xkrvec, xmax,   xnfev,  XXKERX
      character*8       name
      logical           found,  pos
      dimension         tol(*), wmach(*), xkrvec(maxfev), fail(*)
      external          KERX, XKERX, XXKERX
*-----------------------------------------------------------------------
      if (fail(mfail) .gt. 0) return
      call SCRNM(16, nscr, msgscr, icha, 0.d+0, name)
      rmin  = wmach(5)
      tmean = 0.d00
      tvar  = 0.d00
      tarea = 0.d00
      diam  = dmin
c
c ##now we represent the kernel function as a vector, and find an
c ##approximate maximum
c 
      xnfev = dble(nfev)
      do 100 iker = 1, nfev
        xkrvec(iker) = KERX(diam)
        diam         = diam + dadd
c
c ##  CHECK IF THE KERNEL FUNCTION HAS NEGATIVE RESPONSE
c
        if ((-rmin*10.d00) .gt. xkrvec(iker)) then
          ediam  = dexp(diam)
          call ERROR(20, nscr, msgscr, icha, ediam, 
     &              name, fail(mfail))
          return
        end if
100   continue
      imax = IDAMAX(nfev, xkrvec, 1)
      xmax = xkrvec(imax)
c
c ##MAKE SURE THE INSTRUMENT RESPONSE IS REASONABLE
c
      if (xmax .lt. (10.d00*rmin)) then
        call ERROR(44, nscr, msgscr, icha, xmax, name, fail(mfail))
        return
      end if
c
c ##before we find all the intervals of support, area's etc., we need
c ##to check if the support's boundary coincides with the maximum pos-
c ##sible inversion boundary, and we need to initialize the sums to 
c ##zero.
c 
      diam  = dmin
      found = .false.
      if ((xkrvec(1)/xmax) .lt. tol(4)) then
        pos = .false.
      else
        call WARNG(14, nscr, msgscr, icha, xkrvec(1), name)
        pos  = .true.
        dlow = dmin
      end if
c
c ##now we start looking for intervals where the kernel functions are 
c ##nonzero
c
      djerr = -0.5d+0
      do 200 jker = 2, nfev
        diam = diam + dadd
c
c ##if the kernel is nonzero we make sure we're not at the end of the
c ##maximum inversion interval.
c
        if ((xkrvec(jker)/xmax) .gt. tol(4)) then
          if (.not.pos) call REFINE(dadd, dlow, diam, tol, xmax, pos)
          if (jker .eq. nfev) then
            call WARNG(15, nscr, msgscr, icha, xkrvec(nfev), name)
            dup   = diam
            found = .true.
          end if
        else
c
c ##if the kernel is 'zero' we check if we're at the end of an interval
c ##and if so we sharpen the bound and set the integration flag
c
          if (pos) then
            found = .true.
            call REFINE(dadd, dup, diam, tol, xmax, pos)
          end if
        end if
c
c ##if we found an interval (dlow,dup) then we integrate the kernel
c ##function over this interval and record the results in the micron.kri
c ##data file, and add to the totals.
c
        if (found) then
          ibnd  = ibnd + 1
          found = .false.
          call DINTG(KERX,  ierr,  mintg, parea, dlow, dup, tol(6), 
     &               pmax, rmin)
          call DINTG(XKERX, ierx,  mintg, pmean, dlow, dup, tol(6), 
     &               garb, rmin)
          call DINTG(XXKERX, iexx, mintg, pvar,  dlow, dup, tol(6), 
     &               garb, rmin)
          if (ierr .ne. 0) then
            djerr = dble(ierr)
          else if (ierx .ne. 0) then
            djerr = dble(ierx)
          else if (iexx .ne. 0) then
            djerr = dble(iexx)
          end if
          if (xmax .lt. pmax) xmax = pmax
          tarea = tarea + parea
          tmean = tmean + pmean
          tvar  = tvar  + pvar
          pmean = pmean/parea
          pvar  = pvar/parea
          pvar  = pvar  - pmean*pmean
          write(ninf, rec = ibndrc + ibnd) dlow, dup, parea, pmax, 
     &      pmean, pvar
        end if
200   continue
      if (djerr .gt. 0.d+0) 
     & call WARNG (19, nscr, msgscr, icha, djerr, name)
c
c ##after we've checked the entire interval we record the totals
c
      tmean  = tmean/tarea
      tvar   = tvar/tarea
      tvar   = tvar - tmean*tmean
      ichrec = ichrec + 1
      write(ninf, rec = ichrec) ibndrc + 1, ibnd, tarea, xmax
      ichrec = ichrec + 1
      write(ninf, rec = ichrec) tmean, tvar
      return
c
*___end of CHNINF___*
c
      end
      subroutine REFINE(dadd,   dedge,  diam,   tol,    xmax,   pos)
************************************************************************
*
* Given two diameters which bound one end of an interval of kernel 
* support, REFINE returns a diameter which is just outside the interval 
* of support. 
*
* INPUTS:
* pos    - true on input if REFINE is finding the upper diameter of an
*          interval of support, and false otherwise.  The logical value
*          of pos is changed when REFINE is exited.
* diam   - the lower diameter of the subinterval to be refined
* dadd   - the logarithmic length of the subintervals
* xmax   - the estimated maximum of the kernel functions
* tol(4) - tolerance used to decide when the kernel function is zero
* tol(5) - the tolerance used to decide how accurately the edge of
*          the kernel function support needs to be found.
*
* OUTPUT:
* dedge  - the estimated edge of the support.  dedge lies just outside
*          the kernel function support.
*
* LAST MODIFIED: 30 OCT 1987
************************************************************************
      integer           iter,   niter
      double precision  dadd,   dedge,  diam,   dneg,   dpos,   dtemp,  
     &                  fnc,    KERX,   tol,    xiter,  xmax 
      logical           pos
      dimension         tol(*)
*-----------------------------------------------------------------------
c ##decide how many bisections are needed to refine the subinterval
c
      xiter = dlog(dadd/tol(5))/dlog(2.d00)
      niter = int(xiter) + 1
c
      if (pos) then
        dpos = diam - dadd
        dneg = diam
        pos = .false.
      else
        dpos = diam
        dneg = diam - dadd
        pos = .true.
      end if
c
c ##perform bisections to cut the bounds down to the specified width
c
      do 100 iter = 1, niter
        dtemp = (dpos + dneg)/2.d00
        fnc = KERX(dtemp)
        if ((fnc/xmax) .gt. tol(4)) then
          dpos = dtemp
        else 
          dneg = dtemp
        end if
100   continue
      dedge = dneg
      return
c
*___end of REFINE___*
c
      end
      subroutine KERECH(fail,   mfail,  mpar,   mrecl,  msgscr, nscr, 
     &                  par,    cfiles, parnm)
      integer           fail,   mfail,  mpar,   mrecl,  msgscr, nscr
      character*(*)     cfiles, parnm
      dimension         cfiles(*), par(*), parnm(*), fail(*)
************************************************************************
c
c KERECH echos the results stored in the unformatted direct access data 
c file, micron.kri, and places them in the formatted data file kernel.out.
c
c LAST MODIFIED: 30 OCT 1987
c
************************************************************************
      integer           ibnd,   ibndpt, icha,   ichpnt, inst,   ioerr
      integer           ipar,   jbnd,   jchpnt, joerr,  minst
      integer           ncha,   nfor,   npar,   nunf
      double precision  area,   dlow,   dmean,  dup,    dvar,   par
      double precision  var,    xmax,   xmean 
      character*8       name
c
      call SCRNM(17, nscr, msgscr, 0, 0.d+0, ' ')
      call OPNKRI(ioerr, 2, mrecl, nunf,   cfiles)
      if (ioerr .ne. 0) go to 1000
      call FUNIT(nfor)
      open(nfor, file = cfiles(5), status = 'unknown')
      write(nfor, 8)
      read(nunf, rec = 1, err = 1000, iostat = joerr)
     &  ichpnt, ncha, name, npar, minst
      write(nfor, 13) minst
      do 100 inst = 1, minst
        read(nunf, rec = inst, err = 1000, iostat = joerr) 
     &    ichpnt, ncha, name, npar
c
c NEED to make sure npar .le. mpar
c
        if (mpar .lt. npar) npar = mpar
        write(nfor,1) name, ncha
        if (npar .gt. 0) then
          read(nunf, rec = minst + inst, err = 1000, iostat = joerr)
     &    (par(ipar), ipar = 1, npar)  
          read(nunf, rec = minst + minst + inst, err = 1000, 
     &         iostat = joerr) (parnm(ipar)(1:8), ipar = 1, npar)
          write(nfor, 5)
          do 120 ipar = 1, npar
            write(nfor, 6) parnm(ipar)(1:8), par(ipar)
120       continue
          else
          write(nfor, 11)
        end if
        do 110 icha = 1, ncha
          jchpnt = (icha - 1)*2 + ichpnt
          read(nunf, rec = jchpnt, err = 1000, iostat = joerr) 
     &    ibndpt, ibnd, area, xmax
          write(nfor, 2) icha, ibnd, area
          read(nunf, rec = jchpnt + 1, err = 1000, iostat = joerr) 
     &      xmean, var
          xmean = dexp(xmean)
          write(nfor, 3) xmax, xmean, var
          do 111 jbnd = ibndpt, ibndpt + ibnd - 1
            read(nunf, rec = jbnd, err = 1000, iostat = joerr) 
     &      dlow, dup, area, xmax, dmean, dvar
            dlow =  dexp(dlow)
            dup =   dexp(dup)
            dmean = dexp(dmean)
            write(nfor, 4) dlow, dup,   area
            write(nfor,12) xmax, dmean, dvar
111       continue
          write(nfor, 7)
110     continue
        write(nfor, 7)
100   continue
c
1000  if (ioerr .ne. 0) then
        call ERROR(42, nscr, msgscr, ioerr, 0.d+0, 'open', fail(mfail))
      else 
        close(nunf)
        close(nfor)
        if (joerr .ne. 0) 
     &    call ERROR(42, nscr, msgscr, joerr, 0.d+0, 'read', 
     &               fail(mfail))
      end if
      return
1     format(' Instrument: ', a8, '   Number of channels: ', i3,//)
2     format(' Channel: ',t12,i3,t28,'  # of intrvls: ',t43,i2,t58,
     &' integral=  ',t71,d10.3)
3     format(' maximum=  ',t12,d10.3,t28,'  ln-mean dia.= ',t43,d10.3,
     &t58,' variance= ',t71,d10.3,//,' Intervals: ')
4     format(' Lower dia.= ',t12,d10.3,t28,'  Upper dia.= ',t43,d10.3,
     &t58,' integral= ',t71,d10.3)
5     format(' Parameters:')
6     format(1x,a8,' : ',d10.3//)
7     format(//)
8     format(' This file echos the information contained in the kernel.i
     &nf data file.',///)
11    format(//)
12    format(' maximum=    ',t12,d10.3,t28,'  ln-mean dia.= ',t43,d10.3,
     &t58,' Variance= ',t71,d10.3,/)
13    format(/,' Number of instruments = ',i3,/)
c
*___end of KERECH___*
c
      end
      subroutine UPDATE(fail,   maxfev, mfail,  minfev, minst,  mintg, 
     &                  mpar,   mrecl,  msgscr, nch,    npar,   nscr,   
     &                  decdiv, dspc,   par,    tol,    wmach,  xkrvec, 
     &                  cfiles, inme,   parnm)
      integer           fail,   maxfev, mfail,  minfev, minst,  mintg,
     &                  mpar,   mrecl,  msgscr, nch,    npar,   nscr
      double precision  decdiv, dspc,   par,    tol,    wmach,  xkrvec
      character*(*)     cfiles, inme,   parnm
      dimension         fail(*), npar(*), nch(*), dspc(*), par(minst,*), 
     &                  tol(*), wmach(*), xkrvec(*), cfiles(*), inme(*), 
     &                  parnm(*)
************************************************************************
*
* UPDATE determines the following information about the kernel functions
* from a single instrument:
* 1  the number of unconnected intervals which make up the support,
*    their location, boundaries, integral, and maximum height.
* 2  the log-mean diameter
* 3  integral
* 4  maximum height
* 5  variance
*
* UPDATE is 
* similar to the subroutine KERINF except that UPDATE only perform the
* computation for one of the existing instruments and updates the 
* information contained in the micron.kri data file for this instrument. 
* UPDATE may be useful if one instrument depends strongly on a parameter
* and many inversions are to be performed with the parameter changing.
*
*
* LAST MODIFIED:14 FEB 1989
************************************************************************
      integer           ibnd,   ibndrc, icha,   ichrec, iinst,  ioerr,  
     &                  ipar,   ichsve, itrec,  ncha,   nfev,   ninf,   
     &                  ntmp
      double precision  dadd,   dmax,   dmin
      common /cmkerf/   iinst,  icha,   ncha
      save   /cmkerf/
*-----------------------------------------------------------------------
c
c ##first we need to find out what instrument we need to update and
c ##what the parameters are
c
      ichsve = icha
      dmin   = dspc(3)
      dmax   = dspc(4)
      call OPNKRI(ioerr, 2, mrecl, ninf,   cfiles)
      if (ioerr .ne. 0) go to 1000
      call SCRNM(20, nscr, msgscr, 0, 0.d00, inme(iinst)) 
c
c ##and initialize the pointers used in the scratch file
c
      call FUNIT(ntmp)
      open(ntmp, status = 'scratch', access = 'direct', 
     &     form = 'unformatted',
     &     recl = mrecl)
      ncha = nch(iinst)   
      ichrec = 3
      ibndrc = 3 + 2*ncha
c
c ##and determine a reasonable number of function evaluations
c 
      call FNFEN (maxfev, minfev, nfev, dadd, decdiv, dmax, dmin)
c
c ##begin looking at all the instruments
c
      if (nch(iinst) .gt. 0) then
        itrec = ichrec + 1
        write(ntmp, rec = 1) itrec, nch(iinst), inme(iinst), npar(iinst)
        if (npar(iinst) .gt. 0) then
          write(ntmp, rec = 2) 
     &      (par(iinst, ipar),   ipar = 1, npar(iinst))
          write(ntmp, rec = 3) 
     &      (parnm(ipar), ipar = iinst, (npar(iinst) - 1)*minst + iinst, 
     &      minst)
        end if
c
c ##begin looking at each channel in the instrument
c
        do 200 icha = 1, ncha
          if (fail(mfail) .gt. 0) go to 1000
          ibnd = 0
          call CHNINF(fail, ibnd,   ibndrc, ichrec, icha,   maxfev, 
     &              mfail,  mintg,  msgscr, nfev,   ntmp,   nscr, 
     &              dadd,   dmin,   tol,    wmach,  xkrvec, inme(iinst))
          ibndrc = ibndrc + ibnd
200     continue
      else
        write(ntmp, rec=1) 0, 0, inme(iinst), 0
      end if
      call SCRNM (1, nscr, msgscr, 0, 0.d+0, ' ')
      call CMBFIL(fail, iinst, mfail, minst, mpar, 
     &                 mrecl, ninf, ntmp, par, parnm)
c
c ##CLOSE THE APPROPRIATE FILES AND CALL THE NEEDED ERROR MESSAGES 
c ##BEFORE LEAVING
c
1000  if (ioerr .ne. 0) then
        call ERROR(42, nscr, msgscr, ioerr, 0.d+0, 'open', fail(mfail))
      else
        close(ninf)
        close(ntmp)
      end if
      icha = ichsve
      return
*
*___end of UPDATE___*
*
      end
      subroutine CMBFIL(fail,   inst,   mfail,  minst,  mpar,   mrecl, 
     &                  ninf,   ntmp,   par,    parnm)
************************************************************************
* CMBFIL is callled by UPDATE of replace all of the information in
* micron.kri concerning the instrument listed in position inst in 
* micron.dat.  The replacement information is read from the newly formed 
* scratch file (unit = ntmp).
*
*
* LAST MODIFIED:18 NOV 1987
*
************************************************************************
      integer           fail,   ichrec, inst,   jchrec, kbnrec, kchrec, 
     &                  mfail,  minst,  mpar,   mrecl,  nbll,   nbln,   
     &                  nbm,    nbpl,   nbpn,   nchdif, nchn,   ncho,   
     &                  ninf,   ntmp
      double precision  par
      character*8       name,   parnm
      dimension         par(*), parnm(*), fail(*)
*-----------------------------------------------------------------------
c
c ##here we find out how much we have to move the trailing set
c ##of records in unit ninf due to a difference in the number of
c ##intervals,  and move them
c
      if (fail(mfail) .gt. 0) return
      call FBNBL(1, 1, ntmp, nbpn, nbln)
      call FBNBL(inst, inst, ninf, nbpl, nbll)
      nbm = (nbln - nbpn) - (nbll - nbpl)
      if (nbpn .gt. 0) nbm = nbm + 1
      if (nbpl .gt. 0) nbm = nbm - 1
      call MVB(inst + 1, minst, nbm, ninf)
c
c ##we find out how much we have to move the records due to a difference
c ##in the number of channels
c
      read(ntmp, rec = 1)    jchrec, nchn, name
      read(ninf, rec = inst) ichrec, ncho
      nchdif = nchn - ncho
c
c ##and move them
c
      if (nchdif .gt. 0) then 
        call MVB(inst + 1, minst,    2*nchdif, ninf)
        call MVB(1,        inst - 1, 2*nchdif, ninf)
        call MVC(inst + 1, minst,    2*nchdif, ninf)
      else if (nchdif .lt. 0) then
        call MVC(inst + 1, minst,    2*nchdif, ninf)
        call MVB(1, inst - 1,        2*nchdif, ninf)
        call MVB(inst + 1, minst,    2*nchdif, ninf)
      end if
      if (nchn .le. 0) then
        write(ninf, rec = inst) 0, 0, name, 0
        go to 1000
      end if
c
c ##now we find the beginning records of the blocks containing
c ##channel and interval information
c
      call FNP(inst,   kbnrec, kchrec, minst,  ninf)
c
c ##and here we write in the new information
c
      call RNK(inst,   kbnrec, kchrec, minst,  mpar,   mrecl,  
     &         ntmp,   ninf,   par,    parnm)
1000  return
*
*___end of CMBFIL___*
*
      end
      subroutine RNK   (inst,   kbnrec, kchrec, minst,  mpar,   mrecl,  
     &                  ntmp,   ninf,   par,    parnm)
************************************************************************
*
* RNK reads the scratch file containing the updated information for the
* instrument listed in position inst in micron.dat, and writes the 
* information to the permanent data file micron.kri.
*
* NEED to worry about read and write errors possibly.
*
* LAST MODIFIED: 31 OCT 1988
************************************************************************
      integer           ibnd,   icha,   inst,   ipar,   jbnd,   jbndpt
      integer           jchrec, kbnrec, kchrec, minst,  mpar,   mrecl,  
     &                  nchn,   npar,   ntmp,   ninf 
      double precision  area,   dlow,   dup,    xmax,   xmean,  var
      double precision  par
      character*8       nnew,   parnm
      dimension         par(*), parnm(*)
*-----------------------------------------------------------------------
      read(ntmp,  rec = 1)    jchrec, nchn, nnew, npar
      if (inst .eq. 1) then
        write(ninf, rec = inst) kchrec, nchn, nnew, npar, minst, mrecl
      else
        write(ninf, rec = inst) kchrec, nchn, nnew, npar
      end if
      if (npar .gt. 0) then
        read(ntmp,  rec = 2)            
     &   (par(ipar), ipar = 1, 1 + minst*(npar - 1), minst)
        write(ninf, rec = inst + minst) 
     &   (par(ipar), ipar = 1, 1 + minst*(npar - 1), minst)
        read(ntmp,  rec = 3)            
     &   (parnm(ipar), ipar = 1, 1 + minst*(npar - 1), minst)
        write(ninf, rec = inst + minst + minst) 
     &   (parnm(ipar), ipar = 1, 1 + minst*(npar - 1), minst)
      end if
      do 100 icha = 1, nchn
        read(ntmp, rec  = jchrec) jbndpt, jbnd, area, xmax
        write(ninf, rec = kchrec) kbnrec, jbnd, area, xmax
        read(ntmp, rec  = jchrec + 1) xmean, var
        write(ninf, rec = kchrec + 1) xmean, var
        do 110 ibnd = 1, jbnd
          read(ntmp,  rec = jbndpt) dlow, dup, area, xmax, xmean, var
          write(ninf, rec = kbnrec) dlow, dup, area, xmax, xmean, var
          kbnrec = kbnrec + 1
          jbndpt = jbndpt + 1
110     continue
        jchrec = jchrec + 2
        kchrec = kchrec + 2
100   continue
      return
*
*___end of RNK___*
*
      end
      subroutine FNP   (inst,   kbnrec, kchrec, minst,  nunit)
************************************************************************
*
* FNP record indices corresponding to the beginning of the blocks con-
* taining channel information and interval of support information for
* the instrument listed in position inst in micron.dat.
*
* INPUTS:
* inst   - the instrument for which the record positions are desired
* nunit  - the unit number to which micron.dat is connected
*
* OUTPUTS:
* kchrec - the index of the first record containing information on the
*          channels for the instrument inst
* kbnrec - the index of the first record containing information on the
*          intervals of support for instrument inst
*
* LAST MODIFIED: 31 OCT 1987
************************************************************************
      integer           inst,   kbnrec, kchrec, minst,  nbpp,   ncpp
      integer           nunit
c
      call FCHBL(1, inst - 1, nunit, ncpp, kchrec)
      if (kchrec .eq. 0) kchrec = minst*3
      kchrec = kchrec + 1
      call FBNBL(1, inst - 1, nunit, nbpp, kbnrec)
      if (kbnrec .eq. 0) call FCHBL(1, minst, nunit, ncpp, kbnrec)
      kbnrec = kbnrec + 1
      return
*
*___end of FNP___*
*
      end
      subroutine FCHBL (instu,  instl,  nunit, nchbu,  nchbl)
************************************************************************
*
* FCHBL finds the first and last records in micron.kri which contain 
* information about the channels for the instruments listed between
* (and including) positions instu and instl in micron.dat.
*
* LAST MODIFIED: 31 OCT 1987
************************************************************************
      integer           ichrec, inst,   instl,  instu,  jnst,   ncha
      integer           nchbl,  nchbu,  nunit
c
      nchbu = 0
      nchbl = 0
c
c ##look for an instrument which has more than zero channels and
c ##find the first channel record of that instrument,  or exit.
c
      do 100 inst = instu, instl
        read(nunit, rec = inst) ichrec, ncha
        if (ncha .gt. 0) then
          nchbu = ichrec
          jnst  = inst
          go to 1000
        end if
100   continue
      return
c
c ##now that we found the first channel record, look at the remaining
c ##instruments for the last channel record.
c
1000  do 200 inst = jnst, instl
        read(nunit, rec = inst) ichrec, ncha
        if (ncha .gt. 0) nchbl = ichrec + 2*(ncha - 1) + 1
200   continue
      return
*
*___end of FCHBL___*
*
      end
      subroutine FBNBL (instu,  instl,  nunit,  nbnbu,  nbnbl)
************************************************************************
*
* FBNBL finds the first and last record in micron.kri containing infor-
* mation on the intervals of support of the kernel functions
* for all the instruments listed between (and including) positions instl 
* and instu in micron.dat.
*
*
* LAST MODIFIED: 31 OCT 1987
************************************************************************
      integer           ibnd,   instl,  instu,  nunit,  nbnbl,  nbnbu
      integer           nchbl,  nchbu 
c
      nbnbu = 0
      nbnbl = 0
c
c ##FIRST FIND THE FIRST AND LAST CHANNEL RECORDS AND READ
c ##THE POINTERS IN THESE RECORDS
c
      call FCHBL(instu,  instl,  nunit, nchbu, nchbl)
      if (nchbu .eq. 0) return
      read(nunit, rec = nchbu) nbnbu
      read(nunit, rec = nchbl - 1) nbnbl, ibnd
      nbnbl = nbnbl + ibnd - 1    
      return
*
*___end of FBNBL___*
*
      end
      subroutine MVB   (jnst,   knst,   nbm,    nunit)
************************************************************************
*
* MVB shifts the records containing information on the intervals of 
* support for the instruments listed between jnst and knst in micron.dat.  
* Special care is taken to insure the records are not writing over
* records which will be later be moved.  As the information in the 
* records is moved the pointers to the information are adjusted also.
*
* INPUTS:
* nbm    - the amount the records are to be moved
*
* LAST MODIFIED: 31 OCT 1987
************************************************************************
      integer           ibnd,   ibndpt, icha,   ichrec, inst,   jbnd
      integer           jnst,   lshf,   lstp,   lstr,   knst,  mrec
      integer           mshf,   mstr,   mstp,   nbm,    ncha,   nrec
      integer           nshf,   nstr,   nstp,   nunit
      double precision  areab,  areac,  dlow,   dup,    var,    xmaxb,  
     &                  xmaxc,  xmean
*-----------------------------------------------------------------------
      if (nbm .eq. 0) return
      call DLPVR(jnst, knst, nbm, lstr, lstp, lshf)
      do 100 inst = lstr, lstp, lshf
        read(nunit, rec = inst) ichrec, ncha
        if (ncha .gt. 0) then
        call DLPVR(1, ncha, nbm, nstr, nstp, nshf)
        do 110 icha = nstr, nstp, nshf
          nrec = ichrec + 2*(icha - 1)
          read(nunit, rec = nrec) ibndpt, ibnd, areac, xmaxc
          call DLPVR(1, ibnd, nbm, mstr, mstp, mshf)
          do 111 jbnd = mstr, mstp, mshf
            mrec = ibndpt + jbnd - 1
            read(nunit, rec = mrec) 
     &        dlow, dup, areab, xmaxb, xmean, var
            write(nunit, rec = mrec + nbm) 
     &        dlow, dup, areab, xmaxb, xmean, var
111       continue
          write(nunit, rec = nrec) ibndpt + nbm, ibnd, areac, xmaxc
110     continue
        end if
100   continue
      return
c
*___end of MVB___*
c
      end
      subroutine MVC   (jnst,   knst,   nbm,    nunit)
************************************************************************
c MVC shifts the records containing channel information 
c for the instruments listed between jnst and knst in micron.dat.  
c Special care is taken to insure the records are not writing over
c records which will be later be moved.  As the information in the 
c records is moved the pointers to the information are adjusted also.
c
c INPUTS:
c nbm    - the amount the records are to be moved
c nunit  - the unit number of the file being adjusted
c
c LAST MODIFIED: 31 OCT 1987
c
************************************************************************
      integer           ibndpt, ibnd,   ichrec, icha,   inst,   jnst
      integer           knst,   mshf,   mstp,   mstr,   nbm,    ncha
      integer           npar,   nrec,   nshf,   nstp,   nstr,   nunit 
      double precision  area,   var,    xmax,   xmean 
      character*8       name
c
      if (nbm .eq. 0) return
      call DLPVR(jnst, knst, nbm, mstr, mstp, mshf)
      do 100 inst = mstr, mstp, mshf
        read(nunit, rec = inst) ichrec, ncha, name, npar
        if (ncha .gt. 0) then
        call DLPVR(1, ncha, nbm, nstr, nstp, nshf)
        do 110 icha = nstr, nstp, nshf
          nrec = ichrec + 2*(icha - 1)
          read(nunit,  rec = nrec)           ibndpt, ibnd, area, xmax
          write(nunit, rec = nrec + nbm)     ibndpt, ibnd, area, xmax
          read(nunit,  rec = nrec + 1)       xmean, var
          write(nunit, rec = nrec + 1 + nbm) xmean, var
110     continue
        write(nunit, rec = inst) ichrec + nbm, ncha, name, npar
        end if
100   continue
      return
c
*___end of MCV___*
c
      end
      subroutine DLPVR (nchbu,  nchbl,  nbm,    nstr,   nstp,   nshf)
************************************************************************
c DLPVR determines the end points and direction of the do loop used to
c move a block of records.  The endpoints and direction of the do loop
c are chosen to prevent writing over records that will be moved by the
c subroutine calling DLPRV.
c
c INPUTS:
c nchbl -beginning record of block to be moved
c nchbu -ending record of block to be moved
c nbm   -number of records the block will be shifted
c
c OUTPUTS:
c nshf  -do loop counter
c nstr  -beginning do loop index
c nstp  -ending do loop index
c
c LAST MODIFIED: 31 OCT 1987
c
************************************************************************
      integer           nbm,    nchbl,  nchbu,  nshf,   nstp,   nstr
c
      if (nbm .eq. 0) return
      if (nbm .gt. 0) then
        nstr = nchbl
        nstp = nchbu
        nshf = -1
      else
        nstr = nchbu
        nstp = nchbl
        nshf = 1
      end if
      return
c
*___end of DLPVR___*
c
      end
      double precision function XKERX(diam)
************************************************************************
c XKERX is used to determine the mean diameter of the kernel function.
c
c LAST MODIFIED: 5 NOV 1987
c
************************************************************************
      double precision diam,   KERX
c
      XKERX = KERX(diam)*diam
c
      return
      end
      double precision function XXKERX(diam)
************************************************************************
c XXKERX is used to determine the variance of the kernel function.
c
c LAST MODIFIED: 5 NOV 1987
c
************************************************************************
      double precision diam,   KERX
c
      XXKERX = KERX(diam)*(diam**2)
c
      return
      end
