      subroutine SIMDAT(fail,   ichpnt, inspnt, isumfn, maxfev, maxtok, 
     &                  mdat,   mfail,  mfun,   minfev, minsiv, minst,  
     &                  mintg,  mpar,   mrecl,  msgscr, msum,   nch,    
     &                  ninst,  npar,   npdat,  nscr,   nseed,  ntfpar, 
     &                  decdiv, dpar,   dspc,   e,      errmin, errpct, 
     &                  par,    tfpar,  tol,    tpar,   ttfpar, wmach,  
     &                  work,   y,      cfiles, inme,   parnm,  smfile, 
     &                  tdstnm, tfnme,  tfpnm,  tokens)
      integer           fail,   ichpnt, inspnt, isumfn, maxfev, maxtok, 
     &                  mdat,   mfail,  mfun,   minfev, minsiv, minst,  
     &                  mintg,  mpar,   mrecl,  msgscr, msum,   nch,    
     &                  ninst,  npar,   npdat,  nscr,   nseed,  ntfpar
      double precision  decdiv, dpar,   dspc,   e,      errmin, errpct, 
     &                  par,    tfpar,  tol,    tpar,   ttfpar, wmach,  
     &                  work,   y
      character*(*)     cfiles, inme,   parnm,  smfile, tdstnm, tfnme,  
     &                  tfpnm,  tokens
      dimension         fail(*), ichpnt(*), inspnt(*), isumfn(*), 
     &                  nch(*), npar(*), ntfpar(*), dpar(*), dspc(*), 
     &                  e(mdat,*), errmin(*), errpct(*), par(minst,*), 
     &                  tol(*), tfpar(msum,*), tfpnm(minst,*), 
     &                  tpar(minsiv,*), ttfpar(mfun,*), wmach(*), 
     &                  work(*), y(*), cfiles(*), inme(*), 
     &                  parnm(minst,*), tfnme(*), tokens(*)
************************************************************************
*
* SIMDAT will simulate noisy data from a set of instruments and 
* create a data file which can be read by MICRON and inverted.
*
* INPUTS:
* errmin - is used as the baseline error in determining the standard 
*          deviation to which is used to make the .inp file
* errpct: same as errmin except errpct is error which is a fraction
*         of the data.
*
* LAST MODIFIED: 6 FEB 1988
*                2 July 1996 by R.J.Weber
************************************************************************
      integer           ibnd,   icha,   idat,   iierr,  iinst,  ioerr,  
     &                  ipar,   iplot,  ipnt,   jbnd,   joerr,  kinst,  
     &                  ncha,   nidg,   ninf,   ninp,   nplot,  nplt
      double precision  delta,  dplot,  eplot,  eps100, rmin,   SIMPRD, 
     &                  TSTFN,  ytemp,  xintge, rwtemp
      character*8       name
      character*12      ofname, filin, tstplt
      logical           adderr, new,    UPDCHK
      common /cmkerf/   iinst,  icha,   ncha
      save   /cmkerf/
      save              adderr
      external          SIMPRD
      data              adderr /.false./
*-----------------------------------------------------------------------
      if (fail(mfail) .gt. 0) return
      rmin   = wmach(5)
      eps100 = wmach(3) * 100.d+0
      call FUNIT(nidg)
      call MAKEFN(smfile, '.sim', filin)
      open(nidg, file = filin, form = 'formatted', status = 'old', 
     &  err = 1000, iostat = ioerr)
      rewind(nidg)
      call SCRNM(21, nscr, msgscr, 0, 0.d+0, filin)
c
c ##FIRST WE NEED TO FIND OUT WHAT DATA IS TO BE SIMULATED
c
      call SIMRD (fail,   ichpnt, inspnt, isumfn, maxtok, mdat,   
     &            mfail,  mfun,   minsiv, minst,  msgscr, msum,   
     &            nch,    nidg,   ninst,  npar,   npdat,  nplot,  
     &            nscr,   nseed,  ntfpar, errmin, errpct, par,    
     &            tfpar,  tpar,   ttfpar, filin,  inme,   parnm,  
     &            tdstnm, tfnme,  tfpnm,  tokens, adderr)
      close(nidg)
      if (fail(mfail) .gt. 0) return
c
c ##MAKE A DATA FILE OF TEST FUNCTION VALUES IF NECESSARY
c
      if (nplot .gt. 0) then
        call MAKEFN(tdstnm, '.plt', tstplt)
        call OPUPPR(tstplt, 'L')
        call FUNIT(nplt)
        open(nplt, file = tstplt, status = 'unknown')
        delta = (dspc(4) - dspc(3))/dble(nplot - 1)
        dplot = dspc(3)
        do 100 iplot = 1, nplot
          eplot = dexp(dplot)
          write(nplt,8) eplot, TSTFN(dplot)
          dplot = dplot + delta
100     continue
        close(nplt)
      end if
c
c ##AND THEN OPEN THE .inp DATA FILE
c
      call OPNKRI(joerr, 2, mrecl, ninf, cfiles)
      if (joerr .ne. 0) go to 1000 
      call MAKEFN(smfile, '.inp', ofname)
      call FUNIT(ninp)
      open(ninp, file = ofname, status = 'unknown')
      write(ninp, *) ' #This data has been simulated by MICRON.'
      write(ninp, *) ' #The test distribution is ', tdstnm,'.'
c
c ##MAKE THE NOISE VECTOR
c
      if (adderr) 
     &  call DNLOAD(npdat, nseed, work, 1)
c
c ##THIS IS THE LOOP WHICH ACTUALLY GENERATES THE DATA
c
      kinst = 0
      do 200 idat = 1, npdat
        iinst = inspnt(idat)
        icha  = ichpnt(idat)
        name  = inme(iinst)
        if (iinst .eq. 0) go to 1000
c
c ##DETERMINE IF WE ARE SIMULATING DATA FOR A NEW INSTRUMENT
c
        if (idat .eq. 1) then
          new = .true.
        else if (iinst .ne. inspnt(idat - 1)) then
          new = .true.
        else
          new = .false.
        end if
        if (new) then
c
c ##IF WE ARE SIMULATING A DATA FOR A NEW INSTRUMENT, WRITE IT'S NAME 
c ##AND PARAMETERS 
c
          kinst = kinst + 1
          call DCOPY(npar(iinst), tpar(kinst,1), minsiv, 
     &               par(iinst,1), minst)
          write(ninp, 6)
          write(ninp,4) name
          do 210 ipar = 1, npar(iinst)
            write(ninp, 3) parnm(iinst,ipar), par(iinst,ipar)
210       continue
          call SCRNM(13, nscr, msgscr, 0, 0.d00, name)
          ncha = nch(iinst)
c
c ##CHECK IF WE NEED TO UPDATE micron.kri
c
          if (UPDCHK(fail,   iinst,  npar,   mfail,  minst,  mrecl,  
     &               msgscr, nscr,   par,    dpar,   work(mdat+1),   
     &               cfiles))
     &      call UPDATE(fail,   maxfev, mfail,  minfev, minst,  mintg, 
     &                  mpar,   mrecl,  msgscr, nch,    npar,   nscr,   
     &                  decdiv, dspc,   par,    tol,    wmach,  
     &                  work(mdat+1),   cfiles, inme,   parnm)
        end if
        y(idat) = rmin
c
c ##READ micron.kri TO FIND OUT WHERE THE KERNEL FUNCTION IS NONZERO
c ##AND "PAD" THE INTERVALS IF NECESSARY
c
        call RCHB  (fail,   ibnd,   iinst,  icha,   mfail, 
     &              msgscr, nscr,   ninf,   work(mdat+1))
        if (fail(mfail) .gt. 0) then
          go to 1000
        else if (ibnd .gt. 0) then
          call INTVPD(ibnd,  tol(8), work(mdat+1), 2)
        end if
        xintge = -1.d+0
        ipnt   = mdat + 1
        do 220 jbnd = 1, ibnd
c
c ##  OVER EACH INTERVAL OF SUPPORT, INTEGRATE THE KERNEL FUNCTION WITH
c ##  TEST DISTRIBUTION
c
          call DINTG(SIMPRD, iierr, mintg, ytemp, work(ipnt), 
     &               work(ipnt+1), tol(6), rwtemp, rmin)
          if ((iierr .ne. 0) .and. (xintge .lt. 0.d+0)) then
            call WARNG(19, nscr, msgscr, icha, xintge, name)
            xintge = dble(iierr)
          end if
          y(idat) = y(idat) + ytemp
          ipnt = ipnt + 2
220     continue
        e(idat, idat) = errmin(kinst) + 
     &                  errpct(kinst)*y(idat)
        if (adderr) y(idat) = y(idat) + e(idat,idat) * work(idat)
        write(ninp, 2) y(idat), e(idat, idat) 
200   continue
c 
c ##SET OR COMMENT ON THE VALUES OF THE INVERSION OPTIONS
c
      write(ninp, 7) tdstnm
c
c ##CHECK FOR READING ERRORS AND CLOSE FILES
c
1000  if (ioerr .ne. 0) then
        call ERROR(58, nscr, msgscr, 0, 0.d00, filin, fail(mfail))
      else if (joerr .ne. 0) then
        call ERROR(42, nscr, msgscr, joerr, 0.d+0, 'open', fail(mfail))
      else
        close(ninp)
        close(ninf)
      end if
c
c ##PREVENT INPUT FROM ASSUMING THE INVERSION IS A REPEAT
c
      inspnt(1) = 0
      ichpnt(1) = 0
      return
1     format(2x,e11.4,3x,e11.4,3x,i5,3x,'''',a8,'''')
2     format(100(3x,d11.4))
3     format(3x, a8, ' = ', e12.4)
4     format(/,' INSTRUMENT = ',a8)
5     format(/,i2)
6     format(//)
7     format(/,
     &' HEADER:This data has been simulated by MICRON.',/,
     &' HEADER:The test distribution is ',a8,'.',/,
     &' # skip least squares calculation ',/,
     &' # concentration constraints = ',/,
     &' # inversion interval        = ',/,
     &' # solution dimension        = ',/,
     &' # endpoint constraints      = ',/,
     &' # smoothing parameter       = ',/,
     &' # smoothing iterations      = ',/,
     &' # error model               = diagonal')
8     format(2x,d12.4,2x,d12.4)
*
*___end of SIMDAT___*
*
      end
      subroutine SIMRD (fail,   ichpnt, inspnt, isumfn, maxtok, mdat,   
     &                  mfail,  mfun,   minsiv, minst,  msgscr, msum,   
     &                  nch,    nidg,   ninst,  npar,   npdat,  nplot,  
     &                  nscr,   nseed,  ntfpar, errmin, errpct, par,    
     &                  tfpar,  tpar,   ttfpar, filin,  inme,   parnm,  
     &                  tdstnm, tfnme,  tfpnm,  tokens, adderr)
      integer           fail,   ichpnt, inspnt, isumfn, maxtok, mdat,   
     &                  mfail,  mfun,   minsiv, minst,  msgscr, msum,   
     &                  nch,    ndigit, nidg,   ninst,  npar,   npdat,  
     &                  nplot,  nscr,   nseed,  ntfpar
      double precision  emin,   epct,   errmin, errpct, par,    tfpar,  
     &                  tpar,   ttfpar 
      character*(*)     filin,  inme,   parnm,  tdstnm, tfnme,  tfpnm,  
     &                  tokens
      logical           adderr
      dimension         fail(*), ichpnt(*), inspnt(*), isumfn(*), 
     &                  nch(*), npar(*), ntfpar(*), dvalue(2), 
     &                  errmin(*), errpct(*), par(minst,*), 
     &                  tfpar(msum,*), tfpnm(minst,*), tpar(minsiv,*), 
     &                  ttfpar(mfun,*), inme(*), parnm(minst,*), 
     &                  tfnme(*), tokens(*)
************************************************************************
*
* SIMRD is called by SIMDAT to find out about the test distribution that 
* will be used to generate test data, and to find out which instruments 
* to use when generating test data. 
*
* The instruments used to simulate data are listed at the beginning of
* the data file.  The description of a single instrument is as folllows:
* 1) INSTRUMENT       = 'instrument name'
*    'parameter name' = 'parameter value'
*    ERROR            = 'error parameters'
* 
* After the last instrument description the user needs to describe
* the test distribution which will be used to simulate data.  
*
* The test distribution name is assumed to follow the last instrument
* description.  The test distribution name must not be an abbreviated
* instrument name or MICRON will assume it is an instrument name. 
* If a data file containing test function values is desired, then
* an integer must be on the same line as and following the test 
* distribution name. 
*
* Next, the user needs to describe each test function.  The description
* of a single test function is as follows:
* 1) the test function name is the first test function name recognized 
*    on the next line of data.  
* 2) the test function parameters are assumed to start with the first
*    number following the test function name.  The test function 
*    parameters may be preceeded by one string of words on a single
*    line.
*
* VARIABLES:(type, input/output, array)
*
* LAST MODIFIED: 5 JAN 1989
************************************************************************
      integer           IDSAFE, ierr,   ipnt,   isum,   itfun,  jnst,   
     &                  loc1,   loctie, maxkey, maxtie, ntoken
      double precision  dvalue, zero
      character*132     line
      character*16      key,    keys,   tie,    ties
      logical           newlin, sorted
      parameter        (maxkey = 8)
      parameter        (maxtie = 2)
      parameter        (zero   = 0.d+0)
      parameter        (sorted = .true.)
      dimension         keys(maxkey), ties(maxtie)
      save              emin,   epct,   keys
      data              keys
     &/'ADD_ERROR       ', 'FUNCTION:   TEST', 'INSTRUMENT      ', 
     & 'NAME            ', 'NUMBER          ', 'POINTS:   NUMBER', 
     & 'SEED            ', 'TEST            '/
      data              ties
     &/'NO              ', 'YES             '/
      data              emin,   epct
     &                 /0.d+0, 0.025d+0/
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      call ILOAD(msum, 0, isumfn, 1)
      call ILOAD(mdat, 0, ichpnt, 1)
      call ILOAD(mdat, 0, inspnt, 1)
      line      = ' '
      tokens(1) = ' '
      npdat     = 0
      ninst     = 0
      isum      = 0
      ntoken    = maxtok
      call RDSTR (nidg,   ierr,  line,   ntoken, tokens)
10    if ((fail(mfail) .le. 0) .and. (ntoken .ne. 0)) then
        newlin = .true.
        key    = tokens(1)
        call OPLOOK(maxkey, keys,  1, sorted, key,    loc1)
        if      (key .eq. 'INSTRUMENT') then
          newlin = .false.
c
c ##READ THE INSTRUMENT NAME AND PARAMETERS 
c
          call RDINST(line,   inme,   npar,   par,   parnm,  minst,  
     &                ntoken, tokens, maxtok, nidg,   jnst,  ierr)
c
c ##CHECK IF AN ERROR OCCURRED WHILE READING THE INSTRUMENT INFORMATION
c
          if (ierr .gt. 0) then
            call ERROR(2, nscr, msgscr, ierr, 0.d+0, line, fail(mfail))
c
c ##MAKE SURE THERE IS NOT TOO MUCH DATA
c
          else if (npdat+nch(jnst) .gt. mdat) then
            call ERROR(4, nscr, msgscr, mdat, 0.d+0, filin, fail(mfail))
c
c ##MAKE SURE THERE ARE NOT TOO MANY INSTRUMENTS
c
          else if (ninst+1 .gt. minsiv) then
            call ERROR(5, nscr, msgscr, minsiv, 0.d+0, filin, 
     &                 fail(mfail))
c
c ##ELSE WE HAD A SUCCESSFUL READ
c
          else 
            ninst = ninst + 1
            call DCOPY(npar(jnst), par(jnst,1), minst, 
     &                 tpar(ninst,1), minsiv)
            do 110 ipnt = npdat + 1, npdat + nch(jnst)
              ichpnt(ipnt) = ipnt - npdat
              inspnt(ipnt) = jnst
110         continue
            npdat  = npdat + nch(jnst)
            errpct(ninst) = epct
            errmin(ninst) = emin
          end if
        else if (key .eq. 'TEST') then
          newlin = .false.
          call RDINST(line,   tfnme,  ntfpar, ttfpar, tfpnm,  mfun,
     &                ntoken, tokens, maxtok, nidg,   itfun,  ierr)
c
c ##CHECK IF AN ERROR OCCURRED WHILE READING THE TEST FUNCTION 
c ##INFORMATION
c
          if (ierr .gt. 0) then
            call ERROR(6, nscr, msgscr, ierr, 0.d+0, line, fail(mfail))
          else if (isum+1 .gt. msum) then
            call ERROR(10,nscr, msgscr, msum, 0.d+0, filin, fail(mfail))
          else 
            isum         = isum  + 1
            isumfn(isum) = itfun
c
c ##ttpar IS A TEMPORARY HOLDER OF THE TEST FUNCTION PARAMETERS
c
            call DCOPY(ntfpar(itfun), ttfpar(itfun,1), mfun, 
     &                 tfpar(isum,1), msum)
          end if
        else if (key .eq. 'SEED') then
          ndigit = 1
          call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
          nseed = -abs(IDSAFE(dvalue(1)))
        else if (key .eq. 'NUMBER') then
          ndigit = 1
          call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
          nplot = IDSAFE(dvalue(1))
        else if (key .eq. 'NAME') then
          tdstnm = tokens(ntoken)
        else if (key .eq. 'ADD_ERROR') then
          adderr = .true.
          call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie, 
     &                sorted)
          if (tie .eq. 'NO') adderr = .false.
        else if (key .eq. 'ERROR') then
c
c ##RESET THE TWO ERROR PARAMETERS
c
          ndigit = 2
          call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
          emin   = dvalue(1)
          epct   = dvalue(2)
        end if
        if (newlin) then
          ntoken = maxtok
          call RDSTR (nidg,   ierr,  line,   ntoken, tokens)
        end if
        go to 10
      end if
      if (fail(mfail) .le. 0) then
        if (ninst .le. 0) call WARNG(26, nscr, msgscr, 0, 0.d+0, filin)
        if (isum .eq. 0)
     &    call ERROR(79, nscr, msgscr, 0, 0.d+0, filin, fail(mfail))
      end if
c
      return
      end
*
*___end of SIMRD___*
*
      subroutine SIMINI(fail,   maxtok, mfail,   mfun,   msgscr, mtfpar, 
     &                  nscr,   ntfpar, ttfpar,  cfiles, tfnme,  tfpnm,  
     &                  tokens)
      integer           fail,   maxtok, mfail,   mfun,   msgscr, mtfpar, 
     &                  nscr,   ntfpar
      double precision  ttfpar
      character*(*)     cfiles, tfnme,  tfpnm,   tokens
      dimension         fail(*), ntfpar(*), ttfpar(mfun,*), cfiles(*), 
     &                  tfnme(*), tfpnm(mfun,*), tokens(*)
************************************************************************
*
* SIMINI reads the micron.sim data file to find out the name, parameter
* names, and default parameter values of each test function.
*
* VARIABLES:
*
* LAST MODIFIED:25 JAN 1989
************************************************************************
      integer           IDSAFE, ierr,   ifun,   IOUTSD, itfp,   joerr,  
     &                  loc1,   loc2,   maxkey, ndigit, nfun,   nini,   
     &                  nread,  ntoken
      double precision  dvalue
      character*132     line
      character*16      key,    keys
      logical           lfirst, sorted
      parameter        (maxkey = 2)
      parameter        (sorted = .true.)
      dimension         dvalue(2), keys(2)
      save              lfirst
      data lfirst       /.true./
      data              keys
     &/'FUNCTION   :TEST', 'TEST            '/
*-----------------------------------------------------------------------
c
c ##IF THE micron.sim DATA FILE HAS BEEN READ THEN RETURN
c
      if (lfirst) then
        lfirst = .false.
      else
        return
      end if
      line = ' '
      call ILOAD(mfun, 0, ntfpar, 1)
      call CHLOAD(mfun, ' ', tfnme, 1)
      call CHLOAD(mfun*mtfpar, ' ', tfpnm, 1)
      call FUNIT(nini)
      open(nini, file = cfiles(7), status = 'old', err = 1000,
     &     iostat = joerr)
      ntoken = maxtok
      ierr   = 0
      nfun   = 0
c
c ##BEGIN READING THE DATA FILE
c
      call RDSTR (nini,   ierr,  line,   ntoken, tokens)
10    if ((fail(mfail) .le. 0) .and. (ntoken .gt. 0)) then
        key = tokens(1)
        call OPLOOK(maxkey, keys,  1, sorted, key,    loc1)
        if (key .eq. 'TEST') then
          nfun  = nfun + 1
          nread = 1
          call TOKNUM(nread, ntoken-1, tokens(2), dvalue, 1)
          ifun = IDSAFE(dvalue(1))
          if (IOUTSD(ifun,1,mfun)) then
            call ERROR(11, nscr, msgscr, mfun, 0.d+0, ' ', fail(mfail))
c
c ##THE THIRD WORD IS THE TEST FUNCTION NAME
c
          else
            tfnme(ifun) = tokens(3)
c
c ##READ THE PARAMETER NAMES AND DEFAULT VALUES
c ##THE DESIRED FORMAT FOR THE PARAMETER NAMES IS 'NAME = VALUE ' 
c ##AN ERROR IS ASSUMED IF (1) THE NUBMER OF VALUES FOUND IS NOT 1
c
            itfp = 0
20          if ((itfp .le. mtfpar) .and. (fail(mfail) .le. 0)) then
              ntoken = maxtok
              call RDSTR (nini,   ierr,  line,   ntoken, tokens)
              key   = tokens(1)
              call OPLOOK(maxkey, keys,  1, sorted, key,    loc2)
              if ((ntoken .gt. 0) .and. (key .ne. 'TEST')) then
                itfp = itfp + 1
                tfpnm(ifun,itfp) = key
                ndigit = 2
                call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
                ttfpar(ifun,itfp) = dvalue(1)
                if (ndigit .ne. 1) call ERROR(12, nscr, msgscr, ndigit, 
     &                                        0.d+0, key, fail(mfail))
                go to 20
              end if
              ntfpar(ifun) = itfp
            end if
          end if
        else
c
c ##ELSE KEYWORD IS NOT RECOGNIZED
c
          call ERROR(15, nscr, msgscr, 0, 0.d+0, key, fail(mfail))
        end if
        go to 10
      end if
      if ((nfun .eq. 0) .and. (fail(mfail) .le. 0)) 
     &  call ERROR(60, nscr, msgscr, 0, 0.d+0, ' ', fail(mfail))
c
c
1000  if (joerr .ne. 0) then
        call ERROR(72, nscr, msgscr, joerr, 0.d+0, ' ', fail(mfail))
      else
        close(nini)
      end if
      return
*
*___end of SIMINI___*
*
      end 
      double precision function SIMPRD(diam)
************************************************************************
*
* SIMPRD is the product of the test distribution function and the kernel
* function;  it is integrated to simulate noiseless data from an 
* instrument.  
*
* LAST MODIFIED:12 JAN 1989
************************************************************************
      double precision diam,   KERX,   TSTFN
c
      SIMPRD = KERX(diam)*TSTFN(diam)
c
      return
      end
      double precision  function TSTMOM(diam)
************************************************************************
*
* TSTMOM caclulates TSTFN(diam) * diam**ipower
*
* LAST MODIFIED:11 JAN 1988
*
************************************************************************
      integer           ipower
      double precision  diam,   TSTFN
      common /cmipow/   ipower
      save   /cmipow/
c
      TSTMOM = TSTFN(diam) * diam**ipower
c
      return
      end
