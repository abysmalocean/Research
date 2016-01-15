      subroutine MICDAT(fail,   maxtok, mdat,   mfail,  minst,  mpar,   
     &                  mrecl,  msglev, msmth,  msuffx, mword,  nch,    
     &                  nlog,   npar,   nterm,  opchr,  dpar,   dspc,   
     &                  rmax,   rmin,   par,    tol,    wmach,  cfiles, 
     &                  fsuffx, inme,   parnm,  smooth, tokens)
      integer           fail,   maxtok, mdat,   mfail,  minst,  mpar,   
     &                  mrecl,  msglev, msmth,  msuffx, mword,  nch,    
     &                  nlog,   npar,   nterm,  opchr
      double precision  dpar,   dspc,   par,    rmax,   rmin,   tol,    
     &                  wmach
      character*(*)     cfiles, fsuffx, inme,   parnm,  smooth, tokens
      dimension         fail(*), nch(*), msglev(*), nlog(*), npar(*),
     &                  opchr(*), dpar(minst,*), dspc(*), par(minst,*), 
     &                  tol(*), wmach(*), cfiles(*), inme(*), 
     &                  parnm(minst,*), tokens(*)
************************************************************************
*
* MICDAT reads the data file micron.dat to learn about the instruments,
* tolerances, and options that will be used in the calculations.
* micron.dat is created by the user prior to data inversion.
*
* LAST REVISION: 1 FEB 1989
************************************************************************
      integer           IDSAFE, ierr,   iinst,  ioerr,  IOUTSD, ipar,   
     &                  itoken, loc1,   loc2,   locmsg, loctie, maxdat, 
     &                  maxkey, maxtie, mpart,  ndat,   ndigit, nmic,   
     &                  nread,  ntoken
      double precision  dvalue, zero
      character*132     line
      character*16      key,    keys,   tie,    ties
      logical           cmpmch, newlin, sorted
c
c ##maxdat is the maximum of ndat
c
      parameter        (maxdat = 4)
      parameter        (maxkey = 15)
      parameter        (maxtie = 2)
      parameter        (zero   = 0.d+0)
      parameter        (sorted = .true.)
      dimension         ndat(maxkey), dvalue(maxdat), keys(maxkey),
     &                  ties(maxtie)
      data              cmpmch /.true./
      data              ndat
     & /0,                  0,                  2,                  
     &  2,                  4,                  0,                  
     &  1,                  1,                  1,                  
     &  0,                  0,                  0,                  
     &  1,                  0,                  0 /
      data              keys
     &/'CONSTRAIN       ', 'INPSPC          ', 'INSTRUMENT      ', 
     & 'LARGEST         ', 'MACHINE         ', 'MASS            ', 
     & 'MAXIMUM         ', 'MESSAGE         ', 'NORM            ', 
     & 'NUMBER          ', 'OUTSPC          ', 'SMOOTHING       ', 
     & 'TOLERANCE       ', 'VOLUME:     MASS', 'WEAK            '/
      data              ties
     &/'NO              ', 'YES             '/
*-----------------------------------------------------------------------
      mpart = mpar * minst
      call ILOAD(mfail,   0, fail,   1)
      call ILOAD(minst,   0, nch,    1)
      call ILOAD(minst,   0, npar,   1)
      call CHLOAD(minst, ' ', inme,  1)
      call CHLOAD(mpart, ' ', parnm, 1)
      call FUNIT(nmic)
      open(nmic, err = 2000, file = cfiles(1), status = 'old',
     &     iostat = ioerr)
2000  if (ioerr .ne. 0) then
        write(nlog(5), 2) cfiles(1), ioerr
        call UINDEX(fail(mfail))
        return
      end if
      rewind(nmic)
      ntoken = maxtok
      call RDSTR (nmic,   ierr,  line,   ntoken, tokens)
c
c ##BEGIN LOOKING FOR KEYWORDS
c
20    if ((ntoken .ge. 1) .and. (fail(mfail) .le. 0)) then
        newlin = .true.
        key    = tokens(1)
        call OPLOOK(maxkey, keys, 1, sorted, key, loc1)
        loc1   = abs(loc1)
        if (loc1 .ne. 0) then
          nread = ndat(loc1)
          call TOKNUM(nread, ntoken-1, tokens(2), dvalue, 1)
          if (nread .eq. ndat(loc1)) then
            if (key .eq. 'LARGEST') then
              dspc(3) = dvalue(1)
              dspc(4) = dvalue(2)
            else if (key .eq. 'INSTRUMENT') then
              iinst = IDSAFE(dvalue(1))
              if (IOUTSD(iinst, 1, minst)) then
                call ERROR(41, nlog(5), msglev(5), minst, 0.d+0, ' ', 
     &                     fail(mfail))
              else
                inme(iinst) = tokens(3)
                nch(iinst)  = IDSAFE(dvalue(2))
                if (IOUTSD(nch(iinst), 1, mdat)) call ERROR(51, nlog(5), 
     &            msglev(5), mdat, 0.d+0, tokens(3), fail(mfail))
              end if
c
c ##START READING PARAMETERS
c
              ipar = 0
30            if ((ipar .lt. mpar) .and. (fail(mfail) .le. 0)) then
                ntoken = maxtok
                call RDSTR (nmic,   ierr,  line,   ntoken, tokens)
                key = tokens(1)
                call OPLOOK(maxkey, keys,  1, sorted, key,    loc2)
                if ((ntoken .gt. 0) .and. (loc2 .eq. 0)) then
                  ipar = ipar + 1
                  parnm(iinst,ipar) = tokens(1)
                  ndigit  = 2
                  call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
                  if (ndigit .eq. 0) then
                    call ERROR(53, nlog(5), msglev(5), ipar, 0.d+0, 
     &                         inme(iinst), fail(mfail))
                  else if (ndigit .eq. 1) then
                    par(iinst,ipar)  = dvalue(1)
                    dpar(iinst,ipar) = dvalue(1) * 0.02d+0
                  else
                    par(iinst,ipar)  = dvalue(1)
                    dpar(iinst,ipar) = dvalue(2)
                  end if
                  go to 30
                else
                  newlin = .false.
                end if
              end if
              npar(iinst) = ipar
            else if (key .eq. 'CONSTRAIN') then
              opchr(2) = 1
            else if (key .eq. 'MACHINE') then
              cmpmch    = .false.
              if (dvalue(1) .gt. zero) wmach(3)  = dvalue(1)
              if (dvalue(2) .gt. zero) wmach(5)  = dvalue(2)
              if (dvalue(3) .gt. zero) wmach(7)  = dvalue(3)
              if (dvalue(4) .gt. zero) wmach(11) = dvalue(4)
              wmach(4)  = dsqrt(wmach(3))
              wmach(6)  = dsqrt(wmach(5))
              wmach(8)  = dsqrt(wmach(7))
            else if (key .eq. 'MAXIMUM') then
              nterm = IDSAFE(dvalue(1))
            else if (key .eq. 'MESSAGE') then
              do 100 itoken = 2, ntoken
                if (tokens(itoken)(1:1) .eq. '.') then
                  tie = tokens(itoken)
                  call OPLOOK(msuffx, fsuffx,  1, sorted, tie, locmsg)
                  if (locmsg .gt. 0) msglev(locmsg) = IDSAFE(dvalue(1))
                end if
100           continue
            else if (key .eq. 'TOLERANCE') then
              if (dvalue(1) .gt. zero) then
                tol(2) = dvalue(1)
                opchr(3) = 0
              else
                opchr(3) = 1
              end if
            else if (key .eq. 'INPSPC') then
              call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie, 
     &                    sorted)
              opchr(14) = 0
              if (tie .eq. 'YES') opchr(14) = 1
            else if (key .eq. 'OUTSPC') then
              call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie, 
     &                    sorted)
              opchr(15) = 0
              if (tie .eq. 'YES') opchr(15) = 1
            else if (key .eq. 'MASS') then
              opchr(11) = 3
            else if (key .eq. 'NORM') then
              opchr(12) = nint(dabs(dvalue(1)))
            else if (key .eq. 'NUMBER') then
              opchr(11) = 1
            else if (key .eq. 'WEAK') then
              call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie, 
     &                    sorted)
              opchr(13) = 1
              if (tie .eq. 'NO') opchr(13) = 0
            else if (key .eq. 'SMOOTHING') then
c
c ##READ THE TYPE OF SMOOTHING OPTIMUM (THE LAST WORD IN THE LINE)
c
              call FNDKEY(smooth, msmth, tokens, ntoken, tie,    loctie, 
     &                    sorted)
              if (loctie .le. 0) then
                call ERROR(74, nlog(5), msglev(5), ntoken, 0.d+0, tie, 
     &                     fail(mfail))
              else
                opchr(1) = loctie
              end if
            end if
          else
c
c ##ELSE THE CORRECT NUMBER OF NUMBERS WERE NOT FOUND
c
            call ERROR(62, nlog(5), msglev(5), nread, 0.d+0, key, 
     &                 fail(mfail))
          end if
        else
c
c ##ELSE THE KEYWORD WAS NOT RECOGNIZED
c 
          call ERROR(63, nlog(5), msglev(5), 0, 0.d+0, key, fail(mfail)) 
        end if 
c
c ##READ THE NEXT LINE OF DATA
c
        if (newlin) then
          ntoken = maxtok
          call RDSTR (nmic,   ierr,  line,   ntoken, tokens)
        end if
        go to 20
      end if
      if ((fail(mfail) .le. 0) .and. (cmpmch)) 
     &  call MCHPAR()
c
c ##MICRON works with ln(diam)
c
      if (fail(mfail) .le. 0) then
        dspc(3) = dlog(dspc(3))
        dspc(4) = dlog(dspc(4))
        rmin    = wmach(5)
        rmax    = wmach(7)
        nlog(5) = IDSAFE(wmach(11))
        dspc(1) = -rmax
        dspc(2) = -rmax
      end if
      close (nmic)
c
      return
2     format(//
     &' An error occurred while trying to open the ',a12,/,
     &' data file.  The iostat was set equal to ', i5,'.')
*
*___end of MICDAT___*
*
      end
      subroutine MICCHK(fail,   mfail,  minst,  mpar,   mrecl,  mword, 
     &                  nlog,   nterm,  dspc,   wmach,  parnm)
      integer           fail,   mfail,  minst,  mpar,   mrecl,  mword, 
     &                  nlog,   nterm
      double precision  dspc,   wmach
      character*(*)     parnm
      dimension         fail(*), nlog(*), dspc(*), wmach(*), parnm(*)
************************************************************************
*
* MICCHK checks to make sure the variables read in MICDAT are reasonable
* e.g. positive diameters, tolerances. 
*
* While checking the inputs we feel the following limits are realistic
*
* 0.d+0      .lt. wmach(5)   .lt. wmach(3)   .lt.   1.d-06
* 1.d+10     .lt. wmach(7) 
* wmach(3)   .lt.   tol(3)   .lt.   1.d-01
* wmach(3)   .lt.   tol(4)   .lt.   1.d-01
* wmach(3)   .lt.   tol(5)   .lt.   1.d-01
* wmach(3)   .lt.   tol(6)   .lt.   1.d-01
*
* These limits should not be chosen arbitrarily, but depend on the
* the data.  Future versions of MICRON should compute some of these.
*
* LAST REVISION: 11 NOV 1987
************************************************************************
      integer           mrt,    mwrite, nout
      double precision  tdmax,  wm5min, wm7min
*-----------------------------------------------------------------------
c ##IF AN ERROR OCCURS IN MICCHK IT IS AUTOMATICALLY REPORTED TO UNIT 6 
c
      parameter        (mwrite = 100)
c
c ##AND HERE ARE THE PRESET LIMITS ON THE TOLERANCES
c
      parameter        (wm5min = 0.d00)
      parameter        (wm7min = 1.d+10)
      parameter        (tdmax  = 1.d+06)
c
      if (fail(mfail) .gt. 0) then
        nout = 6
      else
        nout = nlog(5)
      end if
c
c ##FIRST MAKE SURE THE RECORDS IN MICRON.KRI ARE LONG ENOUGH
c 
      if (len(parnm(1)) .gt. mword) then
        mrt = len(parnm(1)) * mpar
      else
        mrt = mword         * mpar
      end if
      if (mrt .ge. mrecl) 
     &  call ERROR(52, nout, mwrite, mrt, 0.d00, ' ', fail(mfail))
      wmach(5) = dabs(wmach(5))
      if (nterm .lt. 1) nterm = 1
c
c ##FIRST CHECK THE MACHINE PARAMETERS
c
      if (fail(mfail) .le. 0) then
        if (wmach(5) .lt. wm5min)
     &    call ERROR(66, nout, mwrite, 0, wmach(5), ' ', fail(mfail))
        if (wmach(5) .gt. wmach(3))
     &    call ERROR(35, nout, mwrite, 0, wmach(5), ' ', fail(mfail))
        if (wmach(7) .lt. wm7min) 
     &    call ERROR(57, nout, mwrite, 0, wmach(7), ' ', fail(mfail))
        if (dspc(3) .gt. dspc(4)) call DSWAP(1, dspc(3), 1, dspc(4), 1)
        if (dexp(dspc(3)) .le. wmach(5))
     &    call ERROR(38, nout, mwrite, 0, dspc(3), ' ', fail(mfail))
        if (dexp(dspc(4)) .gt. tdmax)
     &    call ERROR(65, nout, mwrite, 0, dspc(4), ' ', fail(mfail))
      end if
      return
*
*___end of MICCHK___*
*      
      end
      subroutine MICECH(fail,   mfail,  minst,  mpar,   msglev, nch,    
     &                  nlog,   npar,   nterm,  opchr,  dspc,   tol,    
     &                  wmach,  cfiles, inme,   parnm,  smooth)
      integer           fail,   mfail,  minst,  mpar,   msglev, nch,    
     &                  nlog,   npar,   nterm,  opchr
      double precision  dspc,   tol,    wmach
      character*(*)     cfiles, inme,   parnm,  smooth
      dimension         fail(*), msglev(*), nch(*), nlog(*), npar(*), 
     &                  opchr(*), dspc(*), tol(*), wmach(*), cfiles(*), 
     &                  inme(*), parnm(*),  smooth(*)
************************************************************************
*
* MICECH is called whenever a failure occurs while reading micron.dat or
* if msgech [msg(1)] is greater than or equal to 20. The output is 
* written to micron.dmp.
*
* LAST REVISION:15 DEC 1987
************************************************************************
      integer           iinst,  ioerr,  imsg,   ipar,   nech
      double precision  tdmax,  tdmin
*-----------------------------------------------------------------------
      ioerr = 0
      if ((fail(mfail) .ne. 0) .or. (msglev(1) .ge. 20)) then
        if (fail(mfail) .ne. 0) then
          tdmax = dspc(4)
          tdmin = dspc(3)
        else
          tdmax = dexp(dspc(4))
          tdmin = dexp(dspc(3))
        end if
        call FUNIT(nech)
        open(nech, err = 2000, iostat = ioerr, file = cfiles(2),
     &       status = 'unknown')
        rewind(nech)
        write(nech, 3) cfiles(1)
c
        do 100 iinst = 1, minst
          if (nch(iinst) .gt. 0)
     &    write(nech, 4) iinst, inme(iinst), nch(iinst), npar(iinst)
100     continue
        write(nech, 2)
        do 200 iinst = 1, minst
          if (npar(iinst) .gt. 0) then
            write(nech, 5) inme(iinst), (parnm(ipar), ipar = iinst, 
     &                     (npar(iinst) - 1)*minst + iinst, minst)
            if (mod(npar(iinst), 5) .eq. 0) backspace(nech)
          end if
200     continue
c
        write(nech, 11) wmach(3), wmach(5), wmach(7), nlog(5)
        write(nech, 7 ) tdmin,    tdmax
        write(nech, 13) nterm
        write(nech, 14) (msglev(imsg), imsg = 1, 5)
        if (opchr(11) .eq. 3) write(nech, 15) 
        write(nech, 16) smooth(iabs(opchr(1)))
        close(nech)
      end if
c
2000  if (ioerr .ne. 0) write(nlog(5), 9)
      return
c
2     format(//,' List of parameter names:',/)
3     format(' This file echos the data read from ',a12,//,
     &' Insrument information:',//, ' position',
     &t13, 'Instrument name', 
     &t31, 'Number of channels', 
     &t52, 'Number of parameters',/)
4     format(t3,i3,t17,a8,t38,i3,t58,i3)
5     format(1x,a8,' parameter names:',100(t28,5(a8,2x),/))
6     format(2x,'tol(',i1,') = ',d12.5)
7     format(//,' Lower diameter limit = ',d10.3,/,
     &          ' Upper diameter limit = ',d10.3,//)
9     format(' Unable to open file micron.dmp.')
10    format(//,' Failure occurred while reading micron.dat',/,
     &' Check micron.dmp for a data echo.')
11    format(//,' Machine precision          = ',d10.3,/,
     &          ' Smallest positive number   = ',d10.3,/,
     &          ' Largest  positive number   = ',d10.3,/,
     &          ' Default output unit number = ',i4,//)
13    format(' Number of successive failures before termination:',i6)
14    format(//,' Message levels:',/,' msgech = ',i6,/,' msglog = ',
     &i6,/,' msgout = ',i6,/,' msgplt = ',i6,/,' msgscr = ',i6,//)
15    format(' The distribution is rescaled.')
16    format(' Smoothing parameter defined by ',a6,'.',//)
*
*___end of MICECH___*
*
      end
      subroutine LBLRD (mlabel, msglev, nlog,   cfiles, label)
      integer           mlabel, msglev, nlog
      character*(*)     cfiles, label
      dimension         msglev(*), nlog(*), cfiles(*), label(*)
************************************************************************
* 
* LBLRD reads micron.lbl to reset the output labels used by MICRON.  
* The labels should be changed whenever the inversion is not the 
* determination of an aerosol number distribution versus particle
* size (microns).
* The labels have the following meaning
*
* label(1)  :
*
* LAST MODIFIED: 6 SEPT 1988
************************************************************************
      integer           i,      ioerr,  nlbl
      logical           ex
*-----------------------------------------------------------------------
      inquire(file = cfiles(6), exist = ex)
      if (.not. ex) return
      call FUNIT(nlbl)
      open(nlbl, file = cfiles(6), status = 'unknown', err = 1000) 
      read(nlbl,*,end = 2000, err = 2000, iostat = ioerr) 
     &    (label(i), i = 1, mlabel)
c
2000  if (ioerr .eq. 0) 
     &   call SCRNM(6, nlog(5), msglev(5), 0, 0.d+0, ' ')
      close(nlbl)
1000  return
*
*___end of LBLRD___*
*
      end
