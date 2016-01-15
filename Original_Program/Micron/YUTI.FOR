      function UPDCHK  (fail,   iinst,  npar,   mfail,  minst,  mrecl,  
     &                  msgscr, nscr,   par,    dpar,   work,   cfiles)
      integer           fail,   iinst,  npar,   mfail,  minst,  mrecl,  
     &                  msgscr, nscr
      double precision  dpar,   par,    work
      character*(*)     cfiles
      logical           UPDCHK
      dimension         fail(*), npar(*), dpar(*), par(*), work(*), 
     &                  cfiles(*)
************************************************************************
*
* UPDCHK checks if the internal parameters are the same as the 
* parameters in micron.kri.
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           ioerr,  ipar,   joerr,  nkri
      logical           NEWVEC
*-----------------------------------------------------------------------
      if (npar(iinst) .eq. 0) then
        UPDCHK = .false.
        return
      end if
      call OPNKRI(ioerr, 2, mrecl, nkri,  cfiles)
      if (ioerr .eq. 0) then
        read(nkri, rec=iinst+minst, err=2000, iostat=joerr) 
     &  (work(ipar),ipar=1,npar(iinst))  
        UPDCHK = NEWVEC(npar(iinst), work, 1, par(iinst), minst, 
     &                  dpar(iinst), minst)
      end if
2000  if (ioerr .ne. 0) then
        call ERROR(42, nscr, msgscr, 0, 0.d+0, 'open', fail(mfail))
        UPDCHK = .false.
      else 
        close(nkri)
        if (joerr .ne. 0) 
     &    call ERROR(42, nscr, msgscr, 0, 0.d+0, 'read', fail(mfail))
      end if
      return
*
*___end of UPDCHK___*
*
      end
      subroutine OPNKRI(ioerr, istat, mrecl, nunit,  cfiles)
      character*(*)     cfiles
      dimension         cfiles(*)
************************************************************************
*
* OPNKRI is a utility used by MICRON to open the micron.kri data file.
*
* INPUTS:
* isat   = 1 if micron.kri does not already exist
*        = 2 if micron.kri does already exist
*        = 3 if micron.kri's existence is unknown
*
* OUTPUTS:
* ioerr  = the value of the iostat after attempting to open micron.kri
*
* LAST MODIFIED:30 NOV 1987
*
************************************************************************
      integer           ioerr,  istat,  lchar,  mrecl,  nunit
      character*7       chstat
      lchar = 3
*-----------------------------------------------------------------------
      if (istat .eq. 1) then
        chstat(1:lchar) = 'new'
      else if (istat .eq. 2) then
        chstat(1:lchar) = 'old'
      else 
        lchar = 7
        chstat(1:lchar) = 'unknown'
      end if
      call FUNIT(nunit)
      open(nunit, file = cfiles(4), status = chstat(1:lchar),
     &     access = 'direct', form = 'unformatted', recl = mrecl, 
     &     iostat = ioerr, err = 1000)
c
1000  return
c
*___end of OPNKRI___*
c
      end
      subroutine RCHB  (fail,   ibnd,   iinst,  jcha,   mfail, msgscr, 
     &                  nscr,   nunit,  xbnd)
      integer           fail,   ibnd,   iinst,  jcha,   mfail, msgscr, 
     &                  nscr,   nunit
      double precision  xbnd
      dimension         xbnd(2,*), fail(*)
************************************************************************
*
* RCHB reads from the micron.kri file to find out where the kernel func-
* tion is nonzero, and places this information in xbnd.
*
* Outputs: xbnd,ibnd
*
* Called by:CWK
*
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      integer           ibndpt, ioerr,  jbnd,   jchpnt
*-----------------------------------------------------------------------
      read(nunit, rec = iinst, err = 1000, iostat = ioerr) jchpnt
      jchpnt = jchpnt + 2*(jcha - 1)
      read(nunit, rec = jchpnt, err = 1000, iostat = ioerr) ibndpt, ibnd
      if ((ibnd .eq. 0) .or. (ibndpt .eq. 0))
     &  call ERROR(46, nscr, msgscr, 0, 0.d00, ' ', fail(mfail))
      do 100 jbnd = 1, ibnd
        if (fail(mfail) .gt. 0) go to 1000
        read(nunit, rec = ibndpt+jbnd-1, err = 1000, iostat = ioerr) 
     &   xbnd(1,jbnd), xbnd(2,jbnd)
        if (xbnd(1,jbnd) .gt. xbnd(2,jbnd)) 
     &    call ERROR(46, nscr, msgscr, 0, 0.d00, ' ', fail(mfail))
100   continue
c
1000  if (ioerr .ne. 0)
     &  call ERROR(42, nscr, msgscr, ioerr, 0.d00, 'read', fail(mfail))
      return
*
*___end of RCHB___*
*
      end
      subroutine UINDEX (intg)
************************************************************************
*
* UINDEX is used to set the integer fail variable to true.
*
* LAST MODIFIED:30 NOV 1987
************************************************************************
      integer          intg
c
      intg = abs(intg) + 1
c
      return
      end
      subroutine WRCHK (ident,  msglev, nuni,   msgchk, merrm,  iwrite)
************************************************************************
*
* WRCHK decides if message "ident" should be printed.
*
* INPUTS:
* ident  - message identifier
* msglev - integer parameter which controls the amount of output
*
* OUTPUTS:
* iwrite = 0 if the message should not be printed
*
* LAST MODIFIED: 6 JAN 1988
************************************************************************
      integer           ident,  iwrite, merrm,  msgchk, msglev, nuni
      dimension         msgchk(merrm)
c
      if ((ident  .le. 0) .or. (ident .gt. merrm)) then
        iwrite = 1
      else if (msglev .ge. msgchk(ident)) then
        iwrite = 1
      else
        iwrite = 0
      end if
c
      return
c
*___end of WRCHK___*
c
      end
      subroutine PARCHK(idim,   jdim,   varnm,  subnm)
************************************************************************
*
* PARCHK is used by MICRON to determine if the dimension of variable 
* declared in a common block is the same as the dimension declared in 
* the main program.
*
* LAST MODIFIED:26 JAN 1988
************************************************************************
      integer           idim,   jdim
      character*(*)     varnm,  subnm
*-----------------------------------------------------------------------
      if (idim .ne. jdim) then
        write(6, 1) varnm, subnm, jdim
        stop
      end if
c
      return
1     format(//,' PROGRAM ABORT:',/,
     &' The parameter (',a6,') in the subroutine (',a6,') needs',/,
     &' to be changed to ',i6,'.',//)
*
*___end of PARCHK___*
*
      end
      subroutine RDINST(string, keys,   ikeys,  dkeys,  subkey, maxkey, 
     &                  ntoken, tokens, maxtok, ninp,   loc1,   ierr)
      integer           ierr,   ikeys,  loc1,   maxkey, maxtok, ninp,   
     &                  ntoken
      double precision  dkeys
      character*(*)     keys,   string, subkey, tokens
      dimension         ikeys(*), dkeys(maxkey,*), keys(*), 
     &                  subkey(maxkey,*), tokens(*)
***********************************************************************
*
* RDINST searches a line of data for a key word which is contained 
* the array "keys".  The search starts at the end of the string 
* and works backwords.  After a key word is found, successive lines
* are searched for subkeys to initialize the variable associated with
* the subkey. The subkey is assumed to be the first word on the line
* of data.
*
* VARIABLES:(type, input/output, array)
*
* string    (C,I/O) - on input the character string which is searched
*                   - for a key word, and on output the last string
*                     read
* keys      (C,I,A) - the array containing the key words.
* ikeys     (I,I,A) - the array containing the number of parameters
*                   - associated with each key word.
* dkeys     (D,O,A) - the array containing the parameters associated
*                   - with each key word.
* mkeys     (I,I)   - the dimension of keys, and ikeys
* ninp      (I,I)   - the unit number of the file which contains
*                   - the parameters.
* loc1      (I,O)   - the position in the dictionary of the key word.
* ierr      (I,O)   - an index which reports on the status of the 
*                     search;
*                     ierr = 0  :no errors.
*                          = 1  :key word not found
*                          = 2  :a single piece of numeric data was not
*                                found for a subkey
* tokens    (C,-,A) - used to store words found while reading the line
*                   - of data containing the keyword.
*
* LAST MODIFIED: 7 FEB 1988
***********************************************************************
      integer           itoken, loc2,   ndigit
      double precision  dvalue
      character*16      key
      logical           sorted
      parameter        (sorted = .false.)
      dimension         dvalue(2)
*-----------------------------------------------------------------------
      itoken = ntoken 
      ierr   = 0
      loc1   = 0
c
c ##LOOK AT THE LINE OF DATA FOR A NAME
c
10    if ((loc1 .eq. 0) .and. (itoken .gt. 1)) then
        call OPLOOK(maxkey, keys, 1, sorted, tokens(itoken), loc1)
        itoken = itoken - 1
        go to 10
      end if
      itoken = itoken + 1
      loc1   = iabs (loc1)
      if (loc1 .eq. 0) then
        ierr = 1
      else
c
c ##BEGIN LOOKING FOR SUBKEYS
c
20      ntoken = maxtok
        call RDSTR (ninp,   ierr,  string,   ntoken, tokens)
        key    = tokens(1)
        call OPLOOK(ikeys(loc1), subkey(loc1,1), maxkey,sorted,key,loc2)
        if ((ntoken .ne. 0) .and. (loc2 .ne. 0)) then
          ndigit = 2
          call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue, 1)
          if (ndigit .ne. 1) then
            ierr = 2
          else
            dkeys(loc1,loc2) = dvalue(1)
            go to 20
          end if
        end if
      end if
      return
*
*___end of RDINST___*
*
      end
      subroutine KERPLT(iinst,  jcha,   msgscr, nch,    ncha,   nplot,  
     &                  nscr,   dspc,   work,   inme)
      integer           iinst,  jcha,   msgscr, nch,    ncha,   nplot,  
     &                  nscr
      double precision  dspc,   work
      character*(*)     inme
      dimension         nch(*), dspc(*), work(*), inme(*) 
************************************************************************
*
* KERPLT generates an unformatted data file containing the kernel 
* functions for the specified instrument.
* This subroutine is called when the job 'plot' is listed in the 
* micron.jbs data file. 
*
* LAST MODIFIED: 7 DEC 1988
************************************************************************
      integer           icha,   iplt,   iwrite, kcha,   nkpl
      double precision  dadd,   diam,   ediam,  KERX
      character*3       chanch
      character*8       name
      character*16      filenm
      logical           IOUTSD
*-----------------------------------------------------------------------
c
c ##FIRST DETERMINE FOR WHICH INSTRUMENT AND CHANNELS THE KERNEL
c ##FUNCTIONS NEED TO BE EVALUATED
c
      name   = inme(iinst)
      ncha   = nch(iinst)
      filenm = ' '
      call MERGE(filenm, name)
      if (IOUTSD(jcha, 1, ncha)) then
c
c ##THEN WE'LL COMPUTE THE INSTRUMENT RESPONSE FOR ALL THE CHANNELS
c
        icha = 1
        kcha = ncha
        call SCRNM(12, nscr, msgscr, 0, 0.d00, name)
      else
        icha = jcha
        kcha = jcha
        call SCRNM(14, nscr, msgscr, jcha, 0.d+0, name)
c
c ##IF WE'RE COMPUTING THE INSTRUMENT RESPONSE FOR ONLY ONE CHANNEL
c ##THEN MAKE SURE THE FILE NAME REFLECTS THIS
c
        call INTCHR(jcha, chanch)
        call MERGE(filenm, chanch)
      end if
      call OPUPPR(filenm, 'L')
      call MERGE(filenm,'.plt')
      call FUNIT(nkpl)
      open(nkpl, file = filenm, status = 'unknown')
      diam = dspc(3)
      dadd = (dspc(4) - dspc(3))/dble(nplot - 1)
      do 100 iplt = 1, nplot
        do 110 jcha = icha, kcha
          work(jcha) = KERX(diam)
110     continue
        ediam = dexp(diam)
        write(nkpl,1) ediam, (work(iwrite), iwrite = icha, kcha)
        diam = diam + dadd
100   continue
      close(nkpl)
c
1     format(100(2x,d12.6))
      return
*
*___end of KERPLT___*
*
      end
      double precision  function TSTFN(diam)
************************************************************************
c
c TSTFN computes the value of a test distribution function given a 
c discription of the function passed through the
c common block /cmtfpa/. With TF1 and TF2 we can form n-modal distribu-
c tions as well as power law models.  Note diam is actually log(diameter)
c
c INPUTS:
c msum   - the maximum number of times we can sum up TFi to form a TSTFN
c mtfpar - the maximum number of parameters we can pass to TFi
c isumfn - a flag describing which function is to be added to TSTFN
c tfpar  - parameters for TFi's
c
c LAST MODIFIED:12 NOV 1987
c
************************************************************************
      integer           isum,   isumfn, msum,   msumt,  mtfpar, mtfprt
      double precision  diam,   rmin,   tempfn, TF1,    TF2,    TF3,    
     &                  TF4,    tfpar
      logical           first
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
      common /cmtfpa/   isum,   isumfn(msum), tfpar(msum, mtfpar) 
      save   /cmtfpa/
      common /cmrmin/   rmin
      save   /cmrmin/
      common /cmfun/    msumt,  mtfprt
      save   /cmfun/
      save              first
      data first       /.true./
*-----------------------------------------------------------------------
      if (first) then
        first = .false.
        call PARCHK(msum, msumt, 'msum', 'TSTFN')
        call PARCHK(mtfpar, mtfprt, 'mtfpar', 'TSTFN')
      end if
c
c ##OBSERVE mfun = 3
c
      tempfn = rmin
      do 100 isum = 1, msum
        if      (isumfn(isum) .eq. 1) then
          tempfn = tempfn + TF1(diam)
        else if (isumfn(isum) .eq. 2) then
          tempfn = tempfn + TF2(diam)
        else if (isumfn(isum) .eq. 3) then
          tempfn = tempfn + TF3(diam)
        else if (isumfn(isum) .eq. 4) then
          tempfn = tempfn + TF4(diam)
        end if
100   continue
      TSTFN = tempfn
      return
*
*___end of TSTFN___*
*
      end
