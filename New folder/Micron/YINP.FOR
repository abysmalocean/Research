      subroutine INPUT (fail,   ichpnt, ierrmd, iinst,  inspnt, itnsmp, 
     &                  iwork,  kerr,   maxfev, maxtok, mdat,   mdim,   
     &                  merr,   mfail,  minfev, minsiv, minst,  mintg,  
     &                  miter,  mpar,   mrecl,  msglev, nch,    ndim,   
     &                  ninp,   ninst,  nlog,   npar,   npdat,  nsamp,  
     &                  opchr,  skip,   bl,     bspc,   bu,     decdiv, 
     &                  dpar,   dspc,   e,      par,    tol,    tpar,   
     &                  smp,    wmach,  work,   y,      cfiles, header, 
     &                  inme,   inpnm,  parnm,  tokens)
      integer           fail,   ichpnt, ierrmd, iinst,  inspnt, itnsmp, 
     &                  iwork,  kerr,   maxfev, maxtok, mdat,   mdim,   
     &                  merr,   mfail,  minfev, minsiv, minst,  mintg,  
     &                  miter,  mpar,   mrecl,  msglev, nch,    ndim,   
     &                  ninp,   ninst,  nlog,   npar,   npdat,  nsamp,  
     &                  opchr,  skip
      double precision  bl,     bspc,   bu,     decdiv, dpar,   dspc,   
     &                  e,      par,    tol,    tpar,   smp,    wmach,  
     &                  work,   y
      character*(*)     cfiles, header, inme,   inpnm,  parnm,  tokens
      dimension         fail(*), ichpnt(*), inspnt(*), iwork(*), 
     &                  msglev(*), nch(*), nlog(*), npar(*), opchr(*), 
     &                  bl(*), bspc(*), bu(*), dspc(*), dpar(*), 
     &                  e(mdat,*), par(*), tol(*), tpar(*), wmach(*), 
     &                  work(*), y(*), cfiles(*), header(*), inme(*), 
     &                  parnm(minst,*), tokens(*)
************************************************************************
*
* INP reads the data file containing the instrument readings.
*
* VARIABLES:
*
* LAST REVISION: 2 FEB 1989
************************************************************************
      integer           first,  icha,   idat,   IDSAFE, ierr,   ihead,  
     &                  jnst,   lhead,  loc1,   loc2,   loc3,   
     &                  mark,   maxdat, maxkey, maxtie, msglog, mwarn,  
     &                  nbdat,  ndigit, nlog2,  nread,  ntoken, nwarn
      double precision  dvalue, rmin,   xkinst, rmax,   zero
      character*132     line
      character*16      datkey, keys,   key,    key2,   key3,   ties
      logical           DOUTSD, newlin, newpar, NEWVEC, sorted, UPDCHK
      parameter        (maxdat = 4)
      parameter        (maxkey = 19)
      parameter        (maxtie = 11)
      parameter        (mwarn  = 5)
      parameter        (zero   = 0.d+0)
      parameter        (sorted = .true.)
      dimension         dvalue(maxdat), keys(maxkey), ties(maxtie)
      save              keys,   ties,   first
      data              keys        
     &/'CONCENTRATION   ', 'DIMENSION:  NDIM', 'ENDPOINT        ',
     & 'ENDS:   ENDPOINT', 'ERROR           ', 'HEADER          ',
     & 'INSTRUMENT      ', 'INTERVAL:  INTVL', 'INTVL: INVERSION', 
     & 'INVERSION       ', 'ITERATIONS      ', 'ITNS: ITERATIONS', 
     & 'NDIM:   SOLUTION', 'SAMPLE          ', 'SKIP            ', 
     & 'SMOOTHING       ', 'SMP             ', 'SOLUTION        ', 
     & 'TITLE:    HEADER'/
      data ties
     &/'BOTH            ', 'DENSE           ', 'DIAGONAL        ', 
     & 'ITERATIONS      ', 'ITNS: ITERATIONS', 'LOWER           ', 
     & 'NEITHER:    NONE', 'NONE            ', 'PARAMETER       ',
     & 'SPARSE          ', 'UPPER           '/
      data              first /.true./
*-----------------------------------------------------------------------
      if ((fail(1) .gt. 0) .or. (skip .gt. 0)) return
      call SCRNM(22, nlog(5), msglev(5), 0, 0.d+0, inpnm)
c
c ##INITIALIZE
c
      rmin      = wmach(5)
      rmax      = wmach(7)
      nlog2     = nlog(2)
      msglog    = msglev(2)
      npdat     = 0
      idat      = 0
      icha      = 0
      ninst     = 0
      opchr(7)  = 0
      opchr(10) = 0
      bspc(3)   = -rmax
      bspc(4)   = -rmax
      if (first) then
        first   = .false.
        bspc(1) = -rmax
        bspc(2) = -rmax
      end if
      call DLOAD(mdat*merr, zero, e, 1)
      ntoken  = maxtok
      call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
c
c ##BEGIN LOOKING FOR KEYWORDS.  AUTOMATICALLY LOOK FOR NUMERIC DATA.
c
10    if ((ntoken .ge. 1) .and. (fail(4) .le. 0)) then
        newlin = .true.
        ndigit = 2
        call TOKNUM(ndigit, ntoken, tokens(1), dvalue, 1)
        key    = tokens(1)
        key3   = tokens(ntoken)
        call OPLOOK(maxkey, keys, 1, sorted, key, loc1)
        if (key .eq. 'SMOOTHING') then
          key2 = tokens(2)
          call OPLOOK(maxtie, ties, 1, sorted, key2, loc2)
          if (key2 .eq. 'ITERATIONS') then
            key = 'ITERATIONS'
          else
            key = 'SMP'
          end if
        end if
        if (key .eq. 'INSTRUMENT') then
          newlin = .false.
          call RDINST(line,   inme,   npar,   par,   parnm,  minst,  
     &                ntoken, tokens, maxtok, ninp,   jnst,  ierr)
c
c ##CHECK IF AN ERROR OCCURRED WHILE READING THE INSTRUMENT INFORMATION
c
          if (ierr .ne. 0) then
            call ERROR(14, nlog2, msglog, ierr, 0.d+0, line, fail(4))
c
c ##MAKE SURE THERE IS NOT TOO MUCH DATA
c
          else if (idat+nch(jnst) .gt. mdat) then
            call ERROR(13, nlog2, msglog, mdat, 0.d+0, inpnm, fail(4))
c
c ##MAKE SURE THERE ARE NOT TOO MANY INSTRUMENTS
c
          else if (ninst+1 .gt. minsiv) then
            call ERROR(9, nlog2, msglog, minsiv, 0.d+0, inpnm, fail(4))
c
c ##ELSE WE HAD A SUCCESSFUL READ.  CHECK IF KERNEL FUNCTION 
c ##VALUES SHOULD BE RECOMPUTED
c
          else
            ninst = ninst + 1
            if (opchr(10) .le. 0) then
              if (iabs(inspnt(idat+1)) .ne. jnst) then
                opchr(10) = 1
              else if (npar(jnst) .gt. 0) then
                call DCOPY(npar(jnst), dpar(jnst), minst, work, 1)
                call DSCAL(npar(jnst), 0.1d+0, work, 1)
                if (NEWVEC(npar(jnst), par(jnst), minst, 
     &            tpar(ninst), minsiv, work, 1)) opchr(10) = 1
              end if
            end if
c
c ##CHECK IF micron.kri NEEDS TO BE UPDATED
c
            if ((opchr(10) .eq. 1) .and. (npar(jnst) .gt. 0)) then
              newpar = .true.
              if (iabs(inspnt(idat+1)) .eq. jnst) then
                if (iabs(ichpnt(idat+1)) .eq. 1) then
                  call DCOPY(npar(jnst), dpar(jnst), minst, work, 1)
                  call DSCAL(npar(jnst), 0.1d+0, work, 1)
                  newpar =  NEWVEC(npar(jnst), par(jnst), minst, 
     &                      tpar(ninst), minsiv, work, 1) 
                end if
              end if
              if (newpar) then
                iinst = jnst
                if (UPDCHK(fail,   jnst,   npar,   mfail,  minst,  
     &                     mrecl,  msglev(5), nlog(5),   par,    dpar,   
     &                     work,   cfiles))
     &          call UPDATE(fail,   maxfev, mfail,  minfev, minst,  
     &                      mintg,  mpar,   mrecl,  msglev(5), nch,    
     &                      npar,   nlog(5),   decdiv, dspc,   par,    
     &                      tol,    wmach,  work,   cfiles, inme,   
     &                      parnm)
              end if
            end if
            call DCOPY(npar(jnst), par(jnst), minst, tpar(ninst),minsiv)
c
c       ##READ THE DATA FOR THE INSTRUMENT IN POSITION iinst IN 
c       ##THE DATA FILE.  
c
            icha = 0
20          if ((icha .lt. nch(jnst)) .and. (fail(4) .le. 0)) then
              icha = icha + 1
              idat = idat + 1
              call OPLOOK(3, datkey, 1, sorted, tokens(1), loc1)
              if (ntoken .eq. 0) then
                call ERROR(16,nlog2,msglog,0,0.d+0,inme(jnst),fail(4))
              else if (tokens(1) .eq. 'SKIP') then
                dvalue(1) = -1.d+0
                if (inspnt(idat) .gt. 0) opchr(7) = 1
	        inspnt(idat) = -jnst
                ichpnt(idat) = -icha
              else
                ndigit = 4
                call TOKNUM(ndigit, ntoken, tokens(1), dvalue, 1)
                if ((ndigit .lt. 2) .and. (dvalue(1) .ge. -rmin)) then
                  xkinst = dble(ninst)
                  call ERROR(19, nlog2, msglog, icha, xkinst,inme(jnst),
     &                       fail(4))
                else if (dvalue(1) .lt. -rmin) then
                  if (inspnt(idat) .gt. 0) opchr(7) = 1
	          inspnt(idat) = -jnst
                  ichpnt(idat) = -icha
                else
                  if ((inspnt(idat) .ne. jnst) .or. 
     &                (ichpnt(idat) .ne. icha)) opchr(7) = 1
	          inspnt(idat)   = jnst
                  ichpnt(idat)   = icha
                  npdat          = npdat + 1
                  y(npdat)       = dvalue(1)
                  e(npdat,npdat) = dvalue(2)
                end if
              end if
c
c         ##READ DATA BOUNDS (IF ANY)
c
              nbdat = mdim + npdat
              if (ndigit .eq. 3) then
                if (dvalue(3) .le. dvalue(1)) then
                  bl(nbdat) = dvalue(3)
                  bu(nbdat) = -rmax
                else
                  bl(nbdat) = -rmax
                  bu(nbdat) = dvalue(3)
                end if
              else if (ndigit .eq. 4) then
                bl(nbdat) = dvalue(3)
                bu(nbdat) = dvalue(4)
              else
                bl(nbdat) = -rmax
                bu(nbdat) = -rmax
              end if
c
c         ##SET NEGATIVE VALUES TO NEGATIVE INFINITY AND MAKE SURE 
c         ##IF BOUNDS EXIST THEN bl <= bu.
c
              if (bl(nbdat) .le. -rmin) 
     &          bl(nbdat) = -rmax
              if (bu(nbdat) .le. -rmin) then
                bu(nbdat) = -rmax
              else if (bl(nbdat) .gt. bu(nbdat)) then
                call DSWAP(1, bl(nbdat), 1, bu(nbdat), 1)
              end if
              ntoken = maxtok
              call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
              go to 20
            end if
          end if
c
c ##CONCENTRATION CONSTRAINTS;  IF NO NUMBERS WERE FOUND THEN
c ##ASSUME NO CONSTRAINTS
c
        else if (key .eq. 'SAMPLE') then
          if (ndigit .ge. 1) nsamp = IDSAFE(dvalue(1))
        else if (key .eq. 'CONCENTRATION') then
          if (ndigit .eq. 2) then
            bspc(3) = dvalue(1)
            bspc(4) = dvalue(2)
            if (bspc(3) .le. -rmin) bspc(3) = -rmax
            if (bspc(4) .le.  rmin) then
              bspc(4) = -rmax
            else if (bspc(3) .gt. bspc(4)) then
              call DSWAP(1, bspc(3), 1, bspc(4), 1)
            end if
          else if ((ndigit .eq. 1) .and. (dvalue(1) .gt. rmin)) then
            bspc(3) = dvalue(1)
            bspc(4) = dvalue(1)
          else
            bspc(3) = -rmax
            bspc(4) = -rmax
          end if
c
c ##SOLUTION DIMENSION
c
        else if (key .eq. 'SOLUTION') then
          ndim = IDSAFE(dvalue(1))
          if (ndim .gt. mdim) then
            ndim = mdim
            call WARNG(1, nlog2, msglog, mdim, 0.d+0, ' ')
          end if
c
c ##SMOOTHING PARAMETER
c
        else if (key .eq. 'SMP') then
          smp = -dabs(dvalue(1))
c
c ##ENDPOINT CONSTRAINTS
c
        else if (key .eq. 'ENDPOINT') then
          bspc(1) = -rmax
          bspc(2) = -rmax
          call OPLOOK(maxtie, ties, 1, sorted, key3, loc3)
          if (key3 .eq. 'UPPER') then
            bspc(2) = zero
          else if (key3 .eq. 'LOWER') then
            bspc(1) = zero
          else if (key3 .eq. 'BOTH') then
            bspc(1) = zero
            bspc(2) = zero
          end if
        else if (key .eq. 'HEADER') then
          if (ihead .lt. 2) ihead = ihead + 1
          key2 = tokens(2)
          call OPLOOK(1, 'ERASE           ', 1, sorted, key2, loc2)
          if (key2 .eq. 'ERASE') then
            header(ihead) = ' '
          else if (ntoken .ge. 2) then
            first = 1
            nread = 2
            lhead = len(header(1))
            call FNWORD(line, nread, first,  mark)
            header(ihead) = line(first:lhead)
          end if
        else if (key .eq. 'ITERATIONS') then
c
c ##IF THE NUMBER OF SMOOTHING ITERATIONS SPECIFIED IS A NUMBER
c ##THEN SET itnsmp EQUAL TO THAT NUMBER, ELSE ASSUME THE MAXIMUM
c ##IS DESIRED.
c
          if (ndigit .ge. 1) then
            itnsmp = IDSAFE(dvalue(1))
            if (itnsmp .gt. miter) then
              call WARNG(16, nlog2, msglog, miter, 0.d+0, inpnm)
              itnsmp = miter
            else if (itnsmp .lt. 0) then
              call WARNG(17, nlog2, msglog, 0, 0.d+0, inpnm)
              itnsmp = 0
            end if
          else
            itnsmp = miter
          end if
c
c ##INTERVAL OF INVERSION
c
        else if (key .eq. 'INVERSION') then
          opchr(10) = 1
          dspc(1)  = dvalue(1)
          dspc(2)  = dvalue(2)
          if (dspc(1) .le. rmin) dspc(1) = -rmax
          if (dspc(1) .gt. rmin) dspc(1) = dlog(dspc(1))
          if (dspc(2) .le. rmin) dspc(2) = -rmax
          if (dspc(2) .gt. rmin) then
            dspc(2) = dlog(dspc(2))
            if (dspc(2) .lt. dspc(1)) 
     &        call DSWAP(1, dspc(1), 1, dspc(2), 1)
          end if
          if (DOUTSD(dspc(1), dspc(3), dspc(4))) dspc(1) = -rmax
          if (DOUTSD(dspc(2), dspc(3), dspc(4))) dspc(2) = -rmax
c
c ##HAVE MICRON PERFORM LEAST SQUARES
c
        else if (key .eq. 'SKIP') then
          opchr(8) = 0
c
c ##ERROR MODEL
c
        else if (key .eq. 'ERROR') then
          call OPLOOK(maxtie, ties, 1, sorted, key3, loc3)
          if (key3 .eq. 'DIAGONAL') then
            ierrmd = 1
          else if (key3 .eq. 'SPARSE') then
            newlin = .false.
            ierrmd = 2
          else if (key3 .eq. 'DENSE') then
            newlin = .false.
            ierrmd = 3
          else
            call WARNG(10, nlog2, msglog, 0, 0.d+0, key3)
            ierrmd = 1
          end if
          call RDERR (fail,   ierrmd, kerr,   maxtok, mdat,   merr,   
     &                msglog, ninp,   nlog2,  npdat,  work,   e,      
     &                tokens)
        else
c
c ##KEYWORD IS NOT RECOGNIZED
c
          call WARNG(8, nlog2, msglog, 0, 0.d+0, key)
          nwarn = nwarn + 1
        end if
        if (nwarn .ge. mwarn) then
          call ERROR(1, nlog2, msglog, 0, 0.d+0, inpnm, fail(4))
        else if (newlin) then
          ntoken = maxtok
          call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
        end if
        go to 10
      end if
c
c ##CLEAN THE POINTER ARRAYS AND CHECK THE NUMBER OF POSITIVE DATA
c
      if (idat .lt. mdat) then
        call ILOAD(mdat-idat, 0, ichpnt(idat+1), 1)
        call ILOAD(mdat-idat, 0, inspnt(idat+1), 1)
      end if
      if (npdat .le. 0) then
        call ERROR(7, nlog2, msglog, 0, 0.d+0, inpnm, fail(4))
      else if (npdat .le. 16) then
        if (opchr(1) .lt. 0) opchr(1) = -3
      else
        if (opchr(1) .lt. 0) opchr(1) = -1
      end if
      if (fail(4) .gt. 0) call UINDEX(fail(1))
      close(ninp)
      return
*
*___end of INPUT___*
*
      end
      subroutine RDERR (fail,   ierrmd, kerr,   maxtok, mdat,   merr,   
     &                  msglog, ninp,   nlog2,  npdat,  dvalue, e,      
     &                  tokens)
      integer           fail,   ierrmd, kerr,   maxtok, mdat,   merr,   
     &                  msglog, ninp,   nlog2,  npdat
      double precision  dvalue, e
      character*(*)     tokens
      dimension         fail(*), dvalue(*), e(mdat,*), tokens(*)
************************************************************************
*
* RDERR reads in the array (if any) which describes the error model of 
* the data.  For now the only error model we allow is normal.  e(i,j)
* represents the coefficient of a normal random variable which has zero 
* mean and unit standard deviation.
*
* INPUTS:
* ierrmd = 1 :the error model is diagonal and has already been read.
*
*        = 2 :the error model is sparse;  MICRON looks for the i,j'th
*             indices of the (normal) error model and reads in
*             e(i,j) 
*
*        = 3 :the error model is dense;  read in a full array.
*             For the users benefit the diagonal elements must
*             be included and will be reread. The subroutine tries
*             to read (e(1,1), e(1,2), ... e(1,kerr),
*             e(2,1)... ... e(npdat,kerr)
*
* OUTPUTS:
* kerr   - the number of independent error sources
* e      - the error model
*
* LAST REVISION: 3 FEB 1989
************************************************************************
      integer           icol,   IDSAFE, ierr,   ispars, jspars, loc1,   
     &                  maxkey, ndigit, nrow,   ntoken
      double precision  dmerr
      character*132     line
      character*16      key,    keys
      logical           IOUTSD, sorted
      parameter        (maxkey = 1)
      parameter        (sorted = .true.)
      save              keys
      data keys        /'ROW             '/
*-----------------------------------------------------------------------
      if (fail(5) .gt. 0) return
      dmerr  = dble(merr)
      kerr   = npdat
c
c ##THE ERROR MODEL IS SPARSE.  EVERY ENTRY MUST HAVE ONLY 3 NUMERIC
c ##PIECES OF DATA.
c
      if (ierrmd .eq. 2) then
10      ntoken = maxtok
        call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
        ndigit = ntoken
        call TOKNUM(ndigit, ntoken, tokens(1), dvalue, 1)
        if ((ndigit .eq. ntoken) .and. (ndigit .eq. 3)) then
          ispars = IDSAFE(dvalue(1))
          jspars = IDSAFE(dvalue(2))
c
c ##CHECK TO MAKE SURE THE INDICES ARE REASONABLE
c
          if (IOUTSD(ispars, 1, mdat)) then
            call WARNG(11, nlog2, msglog, ispars, dvalue(3),  'first')
          else if (IOUTSD(jspars, 1, merr)) then
            call WARNG(11, nlog2, msglog, jspars, dvalue(3), 'second')
          else
            e(ispars, jspars) = dvalue(3)
            if (kerr .lt. jspars) kerr = jspars
          end if
          go to 10
        end if
c
c ##ELSE THE ERROR MODEL IS DENSE. 
c
      else if (ierrmd .eq. 3) then
        ntoken = maxtok
        call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
        ndigit = ntoken
        call TOKNUM(ndigit, ntoken, tokens(1), dvalue, 1)
20      key    = tokens(1)
        call OPLOOK(maxkey, keys, 1, sorted, key, loc1)
c
c ##READ THE ROW NUMBER IF THE KEYWORD IS 'ROW'
c
        if ((fail(5) .le. 0) .and. (key .eq. 'ROW')) then
          if (ndigit .ne. 1) then
            call ERROR(21, nlog2, msglog, 0, 0.d+0, line, fail(4))
          else
            nrow = IDSAFE(dvalue(1))
          end if
          if (IOUTSD(nrow, 1, mdat)) 
     &      call ERROR(23, nlog2, msglog, mdat, 0.d+0, ' ', fail(4))
c
c ##BEGIN READING THE ENTRIES
c
          icol = 0
30        if (fail(5) .le. 0) then
            ntoken = maxtok 
            call RDSTR (ninp,   ierr,  line,   ntoken, tokens)
            ndigit = ntoken
            call TOKNUM(ndigit, ntoken, tokens(1), dvalue, 1)
            if ((ndigit .gt. 0) .and. (ndigit .eq. ntoken)) then
              if (icol+ndigit .gt. merr) then
                call ERROR(29, nlog2, msglog, merr, 0.d+0, ' ', fail(4))
              else
                call DCOPY(ndigit, dvalue, 1, e(nrow,icol+1), mdat)
                icol = icol + ndigit
              end if
              if (icol .gt. kerr) kerr = icol
              go to 30
            end if
          end if
          go to 20
        end if
      end if
      return
*
*___end of RDERR___*
*
      end
      subroutine ECHO  (fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,   
     &                  mdat,   mdim,   minsiv, minst,  miter,  msglev, 
     &                  ndim,   ninst,  nlog,   npar,   npdat,  opchr,  
     &                  bl,     bspc,   bu,     dspc,   e,      par,    
     &                  smp,    tpar,   wmach,  y,      header, inme,   
     &                  parnm,  smooth)
      integer           fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,   
     &                  mdat,   mdim,   minsiv, minst,  miter,  msglev, 
     &                  ndim,   ninst,  nlog,   npar,   npdat,  opchr
      double precision  bl,     bspc,   bu,     dspc,   e,      par,    
     &                  smp,    tpar,   wmach,  y
      character*(*)     header, inme,   parnm,  smooth
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nlog(*), npar(*), opchr(*), bl(*), bspc(*), 
     &                  bu(*), dspc(*), e(mdat,*), par(*), tpar(*), 
     &                  wmach(*), y(*), header(*), inme(*), parnm(*), 
     &                  smooth(*)
************************************************************************
*
* ECHO will echo all the variables that have reasonable values if
* a reading error has occurred, or an echo is requested. The output is 
* placed in the root.ech file.
*
*
* LAST REVISION: 2 FEB 1989
************************************************************************
      integer           idat,   inst,   jdat,   jerr,   jinst,  kinst,  
     &                  nech,   nparid
      double precision  ed1,    ed2,    rmin
      character*16      errmod
      character*74      datlin
      logical           failed, IOUTSD, lecho
      dimension         errmod(3)
      save              errmod
      data              errmod
     &                 /'DIAGONAL        ', 'SPARSE          ', 
     &                  'DENSE           '/
*-----------------------------------------------------------------------
      if ((fail(2) .gt. 0) .or. (fail(3) .gt. 0)) return
      rmin   = wmach(5)
      nech   = nlog(1)
      failed = (fail(4) .gt. 0) .or. (fail(5) .gt. 0)
      lecho  = (msglev(1) .ge. 30)  .or. 
     &        ((msglev(1) .ge. 10) .and. (failed))
      if (.not. lecho) then
        if ((msglev(1)  .gt. 0) .and. (nech .gt. 0)) close(nech)
        return
      end if
      call SCRNM(23, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
      write(nech, 12)
      write(nech, '(1x,a75)') header(1)
      write(nech, '(1x,a75)') header(2)
      write(nech, 1)  ninst, npdat
      if (.not. IOUTSD(ndim, 5, mdim)) write(nech, 24) ndim
      if (.not. IOUTSD(itnsmp, 0, miter)) write(nech, 25) itnsmp
      if (dspc(1) .ge. dspc(3)) then
        ed1 = dexp(dspc(1))
        write(nech, 23) 'Lower', ed1
      end if
      if (dspc(2) .gt. dspc(3)) then
        ed2 = dexp(dspc(2))
        write(nech, 23) 'Upper', ed2
      end if
      write(nech,2) smooth(iabs(opchr(1))), smp
      write(nech,3)
      inst  = 0
      kinst = 0
      idat  = 0
      do 100 idat = 1, mdat
        if (inspnt(idat) .gt. 0) then
          if (inspnt(idat) .ne. kinst) then
c
c ##then we're echoing data for a new instrument
c
            inst   = inst + 1
            kinst  = inspnt(idat)
            nparid = npar(kinst)
            call DCOPY(npar(kinst), tpar(inst), minsiv, 
     &                 par(kinst), minst)
            write(nech,5)
            write(nech,6) inme(kinst)
            if (nparid .gt. 0) 
     &        write(nech,4) (parnm(jinst), par(jinst), jinst = kinst, 
     &        (nparid - 1)*minst + kinst, minst)
            write(nech,11)
          end if
c
          datlin = ' '
          write(datlin(3:6), '(i3)') ichpnt(idat)
          write(datlin(10:19), '(e10.3)') y(idat)
          write(datlin(31:40), '(e10.3)') e(idat,idat)
          if (bl(mdim + idat) .ge. -rmin) 
     &      write(datlin(50:59), '(e10.3)') bl(mdim+idat)
          if (bl(mdim + idat) .ge. -rmin) 
     &      write(datlin(64:73), '(e10.3)') bl(mdim+idat)
          write(nech, '(1x,a74)') datlin
        end if
100   continue
c
c ##ECHO THE CONCENTRATION CONSTRAINTS AND ENDPOINT CONSTRAINTS
c
      if (bspc(1) .ge. -rmin) then
        write(nech, 18) bspc(1), 'lower'
      else
        write(nech, 19) 'lower'
      end if
      if (bspc(2) .ge. -rmin) then
        write(nech, 18) bspc(2), 'upper'
      else
        write(nech, 19) 'upper'
      end if
      if (bspc(3) .ge. -rmin) then
        write(nech, 21) 'lower', bspc(3)
      else 
        write(nech, 20) 'Lower'
      end if
      if (bspc(4) .ge. -rmin) then
        write(nech, 21) 'upper', bspc(4)
      else
        write(nech, 20) 'Upper'
      end if
c
      write(nech,9) errmod(ierrmd)
      if (ierrmd .gt. 1) then
        do 200 jdat = 1, npdat
          write(nech,10) jdat
          if (ierrmd .eq. 2) then
            call NZPRN(kerr, rmin, e(jdat,1), mdat, nech)
          else if (ierrmd .eq. 3) then
            write(nech,11) (e(jdat,jerr), jerr = 1, kerr)
          end if
200     continue
      end if
      close(nech)
c
1     format(//
     &' # of instruments',t35,' = ',i4,/,
     &' # of positive data',t35,' = ',i4)
2     format(' smoothing parameter (',a4,')',t35,' = ',d12.5)
3     format(//,' Instruments and data:')
4     format(' parameters: ',100(t15,2(a8,'=',d11.3,5x),/))
5     format(///)
6     format(' instrument:   ',a5)
7     format(1x,i3,7x,e10.3,2x,e10.3,8x,a4,8x,a4)
8     format(1x,i3,7x,e10.3,2x,e10.3,2x,e10.3,2x,e10.3)
9     format(///' Error array:',a16,//)
10    format(//' Row ',i4,' of error array. . .')
11    format(/,' Channel', t10, 'Input data', t23, 
     &'Diagonal error component', t50,
     &'lower bound', t64, 'upper bound')
12    format(//,' This is data file echos a MICRON input data file.',//)
18    format(//' Size distribution constrained to equal ',d10.3,/,
     &         ' at the ',a5,' diameter of the inversion interval.')
19    format(/,' The size distribution at the ',a5,' end of the',/,
     &         ' inversion interval is unconstrained.')
20    format(/,1x,a5,' concentration bound: none')
21    format(/,1x,a5,' concentration bound: ',d10.3, ' #/cm^3')
23    format(/,1x,a5,' interval diameter',t35,' = ',d10.3,' microns')
24    format(/,' dimension of solution',t35,' = ',i4)
25    format(/,' smoothing iterations ',t35,' = ',i4)
26    format(6(2x,d11.4))
      return
*
*___end of ECHO___*
*
      end
