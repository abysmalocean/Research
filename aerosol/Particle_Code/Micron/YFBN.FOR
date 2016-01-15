      subroutine BOUND (fail,   ibnd,   ichpnt, id1,    id2,    igap, 
     &                  inspnt, mdat,   mdim,   mfail,  mrecl,  msglev, 
     &                  nlog,   opchr,  dspc,   tol,    wmach,  ybnd,   
     &                  cfiles)
      integer           fail,   ibnd,   ichpnt, id1,    id2,    igap, 
     &                  inspnt, mdat,   mdim,   mfail,  mrecl,  msglev, 
     &                  nlog,   opchr
      double precision  dspc,   tol,    wmach,  
     &                  ybnd
      character*(*)     cfiles
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nlog(*), opchr(*), dspc(*), tol(*), wmach(*), 
     &                  ybnd(2,*), cfiles(*)
************************************************************************
*
* BOUND 
* 1)  finds the inversion interval based on the support of the kernel
*       functions, (ds1,ds2) 
* 2)  determines if d1 and d2 lie inside or outside the interval of
*       support,   (ds1,ds2).
*
* INPUTS:
* ibnd   - the number of subintervals contained in ybnd
* tol(9) - the tolerance used to decide if there are intervals in 
*          (d1,d2) which are not in the support of the kernel functions.
* tol(10) - the tolerance used to decide if d1 or d2 lie inside
*            (ds1,ds2)
* tol(11) - the tolerance used to decide if d1 or d2 lie outside
*            (ds1,ds2)
* OUTPUTS:
* id1   = -1 if d1 lies outside the interval of support
*       =  0 if d1 lies on the edge of the interval of support
*       = +1 if d1 lies inside the interval of support
* id2   - analagous to id1
* igap  - the number of intervals in (d1,d2) not in the support of the 
*         kernel functions.
*
* LAST MODIFIED: 12 JAN 1989
************************************************************************
      integer           ifind,  jgap,   msglog, nlog2
      double precision  d1,     d2,     diflow, difup,  dmin,   dmax,   
     &                  ds2,    ds1,    dstop,  dtemp,  eds1,   
     &                  eds2,   eps,    pctgap, xlnwdt
      character*5       lower,  upper
      logical           DOUTSD
      parameter         (lower = 'lower')
      parameter         (upper = 'upper')
*-----------------------------------------------------------------------
      if (fail(1)  .gt. 0) return
      if ((opchr(7) .eq. 0) .and. (opchr(10) .eq. 0)) return
      d1   = dspc(1)
      d2   = dspc(2)
      dmin = dspc(3)
      dmax = dspc(4)
      if (DOUTSD(d1, dmin, dmax) .or. DOUTSD(d2, dmin, dmax)) then
	call SCRNM(3, nlog(5), msglev(5), 0, 0.d00, ' ')
      else
	call SCRNM(4, nlog(5), msglev(5), 0, 0.d00, ' ')
      end if
      eps    = wmach(3)
      nlog2  = nlog(2)
      msglog = msglev(2)
c
c ##FIND ALL THE SUBINTERVALS WHICH SUPPORT THE KERNEL FUNCTIONS
c
      call RDBND (fail,   ibnd,   ichpnt, inspnt, mdat,   mfail, 
     &            mrecl,  msglev, nlog,   wmach,  ybnd,   cfiles)
      if (fail(mfail) .gt. 0) return
c
c ##"PAD" EACH SUBINTERVAL AND SORT THE PADDED SUBINTERVALS
c
      call INTVPD(ibnd, tol(8), ybnd, 2)
      call SORTBD(ibnd, ybnd)
      ds2 = ybnd(2, 1)
      ds1 = ybnd(1, 1)
      do 200 ifind = 1, ibnd
	ds2 = dmax1(ds2, ybnd(2,ifind))
200   continue
c
c ##IF THE BOUNDS CALCULATED BY MICRON ARE OUTSIDE THE MAXIMUM
c ##INVERSION INTERVAL.
c ##(DUE TO OVERPADDING OR AN OUT OF DATE micron.kri DATA FILE)
c
      if (DOUTSD(ds1, dmin, dmax)) then
        eds1 = dexp(ds1)
        call WARNG(4, nlog2, msglog, 0, eds1, lower)
        ds1  = dmin
      end if
      if (DOUTSD(ds2, dmin, dmax)) then
        eds2 = dexp(ds2)
        call WARNG(4, nlog2, msglog, 1, eds2, upper)
        ds2  = dmax
      end if
      eds1   = dexp(ds1)
      call SCRNM(5, nlog2, msglog, 0, eds1, lower)
      eds2   = dexp(ds2)
      call SCRNM(5, nlog2, msglog, 0, eds2, upper)
c
c ##IF d1 OR d2 ARE OUTSIDE THE MAXIMUM INVERSION INTERVAL 
c ##OR ARE TOO CLOSE TOGETHER, RESET THEM TO THE VALUE DETERMINED 
c ##BY MICRON.
c 
      if (DOUTSD(d1, dmin, dmax)) d1 = ds1
      if (DOUTSD(d2, dmin, dmax)) d2 = ds2
      if (d1 .gt. d2) call DSWAP(1, d1, 1, d2, 1)
      if (((d2 - d1)/dble(mdim)) .lt. eps) then
        call WARNG(12, nlog2, msglog, 0, 0.d+0, ' ')
        d1 = ds1
        d2 = ds2
      end if
c
c ##NOW CHECK IF d1 OR d2 LIE "INSIDE" (ds1,ds2)
c
      call INCHK(ibnd, id1, 1, d1, tol(10), ybnd)
      call INCHK(ibnd, id2, 2, d2, tol(10), ybnd) 
c
c ##NOW WE CHECK IF d1 OR d2 LIE "OUTSIDE" (ds1,ds2)
c
      xlnwdt = d2 - d1
      if (id1 .eq. 0) then
        diflow = (ds1 - d1)/xlnwdt
        if (diflow .gt. tol(11)) id1 = -1
      end if
      if (id2 .eq. 0) then
        difup  = (d2 - ds2)/xlnwdt
        if (difup .gt. tol(11)) id2  = -1
      end if
c
c ##FIND THE NUMBER OF LARGE SUPPORT GAPS IN THE INTERVAL (d1,d2)
c
      igap  = 0
      dtemp = dmax1(ybnd(2,1), d1)
      dstop = dmin1(d2, ds2)
      do 300 jgap = 1, ibnd - 1
	if (dtemp .ge. (ybnd(1,jgap + 1))) then
	  dtemp = dmax1(dtemp, ybnd(2,jgap + 1))
	else
	  pctgap = ((ybnd(1,jgap + 1) - dtemp)/xlnwdt)*100.d00
	  if (pctgap .gt. tol(9)) igap = igap + 1
	  dtemp = ybnd(2,jgap + 1)
	end if
300   continue
c
c ##WRITE WARNING MESSAGES IF THE INVERSION INTERVAL DOES NOT COMPARE WELL
c ##WITH THE SUPPORT DATA CONTAINED IN micron.k
c
      if (id1 .ne. 0) then
        dtemp = dexp(ds1)
        if (id1 .gt. 0) then
          call WARNG(5,  nlog2, msglog, 0, dtemp, lower)
        else 
          call WARNG(36, nlog2, msglog, 0, dtemp, lower)
        end if
      end if
      if (id2 .ne. 0) then
        dtemp = dexp(ds2)
        if (id2 .gt. 0) then
          call WARNG(5,  nlog2, msglog, 0, dtemp, upper)
        else
          call WARNG(36, nlog2, msglog, 0, dtemp, upper)
        end if
      end if
      if (igap .ne. 0)
     &  call WARNG(7, nlog2, msglog, igap, tol(9), ' ')
      dspc(1) = d1
      dspc(2) = d2
      return
c
*___end of BOUND___*
c
      end
      subroutine SORTBD(ibnd,   ybnd)
************************************************************************
c
c SORTBD sorts the intervals based on their lower diameters using a bubble
c sort. ybnd(1,1) will represent the smallest lower diameter of all
c the intervals.
c
c INPUTS:
c ibnd  - the number of intervals to be sorted
c
c LAST MODIFIED: 5 NOV 1987
c
************************************************************************
      integer           ibnd,   iech,   jech
      double precision  ybnd
      logical           exch
      dimension         ybnd(2,*)
c
      do 100 iech = ibnd, 2, -1
	exch = .false.
	do 110 jech = 1, iech - 1
	  if (ybnd(1, jech) .gt. ybnd(1, jech + 1)) then
	    exch = .true.
            call DSWAP(2, ybnd(1, jech), 1, ybnd(1, jech + 1), 1)
	  end if
110     continue
	if (.not.exch) return
100   continue
c
      return
*
*___end of SORTBD___*
*
      end
      subroutine RDBND (fail,   ibnd,   ichpnt, inspnt, mdat,   mfail,  
     &                  mrecl,  msglev, nlog,   wmach,  ybnd,   cfiles)
      integer           fail,   ibnd,   ichpnt, inspnt, mdat,   mfail, 
     &                  mrecl,  msglev, nlog
      double precision  wmach,  ybnd
      character*(*)     cfiles
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nlog(*), wmach(*), ybnd(2,*), cfiles(*)
************************************************************************
*
* RDBND reads from the data file micron.kri to find out where a kernel
* function is supported and places this information in the vector ybnd
*
* INPUTS:
* inst   - the instrument number of the kernel function of interest
* icha   - the channel number of the kernel function of interest
* ibnd   - the number of bounds contained in the vector ybnd
*
* OUTPUTS:
* ibnd   - the number of bounds contained in the vector ybnd
*
* LAST MODIFIED: 5 FEB 1988
************************************************************************
      integer           icha,   inst,   ioerr,  ipnt,   irchbd, 
     &                  irchnl, iread,  jbnd,   joerr,  nbnd,   ncha,   
     &                  ninf,   nrec   
      double precision  dlgmin
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      dlgmin = dlog(wmach(5))
      call OPNKRI(ioerr, 2, mrecl, ninf, cfiles)
      ibnd = 0
      do 100 ipnt = 1, mdat
        if (ichpnt(ipnt) .gt. 0) then
          icha = ichpnt(ipnt)
          inst = inspnt(ipnt)
          read(ninf, rec = inst, err = 1000) irchnl, ncha
          if (icha .gt. ncha) then
            call ERROR(46, nlog(5), msglev(5), inst, 0.d00, ' ', 
     &                 fail(mfail))
          else
            nrec = irchnl + 2*(icha - 1)
            read(ninf, rec = nrec, err = 1000, iostat = joerr) 
     &        irchbd, jbnd
            nbnd = ibnd + jbnd
            nrec = irchbd
            do 110 iread = ibnd + 1, ibnd + jbnd
              read(ninf, rec = nrec, err = 1000, iostat = joerr) 
     &          ybnd(1,iread), ybnd(2,iread)
              if (ybnd(1,iread) .gt. ybnd(2,iread)) then
                call ERROR(46, nlog(5), msglev(5), inst, 0.d00, ' ', 
     &                    fail(mfail))
              else if (ybnd(1,iread) .lt. dlgmin) then
                call ERROR(46, nlog(5), msglev(5), inst, 0.d00, ' ', 
     &                    fail(mfail))
              end if
              if (fail(mfail) .gt. 0) go to 1000
              nrec = nrec + 1
110         continue
            ibnd = nbnd
          end if
        end if
100   continue
c
1000  if (ioerr .ne. 0) call ERROR(42, nlog(5), msglev(5), ioerr, 0.d00, 
     &  'open', fail(mfail))
      if (joerr .ne. 0) call ERROR(42, nlog(5), msglev(5), joerr, 0.d00, 
     &  'read', fail(mfail))
      close(ninf)
      return
c
*___end of RDBND___*
c
      end
      subroutine INCHK (ibnd,   id,     iend,   dia,    tolin, 
     &                  ybnd)
      integer           ibnd,   id,     iend
      double precision  dia,    tolin,  ybnd
      dimension         ybnd(2,*)
************************************************************************
c
c INCHK is called by BOUND to determine if the endpoints of the interval
c (d1,d2) lie "inside" any of the subintervals (y(1,ibnd), y(2,ibnd))
c which make up the interval (ds1,ds2).  d1 is said to lie inside the 
c subinterval (dy1, dy2) if (d1 - dy1)/(dy2 - dy1) is greater than 
c tolin; similarly, d2 is said to lie inside the interval (dy1, dy2) if 
c (d2 - dy2) /(dy1 - dy2) is greater than tolin
c
c INPUTS:
c dia    - either d1 or d2
c iend   = 1 if dia corresponds to d1
c        = 2 id dia corresponds to d2
c ibnd   - the number of intervals contained in ybnd
c
c OUTPUTS:
c id     = 1 if dia lies inside one of the intervals in ybnd and
c          0 if not
c  
************************************************************************
      integer           ihigh,  ilow,   jbnd
      double precision  tolchk
c-----------------------------------------------------------------------
      ilow  = iend
      ihigh = 3 - ilow
      id = 0
      do 100 jbnd = 1, ibnd
        tolchk = (dia - ybnd(ilow, ibnd)) / 
     &           (ybnd(ihigh, ibnd) - ybnd(ilow, ibnd))
        if (tolchk .gt. tolin) then
          id = 1
          return
        end if
100   continue
      return
c
*___end of INCHK___*
c
      end
