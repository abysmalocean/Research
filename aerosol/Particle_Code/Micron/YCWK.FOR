      subroutine CWK   (fail,   ichpnt, inspnt, mdat,   mdim,   mfail,  
     &                  minsiv, minst,  mintg,  mrecl,  msglev, nch,    
     &                  ndim,   nlog,   npar,   npdat,  opchr,  dspc,   
     &                  diam,   fk,     par,    tol,    tpar,   wmach,  
     &                  work,   cfiles, inme)
      integer           fail,   ichpnt, inspnt, mdat,   mdim,   mfail,  
     &                  minsiv, minst,  mintg,  mrecl,  msglev, nch,    
     &                  ndim,   nlog,   npar,   npdat,  opchr
      double precision  dspc,   diam,   fk,     par,    tol,    tpar,   
     &                  wmach,  work
      character*(*)     cfiles, inme
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nch(*), nlog(*), npar(*), diam(*), dspc(*), 
     &                  fk(mdat,*), par(minst,*), tpar(minsiv,*), 
     &                  tol(*), work(*), wmach(*), inme(*), opchr(*)
************************************************************************
*
* CWK computes the matrix "fk" which satisfies [fk][f]=Y where "f" is 
* the size distribution and Y is the data vector.
*
* INPUTS:
* inspnt  - a vector whose i'th element points to the instrument which
*           generated the data "y"(i)
* ichpnt  - a vector whose i'th element points to the channel which
*           generated the data "y"(i)
* mintg   - the maxixmum number of function evaluations allowed
*           in each call to the integration routine
* npdat   - the number of data (or the number of rows of "fk")
* ndim    - the dimension of the solution vector "f" (or the number
*           of columns of "fk")
*
* OUTPUTS:
* fk     - the matrix which represents the set of kernel functions that
*          correspond to the data Y
*
* LAST MODIFIED: 6 JAN 1988
************************************************************************
      integer           ibnd,   idat,   idia,   idim,   idimp1, ierr,   
     &                  ierx,   iinst,  ipnt,   ioerr,  jcha,   kdat,   
     &                  kinst,  ncha,   ninf   
      double precision  fktemp, garb,   KERX,   krx,    piosix, rmin,   
     &                  scale,  xintge, XKERX,  xkrx,   zero
      character*8       name
      logical           callsc, INBND
      common /cmkerf/   iinst,  jcha,   ncha
      save   /cmkerf/
      parameter        (zero   = 0.d+0)
      parameter        (piosix = 0.5235987d+0)
      external          KERX,   XKERX
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (opchr(10) .le. 0) then
        dspc(7)  = diam(2) - diam(1)
        return
      end if
      call OPNKRI(ioerr, 2, mrecl, ninf,  cfiles)
      if (ioerr .ne. 0) then
        call ERROR(42, nlog(2), msglev(2), ioerr, garb, 
     &      'open', fail(mfail))
        call ERROR(42, nlog(5), msglev(5), ioerr, garb, 
     &      'open', fail(mfail))
        return
      end if
c
c ##FIRST INITIALIZE THE DIAMETER VECTOR 
c
      rmin    = wmach(5)
      dspc(7) = (dspc(2) - dspc(1))/dble(ndim-1)
      diam(1) = dspc(1)
      do 40 idia = 2, ndim
	diam(idia) = diam(idia - 1) + dspc(7)
40    continue
c
      call DLOAD(mdat, zero, fk, 1)
      idat    = 0
      kinst   = 0
      iinst   = 0
      do 100 ipnt = 1, mdat
        if (inspnt(ipnt) .ne. 0) then
          idat   = idat + 1
          jcha   = iabs(ichpnt(ipnt))
          callsc = (idat .eq. 1)
          if (.not. callsc)
     &      callsc = (iinst .ne. iabs(inspnt(ipnt)))
          iinst  = iabs(inspnt(ipnt))
          ncha   = nch(iinst)
          name   = inme(iinst)
          if (callsc) then
c
c ##THEN WE'RE LOOKING AT THE DATA FROM A NEW INSTRUMENT
c
            call SCRNM(2, nlog(5), msglev(5), 0, zero, name)
            kinst = kinst + 1
            call DCOPY(npar(iinst), tpar(kinst,1), minsiv, 
     &                 par(iinst,1), minst)
          end if
          call RCHB  (fail,   ibnd,   iinst,  jcha, mfail, 
     &                msglev(2), nlog(2),   ninf,   work)
c
c ##IF AN ERROR OCCURRED WHILE TRYING TO READ THE SUPPORT BOUNDS
c ##THEN WE NEED TO EXIT, ELSE WE SHOULD "PAD" THE INTERVALS READ
c ##BY RCHB
c
          if (fail(mfail) .gt. 0) then
            go to 1000
          else if (ibnd .gt. 0) then
            call INTVPD(ibnd,  tol(8), work,   2)
          end if
          xintge = -1.d+0
          do 110 idim = 1, ndim - 1
            idimp1 = idim + 1
            if (INBND(ibnd, diam(idimp1), diam(idim), work)) 
     &        then
              call DINTG(KERX,  ierr, mintg, krx,  diam(idim), 
     &                   diam(idimp1), tol(6), garb, rmin)
              call DINTG(XKERX, ierx, mintg, xkrx, diam(idim), 
     &                   diam(idimp1), tol(6), garb, rmin)
              fktemp = (xkrx - krx*diam(idim)) / dspc(7)
              fk(ipnt,idimp1) = fktemp
              fk(ipnt,idim)   = fk(ipnt,idim) + krx - fktemp
              if (ierr .ne. 0) then
                xintge = dble(ierr)
              else if (ierx .ne. 0) then
                xintge = dble(ierx)
              end if
            else
              fk(ipnt,idimp1) = zero
            end if
110       continue
          if (xintge .gt. rmin)
     &      call WARNG(19, nlog(2), msglev(2), jcha, xintge, name)
        end if
100   continue
c
c ##RESCALE THE KERNEL FUNCTIONS IF THE MASS DISTRIBUTION
c ##IS THE DESIRED RESULT
c
      if (opchr(11) .eq. 3) then
        do 200 idim = 1, ndim
          scale = 1.d+0 / piosix / dexp(diam(idim)*3.d+0)
          call DSCAL(mdat, scale, fk(1,idim), 1)
200     continue
      end if
1000  if (fail(mfail) .gt. 0) then
        call UINDEX(fail(1))
      else if (msglev(2) .ge. 60) then
c
c ##OUTPUT MATRIX REPRESENTING KERNEL FUNCTIONS
c
        write(nlog(2), 1) npdat, ndim
        do 500 kdat = 1, mdat
          if (inspnt(kdat) .gt. 0) then
            write(nlog(2), 2) kdat
            write(nlog(2), 3) (fk(kdat,idim), idim = 1, ndim)
          end if
500     continue
      end if
      close(ninf)
      return
1     format(//, ' Matrix which represents kernel function:',/,
     &      ' The kernel function matrix is ',i4,' by ',i4)
2     format(/,' Row ', i3, ' of kernel matrix . . .')
3     format(6(2x,d11.4))
*
*___end of CWK___*
*
      end
      logical function INBND(ibnd,     dhigh,  dlow,   xbnd)
      integer                ibnd
      double precision       xbnd,  dlow,   dhigh
      dimension              xbnd(2,*)
************************************************************************
*
* INBND is a logical function which is false when the interval defined
* by the first two arguments is outside all the intervals in defined in
* xbnd.
*
* Inputs:dlow,dhigh,ibnd,xbnd
*
* Outputs: INBND
*
* Called by :CWK
*
* 
* LAST MODIFIED: 29 OCT 1987
************************************************************************
      integer           jbnd
*-----------------------------------------------------------------------
      do 100 jbnd = 1, ibnd
        if ((dlow .ge. xbnd(2,jbnd)) .or. (dhigh .le. xbnd(1,jbnd))) 
     &    then
          INBND = .false.
        else
          INBND = .true.
          return
        end if
100   continue
      return
*
*___end of INBND___*
*
      end
      double precision function KERX(diam)
************************************************************************
*
* KERX returns the kernel function response for intrument iinst.
*
* INPUTS:
* diam  - the natural log of the partical diameter in microns
* iinst - the instrument listed in position iinst in micron.dat for which
*           the instrument response is desired.
*
* OUTPUTS:
*
* LAST MODIFIED: 5 JAN 1989
************************************************************************
      integer          iinst,  jcha,   minst,  mpar,   ncha 
      double precision diam,   INST1,  INST2,  INST3,  INST4,  INST5
      double precision INST6,  INST7,  INST8,  INST9,  INST10, tdiam
      common /cmkerf/  iinst,  jcha,   ncha
      save   /cmkerf/
      common /cmparm/  minst,  mpar
      save   /cmparm/
*-----------------------------------------------------------------------
      tdiam = dexp(diam)
      KERX  = 0.d+0
      if (iinst .le. 5) then
        if      (iinst .eq. 1) then
          KERX = INST1(tdiam)
        else if (iinst .eq. 2) then
          KERX = INST2(tdiam)
        else if (iinst .eq. 3) then
          KERX = INST3(tdiam)
        else if (iinst .eq. 4) then
          KERX = INST4(tdiam)
        else if (iinst .eq. 5) then
          KERX = INST5(tdiam)
        end if
      else
        if (iinst .eq. 6) then
          KERX = INST6(tdiam)
        else if (iinst .eq. 7) then
          KERX = INST7(tdiam)
        else if (iinst .eq. 8) then
          KERX = INST8(tdiam)
        else if (iinst .eq. 9) then
          KERX = INST9(tdiam)
        else if (iinst .eq. 10) then
          KERX = INST10(tdiam)
        end if
      end if
      return
*
*___end of KERX___*
*
      end
      subroutine CMPDM (fail,   ichpnt, inspnt, ldim,   mdat,   mdim,   
     &                  mfail,  mltdim, mrecl,  msglev, ndim,   nlog,   
     &                  npdat,  dspc,   wmach,  work,   cfiles)
      integer           fail,   ichpnt, inspnt, ldim,   mdat,   mdim,   
     &                  mfail,  mltdim, mrecl,  msglev, ndim,   nlog,   
     &                  npdat  
      double precision  dspc,   wmach,  work
      character*(*)     cfiles
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nlog(*), dspc(*), wmach(*), work(*), cfiles(*)
************************************************************************
*
* CMPDM chooses the dimension of the solution vector by looking at the
* difference of the mean diameters of the kernel functions.
*
* The dimension of the solution vector cannot be too small or
* some of the 'structure' of the size distribution may be lost, and too
* large of a size distribution will be too costly.  
* If "ndim">5 the calculation is skipped.
*
* LAST MODIFIED:10 JAN 1989
************************************************************************
      integer           IDAMIN, imin,   kdim,   msglog, nlog2
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (ndim    .ge. 5) then
        if (ndim .lt. npdat) 
     &    call WARNG(3, nlog2, msglog, npdat, 0.d+0, ' ')
        return
      end if
      msglog = msglev(2)
      nlog2  = nlog(2)
      call RDMEAN(fail,   ichpnt, inspnt, mdat,   mfail,  mrecl,  
     &            msglog, nlog2,  work,   cfiles, 'MEAN')
      if (fail(mfail) .gt. 0) then
        call ERROR(91, nlog2, msglog, 0, 0.d+0, ' ', fail(1))
        return
      end if
      call DBSORT(npdat, work)
      call DXMY(npdat-1, work(2), 1, work, 1) 
      imin = IDAMIN( npdat-1, work, 1 )
      if (dabs(work(imin)) .gt. wmach(5)) then
        kdim = (dspc(2) - dspc(1))/dabs(work(imin)) * mltdim
        call SCRNM(41, nlog2, msglog, kdim, 0.d+0, ' ')
      else
        kdim = mdim
        call SCRNM(42, nlog2, msglog, 0, 0.d+0, ' ')
      end if
      ndim = kdim
      if (ndim .lt. ldim) ndim = ldim
      if (ndim .gt. mdim) ndim = mdim
      if (kdim .ne. ndim) 
     &    call SCRNM(43, nlog2, msglog, ndim, 0.d+0, ' ')
      return
*
*___end of CMPDM___*
*
      end
      subroutine DBSORT(n, x)
      integer           n
      double precision  x
      dimension         x(*)
************************************************************************
*
* DBSORT sorts the vector x with the bubble sort algorithm.
*
* LAST MODIFIED: 5 NOV 1987
************************************************************************
      integer           iech,   jech
      double precision  temp
      logical           exch
*----------------------------------------------------------------------
      do 100 iech = n, 2, -1
	exch = .false.
	do 110 jech = 1, iech - 1
	  if (x(jech) .gt. x(jech + 1)) then
	    exch      = .true.
            temp      = x(jech)
            x(jech)   = x(jech+1)
            x(jech+1) = temp
	  end if
110     continue
	if (.not.exch) return
100   continue
      return
*
*___end of SORTBD___*
*
      end
      subroutine RDMEAN(fail,   ichpnt, inspnt, mdat,   mfail,  mrecl,  
     &                  msglog, nlog,   work,   cfiles, job)
      integer           fail,   ichpnt, inspnt, mdat,   mfail,  mrecl,  
     &                  msglog, nlog
      double precision  work
      character*(*)     cfiles, job
      dimension         fail(*), ichpnt(*), inspnt(*), work(*), 
     &                  cfiles(*)
************************************************************************
*
* RDMEAN reads from the data file micron.kri the log mean diameter of
* each of the kernel functions (job = MEAN), or the integral of
* the kernel functions (job = INTG).
*
* LAST MODIFIED: 5 NOV 1987
************************************************************************
      integer           ibnd,   ibndrc, icha,   idat,   inst,   ipnt,   
     &                  ioerr,  irchnl, joerr,  ncha,   ninf,   nrec
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      call OPNKRI(ioerr, 2, mrecl, ninf,   cfiles)
      idat = 0
      do 100 ipnt = 1, mdat
        if (ichpnt(ipnt) .gt. 0) then
          idat = idat + 1
          inst = inspnt(ipnt)
          icha = ichpnt(ipnt)
          read(ninf, rec = inst, err=1000, iostat=joerr) irchnl,ncha
          if (icha .gt. ncha) then
            call ERROR(46, nlog, msglog, inst, 0.d00, ' ', fail(mfail))
          else if (job(1:4) .eq. 'MEAN') then
            nrec = irchnl + 2*(icha - 1) + 1
            read(ninf, rec = nrec, err = 1000, iostat=joerr) work(idat)
          else if (job(1:4) .eq. 'INTG') then
            nrec = irchnl + 2*(icha - 1) 
            read(ninf, rec = nrec, err = 1000, iostat = joerr) 
     &        ibndrc, ibnd, work(idat)
          end if
        end if
100   continue
c
1000  if (ioerr .ne. 0)
     &  call ERROR(42, nlog, msglog, ioerr, 0.d00, 'open', fail(mfail))
      if (joerr .ne. 0)
     &  call ERROR(42, nlog, msglog, joerr, 0.d00, 'read', fail(mfail))
      close(ninf)
      return
*
*___end of RDMEAN___*
*
      end
