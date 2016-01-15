      subroutine OUTPUT(fail,   ichpnt, iconrc, inspnt, itersv, mctotl, 
     &                  mdat,   mdim,   merr,   mfail,  minsiv, minst,  
     &                  miter,  mpar,   msglev, nch,    ndim,   ninst,  
     &                  nlog,   npar,   npdat,  opchr,  bl,     bspc,   
     &                  bu,     diam,   drvhst, dspc,   einv,   f2,     
     &                  objvec, par,    solspc, sphst,  tol,    tpar,   
     &                  volx,   wmach,  y,      ypctdf, yr,     header, 
     &                  inme,   parnm,  smooth)
      integer           fail,   ichpnt, iconrc, inspnt, itersv, mctotl, 
     &                  mdat,   mdim,   merr,   mfail,  minsiv, minst,  
     &                  miter,  mpar,   msglev, nch,    ndim,   ninst,  
     &                  nlog,   npar,   npdat,  opchr
      double precision  bl,     bspc,   bu,     diam,   drvhst, dspc,   
     &                  einv,   f2,     objvec, par,    solspc, sphst,  
     &                  tol,    tpar,   volx,   wmach,  y,      ypctdf, 
     &                  yr
      character*(*)     header, inme,   parnm,  smooth
      dimension         fail(*), ichpnt(*), iconrc(*), inspnt(*), 
     &                  msglev(*), nch(*), nlog(*), npar(*), opchr(*), 
     &                  bl(*), bspc(*), bu(*), diam(*), drvhst(0:*),  
     &                  dspc(*), einv(mdat,*), f2(*), objvec(0:*), 
     &                  par(*), solspc(*), sphst(0:*), tol(*), 
     &                  tpar(*), volx(*), wmach(*), y(*), ypctdf(*), 
     &                  yr(*), header(*), inme(*), parnm(*), smooth(*)
************************************************************************
*
* OUTPUT lives up to it's name and writes the requested amount of 
* output to the "root".out data file.  The  integer  MSGOUT  controls
*  the amount of output according to the following rule: 
*
*      =   0  :the .out file is not opened
*
*      =  10  :header plus solution vector
*
*      =  20  :same as 10 with solution properties (total number, 
*              variance,  mean diamter,  roughness,  discrepancy, 
*              and value of the functional minimized).
*
*      =  30  :same as 20 with an  echo of the  input   instrument 
*              reading, the inverted readings, percent differences,
*              and bounds if any.
*
*      =  40  :same  as 30 with corresponding solution  diameter, 
*              and dV/dln(diam)
*
*      =  50  :same as 40 with an echo of the computation level, 
*              inversion interval, and dimension of solution vector
*
* LAST MODIFIED:16 DEC 1988
************************************************************************
      integer           icon,   icono,  idat,   idim,   inst,   instru, 
     &                  ipnt,   jpar,   kdim,   kinst,  msgout, nclin,  
     &                  nout,   numpar
      double precision  ed,     ed1,    ed2,    edvar1, edvar2 
      character*80      datlin
      character*9       chunit
      character*7       distyp
      logical           DOUTSD, sol
*-----------------------------------------------------------------------
      if ((fail(2) .le. 0) .and. (fail(3) .le. 0)) then
        if (msglev(2) .gt. 0)  close(nlog(2))
        sol = (fail(1)  .le. 0) .or. (fail(15) .gt. 0) .or. 
     &        (fail(16) .gt. 0) .or. (fail(17) .gt. 0) .or. 
     &        (fail(18) .gt. 0)
        msgout = msglev(3)
        nout   = nlog(3)
      else
        return
      end if
      if (opchr(11) .eq. 1) then
        chunit = ' (#/cc)  '
        distyp = ' Number'
      else
        chunit = '(um^3/cc)'
        distyp = ' Volume'
      end if
c
c ##WRITE THE "root".plt FILE
c
      if (sol .and. (msglev(4) .ge. 10)) then
        do 300 idim = 1, ndim
          if (diam(idim) .gt. dspc(6)) then
            kdim = idim - 1
            if (kdim .lt. 1) kdim = 1
          else if (diam(idim) .lt. dspc(5)) then
            kdim = idim + 1
            if (kdim .gt. ndim) kdim = ndim
          else
            kdim = idim
          end if
          if (.not. DOUTSD(diam(kdim), dspc(5), dspc(6))) then
            ed  = dexp(diam(idim))
            if (msglev(4) .ge. 30) then
              write(nlog(4),1) ed, f2(idim), volx(idim)
            else if (msglev(4) .ge. 20) then
              write(nlog(4),1) ed, f2(idim)
            else
              write(nlog(4),1) f2(idim)
            end if
          end if
300     continue 
      end if
c
c ##WRITE THE OUTPUT FILE
c
      if (msgout .ge. 10) then
        call SCRNM(27, nlog(5), msglev(5), 0, 0.d+0, ' ')
        write(nout,104) header(1), header(2)
c
c ##OUTPUT SOLUTION PROPERTIES
c
        if (msgout .ge. 20) then
          write(nout, 201)
          if (sol) then
            write(nout, 202) chunit, solspc(4) 
            write(nout, 206) solspc(1), solspc(2), solspc(3), 
     &                       smooth(iabs(opchr(1))), sphst(itersv), 
     &                       tol(2), drvhst(itersv), objvec(itersv)
          else
            write(nout, 203) smooth(iabs(opchr(1))), sphst(itersv)
          end if
          edvar1 = dexp(dspc(5))
          edvar2 = dexp(dspc(6))
          write(nout, 207) edvar1, edvar2
        end if
c
c ##NOW WE OUTPUT THE SOLUTION VECTOR
c
        if (msgout .ge. 40) then
          if (msgout .ge. 50) then
            ed1 = dexp(diam(1))
            ed2 = dexp(diam(ndim))
            write(nout, 502) ndim
            write(nout, 506) 'Lower', ed1
            if (bspc(1) .gt. -wmach(5))
     &        write(nout, 507) bspc(1), 'lower'
            write(nout, 506) 'Upper', ed2
            if (bspc(2) .gt. -wmach(5))
     &        write(nout, 507) bspc(2), 'upper'
            if (bspc(3) .gt. -wmach(5))
     &        write(nout, 205) 'Lower', chunit, bspc(3)
            if (bspc(4) .gt. -wmach(5))
     &        write(nout, 205) 'Upper', chunit, bspc(4)
          end if
          if (sol) then
            write(nout, 401)
            do 100 idim = 1, ndim
              ed = dexp(diam(idim))
              write(nout, 402) ed, f2(idim), volx(idim)
100         continue
          end if
        else if (sol) then
          write(nout, 101)
          do 150 idim = 1, ndim
            write(nout, 102) f2(idim)
150       continue
        else
          write(nout, 103)
        end if
c
c ##ECHO THE INSTRUMENT DATA AND INVERTED DATA
c
        if (msgout .ge. 30) then
          instru = 0
          nclin  = 0
          kinst  = 0
          idat   = 0
          write(nout, 301)
          do 200 ipnt = 1, mdat
            if (inspnt(ipnt) .gt. 0) then
              inst = inspnt(ipnt)
              idat = idat + 1
c
c ##IF WE HAVE A NEW INSTRUMENT, WRITE IT'S NAME AND PARAMETERS
c
              if (inst .ne. instru) then
                instru = inst
                kinst  = kinst + 1
                call DCOPY(npar(inst), tpar(kinst), minsiv, 
     &                     par(inst), minst)
                write(nout, 302) inme(inst)
                numpar = npar(inst)
                if (numpar .gt. 0) 
     &            write(nout, 303) (parnm(jpar), par(jpar), jpar = inst, 
     &                              (numpar - 1)*minst + inst, minst)
                write(nout, 305)
              end if
c
c ##SO MANY IF STATEMENTS TO CHECK IF THE DATA HAS BOXED CONSTRAINTS
c ##COUPLED WITH AN INCOMPETENT PROGRAMMER
c
              icon = iconrc(idat)
              if (icon .ne. 0) then
                nclin = nclin + 1
                icono = nclin + ndim
              end if
c
c ##WRITE THE INVERTED DATA
c
              datlin = ' '
              write(datlin(3:6), '(i3)') ichpnt(ipnt)
              write(datlin(10:19), '(e10.3)') y(idat)
              if (sol) write(datlin(23:32), '(e10.3)') yr(idat)
              if (sol) write(datlin(36:45), '(e10.3)') ypctdf(idat)
              if ((icon .eq. 1) .or. (icon .eq. 3)) 
     &          write(datlin(51:60), '(e10.3)') bl(idat)
              if (icon .gt. 2) write(datlin(64:73), '(e10.3)') bu(idat)
              write(nout, '(a)') datlin(1:73)
            end if
200       continue
        end if
      end if
      if (msgout .gt. 0) close(nout)
      if (msglev(4) .gt. 0) close(nlog(4))
      return
c-----------------------------------------------------------------------
1     format(3(2x,d13.4))
101   format(//' Number distribution:',/,
     &       ' dN/dlog[dia.]'/,3x,'(1/cm^3)')
102   format(d11.4)
103   format(' Approximate solution not computed.')
104   format(//,1x,a80,/,1x,a80)
201   format(//' Solution Properties:',/)
202   format(' Integral ',a9,t42, 
     &       '=', t50,d11.4)
203   format(' Smoothing parameter (',a4,')',t42, '=', t50,d11.4)
205   format(/,1x,a5,' concentration constraint ',a9,t42,'=', t50, 
     &       d11.4)
206   format(
     &       ' Log-mean particle diameter (microns)',t42, 
     &       '=', t50,d11.4,/,
     &       ' Geometric standard deviation',t42, '=', t50,d11.4,/,
     &       ' Integral of (f''''(x))^2',t42, '=', t50,d11.4,///,
     &       ' Smoothing parameter (',a4,')',t42, '=', t50,d11.4,/,
     &       ' Tolerance used to find target (tol(2))',t42,
     &       '=',t50,d11.4,/,
     &       ' Smoothing choice functional',t42, '=', t50,d11.4,/,
     &       ' Quadratic functional ',t42, '=', t50,d11.4)
207   format(
     &       ' Interval of kernel function variance',t42,'=',
     &       t49,'(',d11.4,',',d11.4,')')
401   format(//' Number distribution and volume distribution:',//,
     &       '   dia. ',t30,'dN/dlog[dia.]',t60,
     &       'dV/dlog[dia.]',/' (microns)',t33,'(1/cm^3)',
     &       t58,'((microns/cm)^3)',/)
402   format(d11.4,t30,d11.4,t60,d11.4)
301   format(//,' Instrument readings:',//)
302   format(/,' Instrument: ',a8)
303   format(' Parameters:',/,
     &       100(1x,a8,' = ',d11.4,10x,a8,' = ',d11.4,/))
305   format(/,' Channel', t10, 'Input data', t22, 'inverted data', t38,
     &'rsdiff.', t51, 'lower bound', t64, 'upper bound')
501   format(//' Level of computation',t42, '=', t50,i11)
502   format(' Dimension of solution vector',t42, '=', t50,i11)
506   format(/,1x,a5,' interval diameter (microns)',t42, '=', t50,d11.4)
507   format(' Size distribution constrained to',t42,'=',t50,d11.4,/,
     &       ' at the ',a5,' interval diameter.')
c
*___end of OUTPUT___*
c
      end
      subroutine PREOUT(fail,   ichpnt, inspnt, istate, itersv, mdat,   
     &                  mdim,   mfail,  miter,  mrecl,  msglev, ndim,   
     &                  nlog,   npdat,  opchr,  p,      diam,   drvhst, 
     &                  dspc,   enorm,  f2,     fk,     fkhat,  objvec, 
     &                  solspc, sphst,  volx,   wmach,  work,   y,      
     &                  yhatm,  yhatr,  ypctdf, yr,     cfiles)
      integer           fail,   ichpnt, inspnt, istate, itersv, mdat,   
     &                  mdim,   mfail,  miter,  mrecl,  msglev, ndim,   
     &                  nlog,   npdat,  opchr,  p
      double precision  diam,   drvhst, dspc,   enorm,  f2,     fk,     
     &                  fkhat,  objvec, sphst,  solspc, volx,   wmach,  
     &                  work,   y,      yhatm,  yhatr,  ypctdf, yr
      character*(*)     cfiles
      dimension         fail(*), ichpnt(*), inspnt(*), istate(*), 
     &                  msglev(*), nlog(*), opchr(*), diam(*), 
     &                  drvhst(0:*), dspc(*), enorm(*), f2(*), 
     &                  fk(*), fkhat(*), objvec(0:*), 
     &                  sphst(0:*), solspc(*), volx(*), wmach(*), 
     &                  work(*), y(*), yhatm(*), yhatr(*), 
     &                  ypctdf(*), yr(*), cfiles(*)
************************************************************************
*
* PREOUT computes some of the solution properties written to the output
* file.
*
* VARIABLES:(type, input/output, array)
* f2       - the size distribution (solution)
* rmin     - the smallest positive machine number
* wmach(7) - the largest  positive machine number
* dlmean   - log-mean diameter of the distribution
* geostd   - geometric standard deviation of the distribution
* tnumbr   - particle concentration (#/cc)
* volume   - particle volume (microns/cm)^3
* volx     - volume size distribution
* yr       - the recovered data
*
* LAST MODIFIED:18 FEB 1988
************************************************************************
      integer           ffrst,  i,      icp,    ipos,   idat,   idat2,  
     &                  IDXGT,  ihigh,  ilow,   j,      mesgnm
      double precision  dltemp, DROUGH, dvart,  fintg,  one,    
     &                  piosix, rmin,   temp,   TRAP,   varian, zero
      logical           sol
      parameter        (one    = 1.d+0)
      parameter        (piosix = 0.52359877560d+0)
      parameter        (zero   = 0.d+0)
*-----------------------------------------------------------------------
      if (opchr(5) .le. 0) return
c
c ##FIRST REPORT THE COMPLETION STATUS OF THE INVERSION
c
      mesgnm = 30
      ffrst = 3
      if (fail(1) .gt. 0) 
     &  mesgnm = mesgnm + (ffrst - 1) + 
     &           IDXGT(mfail - ffrst - 1, 0, fail(ffrst), 1)
      call SCRNM(mesgnm, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
c
      sol = (fail(1)  .le. 0) .or. (fail(15) .gt. 0) .or. 
     &      (fail(16) .gt. 0) .or. (fail(17) .gt. 0) .or. 
     &      (fail(18) .gt. 0)
      if (.not. sol) return
      if (fail(18) .gt. 0) then
        rmin = wmach(5)
c
c ##RECOMPUTE yhatr AND rough
c
        call MMMLT(mdat, mdim, mdat, ndim, 1,  p, fkhat, f2, yhatr)
        solspc(3)  = DROUGH(ndim, dspc(7), f2, work)
      end if
c
c ##DETERMINE THE INTERVAL DEFINED BY THE VARIANCE OF THE KERNEL
c ##FUNCTIONS
c
      if ((opchr(10) .ne. 0) .and. (opchr(7) .ne. 0)) then
        call RDMNVR(fail,    ichpnt, inspnt, mdat, mfail, mrecl,  
     &              msglev(2), nlog(2), wmach,  work,  cfiles)
        temp  = dsqrt(work(2))
        dspc(5) = work(1) - temp
        dspc(6) = work(1) + temp
        do 100 idat = 2, npdat
          idat2 = idat + idat
          temp  = dsqrt(work(idat2)) 
          dvart = work(idat2-1) - temp
          if (dspc(5) .gt. dvart) dspc(5) = dvart
          dvart = work(idat2-1) + temp
          if (dspc(6) .lt. dvart) dspc(6) = dvart
100     continue
      end if
c
c ##COMPUTE THE RECOVERED DATA AND THE SCALED DIFFERENCE VECTOR
c
      idat = 0
      do 200 icp = 1, mdat
        if (inspnt(icp) .gt. 0) then
          idat = idat + 1
          call DCOPY(ndim, fk(icp), mdat, work(idat), mdat)
        end if
200   continue
      call MMMLT(mdat, mdim, mdat, ndim, 1,  npdat, work, f2, yr )
      do 300 idat = 1, npdat
        temp = y(idat) - yr(idat)
        if (temp .eq. zero) then
          ypctdf(idat) = zero
        else if (enorm(idat) .le. rmin) then
          ypctdf(idat) = -wmach(7)
        else
          ypctdf(idat) = (y(idat) - yr(idat)) * enorm(idat)
        end if
300   continue
c
c ##CHECK FOR COLLAPSING SOLUTION
c
      if (opchr(2) .le. 0) then
        ipos = 0
        do 400 idat = 1, npdat
          if (yr(idat) .ge. y(idat)) ipos = ipos + 1
400     continue
        if (ipos .lt. (npdat/10 + 1)) 
     &    call WARNG(28, nlog(2), msglev(2), 0, 0.d+0, ' ')
      end if
c
c ##CHECK FOR LARGE INTERVALS OF ZERO'S AT THE ENDPOINTS
c
      ihigh = ndim
      ilow  = 0
      do 500 i = 1, ndim
        j = ndim - i + 1
        if (istate(j) .eq. 0) ilow  = j - 1
        if (istate(i) .eq. 0) ihigh = i + 1
500   continue
      if ( dble(ilow + ndim - ihigh)/dble(ndim) .gt. 0.15d+0 )
     &  call SCRNM(8, nlog(2), msglev(2), 0, 0.d+0, ' ')
c
c ##MAKE SURE THE SOLUTION HAS NOT REACHED THE UPPER BOUND
c
c
c ##COMPUTE THE DISTRIBUTION INTEGRAL
c
      fintg = TRAP  (dspc(7),  ndim,   f2,     1 )
c
c ##COMPUTE THE LOG-MEAN DIAMETER 
c
      if (fintg .gt. rmin) then
        call     DCOPY (ndim,   f2,     1,      work,   1)
        call     DPSCAL(ndim,   1,      diam,   1,      work,   1)
        dltemp = TRAP  (dspc(7),  ndim,   work,   1 ) / fintg
      else
        dltemp = (diam(1) + diam(ndim))/2.d+0
      end if
      solspc(1) = dexp(dltemp)
c
c ##COMPUTE THE GEOMETRIC STANDARD DEVIATION
c
      if (fintg .gt. rmin) then
        call     DCOPY (ndim,   f2,    1,       work,   1 )
        call     DPSCAL(ndim,   2,     diam,    1,      work,   1)
        varian = TRAP  (dspc(7),  ndim,   work,   1 ) / fintg -
     &           dltemp**2
      else
        varian = ((diam(ndim) - diam(1))**2) / 12.d+0
      end if
      solspc(2) = dexp(dsqrt(varian))
c
c ##DETERMINE TOTAL NUMBER, VOLUME, AND CONVERT f2 TO VOLUME OR
c ##NUMBER DISTRIBUTION
c
      call     DCOPY (ndim,   f2,     1,      volx,   1 )
      if (opchr(11) .eq. 1) then
        call     DESCAL(ndim,   3.d+0,  diam,   1,      volx, 1 )
        call     DSCAL (ndim,   piosix, volx,   1 )
        solspc(4) = fintg
        solspc(5) = TRAP  (dspc(7),  ndim,   volx,   1 )
      else 
        call     DESCAL(ndim,  -3.d+0, diam,    1,      f2, 1)
        call     DSCAL (ndim,  (one/piosix), f2,  1)
        solspc(4) = TRAP  (dspc(7),  ndim,   f2,     1)
        solspc(5) = fintg
      end if
      return
*
*___end of PREOUT___*
*
      end 
      subroutine RDMNVR(fail,   ichpnt, inspnt, mdat,   mfail, mrecl,  
     &                  msglog, nlog,   wmach,  work,   cfiles)
      integer           fail,   ichpnt, inspnt, mdat,   mfail, mrecl,  
     &                  msglog, nlog
      double precision  wmach,  work
      character*(*)     cfiles
      dimension         wmach(*), fail(mfail), ichpnt(*), work(2,*), 
     &                  inspnt(*), cfiles(*)
************************************************************************
*
* RDMNVR reads from micron.kri to find out the mean diameter and 
* variance of the kernel function.  
*
* LAST MODIFIED: 5 FEB 1988
************************************************************************
      integer           icha,   idat,   inst,   ioerr,  ipnt,   irchnl, 
     &                  joerr,  ncha,   ninf,   nrec
*-----------------------------------------------------------------------
      call OPNKRI(ioerr, 2, mrecl, ninf,  cfiles)
      idat = 0
      do 100 ipnt = 1, mdat
        if (inspnt(ipnt) .gt. 0) then
          idat = idat + 1
          icha = ichpnt(ipnt)
          inst = inspnt(ipnt)
          read(ninf, rec = inst, err = 1000) irchnl, ncha
          if (icha .gt. ncha) then
            call ERROR(46, nlog, msglog, inst, 0.d00, ' ', fail(mfail))
          else
            nrec = irchnl + 2*(icha - 1) + 1
            read(ninf, rec = nrec, err = 1000, iostat = joerr) 
     &       work(1,idat), work(2,idat)
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
c
*___end of RDMNVR___*
c
      end
