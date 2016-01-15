      subroutine OUTSPC(fail,   ichpnt, iconrc, inspnt, itersv, kerr,   
     &                  kfiles, mctotl, mdat,   mdim,   merr,   mfail,  
     &                  minsiv, minst,  mintg,  miter,  mpar,   msglev, 
     &                  nch,    ndim,   ninst,  nlog,   npar,   npdat,  
     &                  opchr,  p,      bl,     bspc,   bu,     diam,   
     &                  drvhst, dspc,   e,      einv,   enorm,  f,      
     &                  f2,     fk,     objvec, par,    solspc, sphst,  
     &                  tol,    tpar,   volx,   wmach,  work,   y,      
     &                  yhatm,  yhatr,  yinf,   ypctdf, yr,     fname,  
     &                  header, inme,   inpnm,  parnm,  tdstnm, tfnme)
      integer           fail,   ichpnt, iconrc, inspnt, itersv, kerr,   
     &                  kfiles, mctotl, mdat,   mdim,   merr,   mfail,  
     &                  minsiv, minst,  mintg,  miter,  mpar,   msglev, 
     &                  nch,    ndim,   ninst,  nlog,   npar,   npdat,  
     &                  opchr,  p
      double precision  bl,     bspc,   bu,     diam,   drvhst, dspc,   
     &                  e,      einv,   enorm,  f,      f2,     fk,     
     &                  objvec, par,    solspc, sphst,  tol,    tpar,   
     &                  volx,   wmach,  work,   y,      yhatm,  yhatr,  
     &                  yinf,   ypctdf, yr
      character*(*)     fname,  header, inme,   inpnm,  parnm,  tdstnm, 
     &                  tfnme
      dimension         fail(*), ichpnt(*), iconrc(*), inspnt(*), 
     &                  msglev(*), nch(*), nlog(*), npar(*), opchr(*), 
     &                  bl(*), bspc(*), bu(*), diam(*), drvhst(0:*), 
     &                  dspc(*), e(mdat,*), einv(mdat,*), enorm(*), 
     &                  f(*), f2(*), fk(mdat,*), objvec(0:*), 
     &                  par(minst,*), sphst(0:*), solspc(*), tol(*), 
     &                  tpar(minsiv,*), volx(*), wmach(*), work(*), 
     &                  y(*), yhatm(*), yhatr(*), yinf(*), ypctdf(*), 
     &                  yr(*), inme(*), header(*), parnm(*)
************************************************************************
* OUTSPC is the specialized output subroutine called by MICRON.  The 
* purpose of this subroutine is to provide the user with some means 
* of controlling the output in special circumstances.  
*
* LAST REVISION:11 JAN 1988
************************************************************************
      integer           i,      j,      k(4),  nfiles
      double precision  DDOT,   one,    zero,  scale
      double precision  TSTFN, dmean,   smin,  smax
      character*8       infout(4)
      parameter        (zero = 0.d+0)
      parameter        (one  = 1.d+0)
      save              k, scale
      data   infout    /'idd1.plt', 'idd2.plt', 'idd3.plt', 'idd4.plt'/
*-----------------------------------------------------------------------
      nfiles = 100 
      smin   = dsqrt(0.001d+0)
      smax   = 1.d+0
      scale  = dexp(dble(kfiles-1) * dlog(smax/smin) / dble(nfiles-1) )
     &       * smin
      scale  = scale * scale
      if (kfiles .eq. 1) then
        do 100 i = 1, 4
          call FUNIT(k(i))
          open(k(i), file = infout(i))
100     continue
      end if
      do 300 i = 1, 3
        write(k(i),1) scale, yinf(ndim+3+i)
300   continue
      write(k(4),1) scale, yinf(ndim+1)
1     format(2(2x,d10.3))
*     par(5,5) = par(5,5) * scale
      return
*
*___end of OUTSPC___*
*
      end
