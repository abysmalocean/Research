      subroutine INPSPC(fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,
     &                  kfiles, mctotl, mdat,   mdim,   merr,   mfail,  
     &                  minsiv, minst,  mpar,   msglev, nch,    ndim,   
     &                  ninp,   ninst,  nlog,   npar,   npdat,  opchr,  
     &                  skip,   bl,     bspc,   bu,     dspc,   e,      
     &                  par,    smp,    tpar,   wmach,  work,   y,      
     &                  header, fname,  inme)
      integer           fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,
     &                  kfiles, mctotl, mdat,   mdim,   merr,   mfail,  
     &                  minsiv, minst,  mpar,   msglev, nch,    ndim,   
     &                  ninp,   ninst,  nlog,   npar,   npdat,  opchr,  
     &                  skip
      double precision  bl,     bspc,   bu,     dspc,   e,      par,    
     &                  smp,    tpar,   wmach,  work,   y
      character*(*)     header, fname,  inme
      dimension         fail(*), ichpnt(*), inspnt(*), msglev(*), 
     &                  nch(*), nlog(*), npar(*), opchr(*), bl(*), 
     &                  bspc(*), bu(*), dspc(*), e(mdat,*), 
     &                  par(minst,*), tpar(minsiv,*), wmach(*), 
     &                  work(*), y(*), header(2), inme(*)
************************************************************************
*
*
* LAST MODIFIED:17 AUG 1988
************************************************************************
      integer           idat,   nclear, nseed,   nfiles
      double precision  DDOT,   DNRM2,   DSUM,   one,    scale,  zero
      double precision  smin,   smax,    cscale
      logical           first
      double precision  eind,    ysave(100), beta(12)
      parameter        (zero = 0.d+0)
      parameter        (one  = 1.d+0)
*-----------------------------------------------------------------------
      if (kfiles .eq. 1) call DCOPY(npdat, y, 1, ysave, 1)
      call DCOPY(npdat, ysave, 1, y, 1)
      ndim   = 151
      ierrmd = 2
      kerr   = npdat
      nfiles = 100 
      smin   = dsqrt(0.001d+0)
      smax   = 1.d+0
      scale  = dexp(dble(kfiles-1) * dlog(smax/smin) / dble(nfiles-1) )
     &       * smin
      cscale = dsqrt(smax * smax - scale * scale)
      call DLOAD(mdat*mdim, zero, e, 1)
      call DCOPY(npdat, y, 1,   e, mdat+1)
      call DSCAL(npdat, scale,  e, mdat+1)
      call DCOPY(npdat, y, 1,   e(1,npdat), 1)
      call DSCAL(npdat, cscale, e(1,npdat), 1)
      e(npdat,npdat) = y(npdat)
      return
*
*___end of INPSPC___*
*
      end
