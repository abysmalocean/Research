      subroutine CMINH (fail,   ispc,   istate, iter,   itnsmp, iwork,  
     &                  mdat,   mdim,   msglev, nclin,  ndim,   nlog,   
     &                  opchr,  p,      bl,     bu,     drvhst, dspc,   
     &                  fkhat,  smp,    sphst,  tnrmsq, tol,    wmach,  
     &                  work,   yhatm)
      integer           fail,   ispc,   istate, iter,   itnsmp, iwork,  
     &                  mdat,   mdim,   msglev, nclin,  ndim,   nlog,   
     &                  opchr,  p
      double precision  bl,     bu,     drvhst, dspc,   fkhat,  smp,    
     &                  sphst,  tnrmsq, tol,    wmach,  work,   yhatm
      dimension         fail(*), ispc(*), istate(*), iwork(*), 
     &                  msglev(*), nlog(*), opchr(*), bl(*), bu(*), 
     &                  drvhst(0:*), dspc(*), fkhat(mdat,*), sphst(0:*), 
     &                  tol(*), wmach(*), work(mdim,*), yhatm(*)
************************************************************************
* CMIN computes V(lambda) or R(lambda) given a set of constraints
* on the solution.  If (ispc(1) = 0) (successive substitution) is true
* then CMIN computes the next smoothing parameter for FSPC and stores 
* it in sphst(iter+1).
*
* LAST MODIFIED: 1 MAR 1989
* HOMOGENEOUS VERSION, NCLIN = 0
************************************************************************
      integer           i,      inform, iy,     nfdim,  norm
      double precision  dlgcst, RLAM,   rvmin,  sigsq,  VLAM,   xmin
      logical           sigkn
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (nclin .ne. 0) 
     &  write(nlog(5),*) ' CMIN temporary warning: nclin .ne. 0'
      iy     = mdim + mdim + 5
      nfdim  = 0
      sigsq  = tnrmsq / dble(p)
      norm   = opchr(12)
c
c ##FIRST REMOVE THE BOUNDED COLUMNS FROM fkhat
c
      call DCOPY(p, yhatm, 1, work(1,iy), 1)
      do 100 i = 1, ndim
        if (istate(i) .eq. 0) then
          nfdim = nfdim + 1
          call DCOPY(p, fkhat(1,i), 1, work(1,nfdim), 1)
        else if (istate(i) .eq. 1) then
          call DAXPY(p, bl(i), fkhat(1,i), 1, work(1,iy), 1)
        else
          call DAXPY(p, bu(i), fkhat(1,i), 1, work(1,iy), 1)
        end if
100   continue
      call PRERV (fail,   mdat,   mdim,   msglev, nfdim,  nlog,   
     &            norm,   p,    sigsq,  wmach,  work,   work(1,iy))
      if (fail(1) .gt. 0) return
      sigkn = (iabs(opchr(1)) .eq. 3)
      dlgcst = dlog( dble(p) / (dspc(7)**(norm + norm - 1)))
      work(1,1) = dlog(smp) + dlgcst
      if (sigkn) then
        drvhst(iter) = RLAM(work)
      else
        drvhst(iter) = VLAM(work)
      end if
c
c ##COMPUTE THE NEXT SMOOTHING PARAMETER FOR FSPC, AND STORE IN
c ##sphst(2)
c
      if (ispc(1) .eq. 0) then
        call MINRV( inform, msglev, nlog,   p,      dlgcst, dspc,   
     &              rvmin,  tol,    wmach,  work,   xmin,   sigkn)
        if (inform .ne. 0) ispc(1) = -1
        sphst(iter+1) = dexp(xmin - dlgcst) 
      end if
      return
*
*___end of CMIN___*
*
      end
