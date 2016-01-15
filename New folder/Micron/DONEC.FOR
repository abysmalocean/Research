      subroutine DONEC (fail,   ispc,   iter,   itersv, itnsmp, msglev, 
     &                  ndim,   nlog,   opchr,  drvhst, drvx,   drvy,   
     &                  f,      f2,     sphst,  tol,    done,   weak)
      integer           fail,   ispc,   iter,   itersv, itnsmp, msglev, 
     &                  ndim,   nlog,   opchr
      double precision  drvhst, drvx,   drvy,   f,      f2,     sphst,  
     &                  tol
      logical           done,   weak
      dimension         fail(*), ispc(*), msglev(*), nlog(*), opchr(*), 
     &                  drvhst(0:*), drvx(*), drvy(*), f(*), f2(*), 
     &                  sphst(0:*), tol(*)
************************************************************************
* DONEC determines if the search for the lambda that minimizes 
* V(lambda) and R(lambda) subject to constraints should be terminated.
* Some simple bookkeeping is done for FSPC. 
* If smoothing is chosen by constrained gcv, then weak solutions are
* prevented.
*
* DONEC returns done as true if
* 1) fail(1) > 0
* 2) iter = itnsmp
* 3) delta(lambda) < tol(2)  
* Condition 2 causes a warning to be printed.
*
* LAST MODIFIED:12 DEC 1988
************************************************************************
      double precision  dlgdlt, dlgsp,  rvli,   two,    xm,     zero
      double precision  deltsp, drvmax, smpmax
      logical           bound
      parameter        (two    = 2.d+0)
      parameter        (zero   = 0.d+0)
c-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
c
c ##CHECK IF THE MINIMUM IS BOUND.  WHILE TRYING TO BOUND THE MINIMUM,
c ##drvx(2) IS THE SMALLER LOG(LAMNDA), AND dxrv(3) IS THE LARGER ON 
c ##EXIT
c
      itersv  = iter
      rvli    = drvhst(iter)
      dlgsp   = dlog(sphst(iter))
      drvy(1) = rvli
      done    = .false.
      if (iter .eq. 1) then
        smpmax = sphst(iter)
        drvmax = drvhst(iter)
      else if (sphst(iter) .gt. smpmax) then
        smpmax = sphst(iter)
        drvmax = drvhst(iter)
      end if
c
c ##CHECK IF WE CAN CONTINUE SUCC. SUBST.
c
      if ( ispc(1) .le. 0 ) then
        if ( ispc(1) .ne. 0 ) ispc(1) = 1
        if ( (iter .ge. 2) .and. ( drvhst(iter) .ge. drvhst(iter-1) ) )
     &    ispc(1) = 1
        if ( dlog10(sphst(iter)/sphst(iter+1)) .ge. 2.d+0 ) ispc(1) = 1
      else if (ispc(1) .gt. 0) then
        ispc(1) = ispc(1) + 1
        if (weak .and. (opchr(13) .eq. 0)) then
          smpmax = sphst(iter)
          drvmax = drvhst(iter)
          ispc(1) = 1
          ispc(2) = 0
        end if
      end if
c
c ##IF SUCC. SUBST. IS BEING USED TO FIND THE SMOOTHING PARAMETER
c
      if (ispc(1) .eq. 0) then
        deltsp = dabs(dlog(sphst(iter)) - dlog(sphst(iter+1))) 
        done   = (deltsp .lt. tol(2))
      else 
        if (ispc(1) .eq. 1) then 
          drvx(2) = dlog(smpmax)
          drvy(2) = drvmax
          drvx(3) = dlog(smpmax)
          drvy(3) = drvmax
        else if ( ispc(1) .eq. 2 ) then
          call REORD (dlgsp, drvx,   drvy,   rvli)
          done = (drvx(3) - drvx(2)) .lt. tol(2)
        else if ( ispc(2) .eq. 0 ) then
          bound = ((dlgsp .gt. drvx(3)) .and. (rvli .gt. drvy(3)))
     &      .or.  ((dlgsp .lt. drvx(2)) .and. (rvli .gt. drvy(2)))
          if (bound) then
c
c ##THEN PREPARE FOR THE MINIMIZATION ROUTINE
c
            ispc(2) = 1
            drvx(7) = zero
            drvx(8) = zero
            if (dlgsp .lt. drvx(2)) then
              drvy(4) = drvy(2)
              drvx(4) = drvx(2)
              drvx(5) = dlgsp
              drvx(6) = drvx(3)
            else
              drvy(4) = drvy(3)
              drvx(4) = drvx(3)
              drvx(5) = drvx(2)
              drvx(6) = dlgsp
            end if
            drvx(2) = drvx(4)
            drvx(3) = drvx(4)
            drvy(2) = drvy(4)
            drvy(3) = drvy(4)
          else
            call REORD (dlgsp, drvx,   drvy,   rvli)
            done = (drvx(3) - drvx(2)) .lt. tol(2)
          end if
        else
          ispc(2) = ispc(2) + 1
          call DBMORD(drvx,   drvy)
        end if
c
c ##SEE IF WE'RE DONE
c
        xm     = (drvx(5) + drvx(6)) / two
        dlgdlt = two * tol(2) - (drvx(6) - drvx(5)) / two
        done   = done .or. ( (ispc(2) .ge. 1) .and. 
     &                       (dabs(drvx(4) - xm) .le. dlgdlt) )
      end if
      if ( ( .not. done ) .and. ( iter .ge. itnsmp ) ) then
        done = .true.
        call WARNG(30, nlog(2), msglev(2), itnsmp, 0.d+0, ' ')
      end if
      call DCOPY(ndim, f, 1, f2, 1)
      return
c
*___End of DONEC___*
c
      end
      subroutine REORD (dlgsp,  drvx,   drvy,   rvl)
      double precision  dlgsp,  drvx,   drvy,   rvl
      dimension         drvx(*), drvy(*)
************************************************************************
* REORD is a bookeeping utility used only by DONEC.
*
* LAST MODIFIED:12 DEC 1988
************************************************************************
      if (dlgsp .lt. drvx(2)) then
        drvx(3) = drvx(2)
        drvx(2) = dlgsp
        drvy(3) = drvy(2)
        drvy(2) = rvl
      else if (dlgsp .gt. drvx(3)) then
        drvx(2) = drvx(3)
        drvx(3) = dlgsp
        drvy(2) = drvy(3)
        drvy(3) = rvl
      else if (drvy(2) .gt. drvy(3)) then
        drvy(2) = rvl
        drvx(2) = dlgsp
      else
        drvy(3) = rvl
        drvx(3) = dlgsp
      end if
      return
c
*___end of REORD___*
c
      end
