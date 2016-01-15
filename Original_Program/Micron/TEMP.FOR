      subroutine INVEK (fail,   ichpnt, ierrmd, mdat,   msglev, ndim,   
     &                  nlog,   npdat,  p,      einv,   fk,     fkhat,  
     &                  wmach,  work,   y,      yhatm)
      integer           fail,   ichpnt, ierrmd, mdat,   msglev, ndim,   
     &                  nlog,   npdat,  p
      double precision  einv,   fk,     fkhat,  y,      yhatm,  wmach,  
     &                  work
      dimension         fail(*), ichpnt(*), msglev(*), nlog(*), 
     &                  einv(mdat,*), fk(*), fkhat(*), wmach(*), 
     &                  work(*), y(*), yhatm(*)
************************************************************************
*
* LAST MODIFIED:19 OCT 1989
************************************************************************
      integer           idat,   ip,     isave
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      call SCRNM(25, nlog(5), msglev(5), 0, 0.d+0, ' ')
c
c ##FIRST COPY THE USED ROWS OF fk TO work
c
      isave = 1
      idat  = 1
10    if ((isave .le. npdat) .or. (idat .le. mdat)) then
        if (ichpnt(idat) .gt. 0) then
          call DCOPY(ndim, fk(idat), mdat, work(isave), mdat)
          isave = isave + 1
        end if
        idat = idat + 1
        go to 10
      end if
      if (ierrmd .eq. 1) then
c
c ##  THE COMPUTATIONS ARE SIMPLE WHEN THE ERROR ARRAY IS DIAGONAL
c ##  NOTE WE NEED TO ACCOUNT FOR 0 STANDARD DEVIATIONS
c
        ip = 0
        do 200 idat = 1, npdat
          if (einv(idat,idat) .gt. wmach(5)) then
            ip = ip + 1
            call DCOPY(ndim, work(idat), mdat, fkhat(ip), mdat)
            call DSCAL(ndim, einv(idat,idat),  fkhat(ip), mdat)
            yhatm(ip) = y(idat)*einv(idat,idat)
          end if
200     continue
        if (ip .ne. p) write(6,*) ' SEVERE ERROR: IPNEP'
      else
c
c ##  ELSE WE NEED TO PERFORM SOME MATRIX MULTIPLICATIONS
c
        call MMMLT(mdat, mdat, mdat, npdat, ndim, p, einv, work, 
     &            fkhat)
        call MMMLT(mdat, mdat, mdat, npdat, 1,    p, einv, y, 
     &            yhatm)
      end if
*
*___end of INVEK___*
*
      end
