      program          MICRON
************************************************************************
*
* MICRON (Multi-instrument Inversion using Constrained RegularizatiON)
* inverts aerosol size distribution data  from a combination of linear
* instruments.
*
* This work was supported in part by the National Science Foundation.
*
* VARIABLES:(type, input/output, array)
*
* nch      :number of channels in instrument i
* npar     :number of parameters passed to instrument i
* inme     :name of instrument i
* parnm    :parameter names
*
* wmach    :machine parameters
*         3- the machine floating-point precision
*         5- smallest positive machine number
*         7- largest positive machine number
*        11- the default output unit number
*
* tol      :calculation tolerances
*         1- minimum relative delta discrepancy before termination
*         3- minimum relative delta smoothing parameter before termination
*         4- minimum nonzero value of a kernel function relative to max
*         5- maximum value of abs(ln(d/ds)) where ds is the diameter where
*           the kernel function just becomes nonzero and d is the
*           approximation to ds
*         6- the accuracy with which we perform integrations
*         7-
*         8- the % padding added to the kernel function intervals
*         9- used to decide if there are intervals in (d1,d2)
*              that are not in the support of the kernel functions.
*        10- used to decide if d1 or d2 lie inside  the support interval
*        11- used to decide if d1 or d2 lie outside the support interval
*
* opchr    :(I, ,A)    do if positive
*         1- specifies smoothing search algorithm
*            1:constrained   GCV                   (CGCV)
*            3:constrained   ridge regression      (CRR)
*            4:discrepency                         (DISC)
*            5:unconstrained GCV                   (UGCV)
*            7:unconstrained ridge regression      (URR)
*         2- constrain the sum of the errors to be positive
*         3- recompute tol(2), the accuracy used to find the target
*            discrepancy
*         4- do TGTNRM calculations
*         5- compute solution characteristics
*         6- compute fkhat, a, do QR (do HESLIN calculations)
*         7- 0:the vector ichpnt and inspnt were not changed for
*                  the present inversion
*         8- compute the minimum error
*         9- invert the error array
*        10- compute the matrix of kernel function values
*        11- 1:invert number distribution (number distribution is
*              smoother)
*            3:invert mass distribution  (mass distribution is smoother)
*        12- norm type
*        13- 0:weak solutions are not ok
*        14- 0:do not call INPSPC
*        15- 0:do not call OUTSPC
*
* msglev   :controls the amount of output written to
*         1- echo file     (.ech)
*         2- log  file     (.log)
*         3- solution file (.out)
*         4- plot file     (.plt)
*         5- screen
*
* nlog     :unit numbers for the corresponding to msglev
*
* ispc     :special counters
*         1- # of smp's chosen by FSPC without succ. subs.
*         2- # of smp's chosen by FSPC based on a bounded interval
*
* nspc     :special integer variables
*         1- positive if the
*
* dspc     :diameter variables
*         1- the lower inversion diameter
*         2- the upper inversion diameter
*         3- the minimum possible inversion diameter
*         4- the maximim possible inversion diameter
*         5- the lower diameter defined by the variance of the
*            kernel functions
*         6- the upper diameter defined by the variance of the
*            kernel functions
*         7- the diameter spacing of the solution (delta)
*
* bspc     :bounds variables
*         1- the constrained solution value at the smallest diameter
*         2- the constrained solution value at the largest  diameter
*         3- the lower integral bound
*         4- the upper integral bound
*
* solspc   :solution properties
*         1- the log-mean diameter
*         2- geometric standard deviation
*         3- penalty term (integral of square 2-derivative for example)
*         4- total number
*         5- total volume
*
* drvx     :x values used in the minmization of R(x) or V(x),
*           see the Numerical Recipies routine BRENT,
*         1-8 u, v, w, x, a, b, d, e
*
* drvy     :function values used in the minimization of R(x) or V(x)
*           see Numerical Recipies routine BRENT
*         1-4 fu, fv, fw, fx
*
* LAST MODIFIED:13 JAN 1988
*
* Logic error: updating of micron.kri should not occur in INPUT.
************************************************************************
      integer           fail,   ibnd,   ichpnt, iconrc, id1,    id2,
     &                  IDSAFE, ierr,   ierrmd, ifail,  igap,   inform,
     &                  inspnt, iinst,  ioerr,  ipower, ismth,  ispc,
     &                  istate, isum,   isumfn, iter,   iterqp, itersv,
     &                  itnsmp, iwork,  jcha,   kerr,   kfiles, kx,
     &                  kxr,    ldim,   lenw,   loc1,   loctie, maderr,
     &                  madcns, maxfev, maxkey
      integer           maxtie, maxtok, mbrent, mclin,  mctotl, mdat,
     &                  mdim,   merr,   mfail,  mfiles, mfun,   minf,
     &                  minfev, minsiv, minst,  minstt, mintg,  miter,
     &                  mlabel, mlog,   mltdim, mopchr, mpar,   mpart,
     &                  mrecl,  msglev, msgscr, msmth,  msuffx, msum,
     &                  msumt,  mtfpar, mtfprt, mtol,   mwmach, mword
      integer           ncha,   nch,    nclin,  ndigit, ndim,   nextra,
     &                  nfail,  ninp,   ninst,  njob,   nlf,    nlog,
     &                  npar,   npdat,  nplot,  nranka, nrankr, nsamp,
     &                  nscr,   nseed,  nterm,  ntfpar, ntoken, opchr,
     &                  p,      skip
c
      double precision  a,      ainp,   b,      bl,     bspc,   bu,
     &                  c,      clamda, cvec,   decdiv, diam,   dpar,
     &                  drvhst, drvx,   drvy,   dsd,    dspc,   dvalue,
     &                  e,      einv,   enorm,  errmin, errpct, f,
     &                  f2,     fk,     fkhat,  obj,    objvec, par,
     &                  r,      rmax,   rmin
      double precision  solspc, sp,     spdsda, sphst,  spinfb, sqemin,
     &                  ssq1,   tfpar,  tnrmsq, tol,    tpar,   ttfpar,
     &                  vlfam,  vlfbm,  volx,   work,   wmach,  y,
     &                  yhatm,  yhatr,  yinf,   ypctdf, yr
c
      character*132     line
      character*80      header
      character*16      fsuffx, key,    keys,   smooth, tie,    ties,
     &                  tokens
      character*16      inpnm
      character*16      cfiles, drfile, fname,  label,  smfile, tdstnm,
     &                  tfnme
      character*8       inme,   parnm,  tfpnm
      logical           bnd,    done,   make,   newlin, sorted, weak
c
c ##PARAMETER: DO NOT ADJUST
c ##Parameter: SOME EXTRA EDITING REQUIRED, SEE MANUAL
c ##parameter: ADJUSTABLE, SEE MANUAL
c
      parameter        (ldim   = 40)
      parameter        (madcns = 5)
      parameter        (maderr = 10)
      parameter        (maxfev = 400)
      parameter        (maxkey = 21)
      parameter        (maxtie = 2)
      parameter        (mbrent = 100)
      parameter        (mdat   = 300)      ! changed by S.P. 1/7/94
      parameter        (mdim   = 201)
      parameter        (mclin  = mdat  + madcns)
      parameter        (mctotl = mclin + mdim)
      parameter        (merr   = mdat  + maderr)
c
c ##COMMENT OUT THE PARAMETER STATEMENT THAT RESULTS IN THE
c ##SMALLER VALUE OF lenw
c
      parameter        (lenw   = 2*mdim*mdim + 10*mdim + 6*mclin)
*     parameter        (lenw   = (3*merr*merr + 7*merr)/2 + mdat + 1)
*     parameter        (lenw   = maxfev)
      parameter        (mfail  = 20)
      parameter        (mfiles = 10000)
      parameter        (mfun   = 5)
      parameter        (minf   = mdim + mdim)
      parameter        (minfev = 50)
      parameter        (minsiv = 10)
      Parameter        (minst  = 10)
      parameter        (mintg  = 500)
      parameter        (miter  = 20)
      parameter        (mlabel = 10)
      PARAMETER        (mlog   = 5)
      parameter        (mltdim = 1)
      PARAMETER        (mwmach = 15)
      PARAMETER        (mopchr = 15)
      Parameter        (mpar   = 6)
      parameter        (mrecl  = 100)
      parameter        (msmth  = 7)
      parameter        (msuffx = 5)
      parameter        (msum   = 10)
      parameter        (mtfpar = 5)
      PARAMETER        (mtol   = 11)
      parameter        (maxtok = 30)
      parameter        (mword  = 8)
      parameter        (decdiv = 50.d00)
      parameter        (sorted = .true.)
c
      common /cmfun/   msumt,  mtfprt
      common /cmipow/  ipower
      common /cmkerf/  iinst,  jcha,   ncha
      common /cmp/     p
      common /cmpar/   par(minst, mpar)
      common /cmparm/  minstt, mpart
      common /cmrmin/  rmin
      common /cmsqem/  sqemin
      common /cmtfpa/  isum,   isumfn(msum), tfpar(msum, mtfpar)
      common /solmch/  wmach(mwmach)
      common /cmrmax/  rmax
c
      dimension        fail(mfail),  ichpnt(mdat), iconrc(mclin),
     &                 inspnt(mdat), ispc(2), istate(mctotl),
     &                 iwork(mdim), kx(mdim), kxr(mdim), msglev(mlog),
     &                 nch(minst),nlog(mlog), npar(minst), ntfpar(mfun),
     &                 opchr(mopchr)
c
      dimension        a(mdim,mdim), ainp(mdim,mdim), b(mdat),
     &                 bl(mctotl), bspc(4), bu(mctotl), c(mclin,mdim),
     &                 clamda(mctotl), cvec(mdim), diam(mdim),
     &                 dpar(minst,mpar),drvhst(0:miter), drvx(8),
     &                 drvy(4), dsd(3), dspc(7), dvalue(2),
     &                 e(mdat,merr), einv(mdat,mdat), enorm(mdat),
     &                 errmin(minsiv)
      dimension        errpct(minsiv), f(mdim), f2(mdim), fk(mdat,mdim),
     &                 fkhat(mdat,mdim), objvec(0:miter), r(mdat,mdim),
     &                 sp(3), solspc(5), sphst(0:miter+1), tol(mtol),
     &                 tpar(minsiv,mpar), ttfpar(mfun,mtfpar)
      dimension        vlfam(minf), vlfbm(minf), volx(mdim),
     &                 work(lenw), y(mdat), yhatm(mdat), yhatr(mdat),
     &                 yinf(minf), ypctdf(mdat), yr(mdat)
c
      dimension        header(2), cfiles(8), fsuffx(msuffx),
     &                 inme(minst), keys(maxkey), label(mlabel),
     &                 parnm(minst,mpar), smooth(msmth), tfnme(mfun),
     &                 tfpnm(mfun,mtfpar), ties(maxtie), tokens(maxtok)
c
c ##INITIALIZE VARIABLES
c
      data             msglev
     &                /  20,     40,     100,      0,     30/
      data opchr      /-1,  0,  1,  1,  1,  1,  1,  1,  1,  1,
     &                  1,  3,  1,  0,  0/
      data             dspc(3), dspc(4),  itersv,   ierrmd, make
     &                /1.d-3,   11.d+0,      1,        1,   .true./
      data             nplot, nterm,  nlog(5), skip,      tdstnm
     &                / 100,    5,    6,       0,      'test_dis'/
      data             nsamp, nseed / 0, 111 /
      data tol        /1.d-4,  0.25d+0, 1.d-6,   1.d-3,   1.d-3,  1.d-2,
     &                 1.d-2,  0.d+0,   0.15d+0, 0.15d+0, 0.15d+0/
      data cfiles
     &/'micron.dat      ', 'micron.ech      ', 'micron.jbs      ',
     & 'micron.kri      ', 'micron.cpy      ', 'micron.lbl      ',
     & 'micron.sim      ', 'micron.drp      '/
      data              fsuffx
     &/'.ECH            ', '.LOG           ',  '.OUT           ',
     & '.PLT            ', '.SCR           '/
      data              keys
     &/'CALCULATE       ', 'COPY:       ECHO', 'CREATE          ',
     & 'DROP            ', 'ECHO            ', 'INFORMATION     ',
     & 'INSTRUMENT      ', 'INVERT          ', 'INPSPC          ',
     & 'LEAST           ', 'MAKE:     CREATE', 'MESSAGE         ',
     & 'OUTSPC          ', 'PLOT:  CALCULATE', 'POINTS          ',
     & 'SIMULATE        ', 'RECOVER         ', 'SMOOTHING       ',
     & 'TOLERANCE       ', 'UPDATE          ', 'WEAK            '/
      data              smooth
     &/'CGCV            ', 'CRIDGE:      CRR', 'CRR             ',
     & 'DISCREPANCY     ', 'UGCV            ', 'URIDGE:      URR',
     & 'URR             '/
      data              ties
     &/'NO              ', 'YES             '/
*-----------------------------------------------------------------------

CRAY  call DROPFILE(0)

      msgscr    = msglev(5)
      nscr      = nlog(5)
      minstt    = minst
      mpart     = mpar
      msumt     = msum
      mtfprt    = mtfpar
      header(1) = ' '
      header(2) = ' '
      wmach(7)  = -1.d+0
      itnsmp    = miter
c
c #READ THE DATA FILE micron.dat FOR A DESCRIPTION OF EACH
c #INSTRUMENT (e.g. name, number of channels, number of parameters)
c
      call MICDAT(fail,   maxtok, mdat,   mfail,  minst,  mpar,
     &            mrecl,  msglev, msmth,  msuffx, mword,  nch,
     &            nlog,   npar,   nterm,  opchr,  dpar,   dspc,
     &            rmax,   rmin,   par,    tol,    wmach,  cfiles,
     &            fsuffx, inme,   parnm,  smooth, tokens)
      call MICCHK(fail,   mfail,  minst,  mpar,   mrecl,  mword,
     &            nlog,   nterm,  dspc,   wmach,  parnm)
      call MICECH(fail,   mfail,  minst,  mpar,   msglev, nch,
     &            nlog,   npar,   nterm,  opchr,  dspc,   tol,
     &            wmach,  cfiles, inme,   parnm,  smooth)
c
c #READ THE OUTPUT LABELS (E.G. dN/dLOG(diam) )
c
      call LBLRD(mlabel, msglev, nlog,   cfiles, label)
c
c #SAY HELLO TO THE USER
c
      call HELLO (fail,   mfail,  msglev, nlog)
c
c #OPEN THE FILE OF JOBS
c
      if (fail(mfail) .le. 0) then
        call FUNIT(njob)
        open(njob, err = 1000, iostat = ioerr, file = cfiles(3),
     &       status = 'old')
        rewind(njob)
1000    if (ioerr .ne. 0)
     &    call ERROR(58, 6, 100, ioerr, 0.d00, cfiles(3), fail(mfail))
      end if
c
c ##BEGIN READING micron.jbs
c
      ntoken = maxtok
      call RDSTR (njob,   ierr,   line,   ntoken, tokens)
10    if ((fail(mfail) .le. 0) .and. (ntoken .gt. 0)) then
        newlin = .true.
        if (skip .le. 0) then
          ndigit = 1
          call TOKNUM(ndigit, ntoken-1, tokens(2), dvalue,1)
          key    = tokens(1)
          call OPLOOK(maxkey, keys, 1, sorted, key, loc1)
c
c ##CHECK IF WE ARE TO SKIP READING INPUT
c
          if (key .eq. 'SKIP') then
            skip = IDSAFE(dvalue(1))
            if (skip .le. 0) skip = 1
          end if
        else
          key = 'SKIP'
        end if
        if (((key .eq. 'SKIP') .or. (key .eq. 'INVERT') .or.
     &      (key .eq. 'INFORMATION') .or. (key .eq. 'SIMULATE'))
     &      .and. (make)) then
          make = .false.
          call KRICHK(fail,   iinst,  maxfev, maxkey, maxtok, mfail,
     &                minfev, minst,  mintg,  mpar,   mrecl,  msgscr,
     &                nch,    npar,   nscr,   decdiv, dpar,   dspc,
     &                par,    tol,    wmach,  work,   cfiles, inme,
     &                parnm,  tokens)
          if (fail(mfail) .gt. 0) go to 10
        end if
        if ( (key .eq. 'SKIP') .or. (key .eq. 'INVERT') .or.
     &       (key .eq. 'INFORMATION') ) then
          if (key .eq. 'SKIP') then
            opchr(10) = 0
            opchr(7)  = 0
          else
            fname = tokens(2)
            call OPUPPR(fname, 'L')
            call MAKEFN(fname, '.inp', inpnm)
          end if
          kfiles = kfiles + 1
c
c ##OPEN ALL THE DATA FILES MICRON WILL USE FOR THIS INVERSION
c
          call FILES (fail,   kfiles, mfail,  msglev, nextra, ninp,
     &                nlog,   skip,   fname,  inpnm,  key)
c
c ##READ THE DATA FILE CONTAINING THE INSTRUMENT DATA AND ANY
c ##ADDITIONAL INFORMATION ABOUT THE INVERSION
c
          call INPUT (fail,   ichpnt, ierrmd, iinst,  inspnt, itnsmp,
     &                iwork,  kerr,   maxfev, maxtok, mdat,   mdim,
     &                merr,   mfail,  minfev, minsiv, minst,  mintg,
     &                miter,  mpar,   mrecl,  msglev, nch,    ndim,
     &                ninp,   ninst,  nlog,   npar,   npdat,  nsamp,
     &                opchr,  skip,   bl,     bspc,   bu,     decdiv,
     &                dpar,   dspc,   e,      par,    tol,    tpar,
     &                sp(1),  wmach,  work,   y,      cfiles, header,
     &                inme,   inpnm,  parnm,  tokens)
c
c ##DETERMINE THE INVERSION INTERVAL
c
          call BOUND (fail,   ibnd,   ichpnt, id1,    id2,    igap,
     &                inspnt, mdat,   mdim,   mfail,  mrecl,  msglev,
     &                nlog,   opchr,  dspc,   tol,    wmach,  work,
     &                cfiles)
c
c ##DETERMINE A REASONABLE VALUE FOR THE DIMENSION OF THE SOLUTION
c ##VECTOR
c
          call CMPDM (fail,   ichpnt, inspnt, ldim,   mdat,   mdim,
     &                mfail,  mltdim, mrecl,  msglev, ndim,   nlog,
     &                npdat,  dspc,   wmach,  work,   cfiles)
c
c ##PERFORM ANY COMMANDS THE USER HAS WRITTEN TO MODIFY THE INPUT
c
          if (opchr(14) .ne. 0)
     &    call INPSPC(fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,
     &                kfiles, mctotl, mdat,   mdim,   merr,   mfail,
     &                minsiv, minst,  mpar,   msglev, nch,    ndim,
     &                ninp,   ninst,  nlog,   npar,   npdat,  opchr,
     &                skip,   bl,     bspc,   bu,     dspc,   e,
     &                par,    sp(1),  tpar,   wmach,  work,   y,
     &                header, fname,  inme)
c
c ##ECHO THE INPUT DATA
c
          call ECHO  (fail,   ichpnt, ierrmd, inspnt, itnsmp, kerr,
     &                mdat,   mdim,   minsiv, minst,  miter,  msglev,
     &                ndim,   ninst,  nlog,   npar,   npdat,  opchr,
     &                bl,     bspc,   bu,     dspc,   e,      par,
     &                sp(1),  tpar,   wmach,  y,      header, inme,
     &                parnm,  smooth)
c
c ##CALCULATE THE MATRIX WHICH APPROXIMATES THE KERNEL FUNCTION
c ##RESPONSE
c
          call CWK   (fail,   ichpnt, inspnt, mdat,   mdim,   mfail,
     &                minsiv, minst,  mintg,  mrecl,  msglev, nch,
     &                ndim,   nlog,   npar,   npdat,  opchr,  dspc,
     &                diam,   fk,     par,    tol,    tpar,   wmach,
     &                work,   cfiles, inme)
c
c ##COMPUTE THE INVERSE OF THE ERROR TRANSFORMATION
c
          call SIGTXF(fail,   ierrmd, iwork,  kerr,   mdat,   merr,
     &                msglev, ndim,   npdat,  nlog,   opchr,  p,
     &                e,      einv,   enorm,  fk,     wmach,  work,
     &                y)
          if (key .eq. 'INFORMATION') then
            call INVEK (fail,   ichpnt, ierrmd, mdat,   msglev, ndim,
     &                  nlog,   npdat,  p,      einv,   fk,     fkhat,
     &                  wmach,  work,   y,      yhatm)
            call LOADR (0,      mdim,   ndim,   a,      diam)
            call VMAT  (fail,   kx,     mdat,   mdim,   msglev, ndim,
     &                  nlog,   p,      fkhat,  a,      ainp,   r,
     &                  cvec,   wmach,  work)
            call SAMPLE(fail,   kx,     mdim,   ndim,   nseed,  nsamp,
     &                  a,      diam,   f,      work)
            call CVVEC (fail,   kx,     mdat,   mdim,   ndim,   nlf,
     &                  p,      a,      ainp,   cvec,   diam,   f,
     &                  r,      vlfam,  vlfbm,  work,   yinf)
            call VREP  (fail,   msglev, ndim,   nlf,    nlog,   vlfam,
     &                  vlfbm,  yinf)
          else
c
c ##COMPUTE THE QUADRATIC AND LINEAR TERM IN THE QUADRATIC PROGRAM
c
          call HESLIN(fail,   ichpnt, ierrmd, iwork,  kx,     kxr,
     &                mdat,   mdim,   msglev, ndim,   nlog,   npdat,
     &                nrankr, opchr,  p,      a,      ainp,   b,
     &                cvec,   einv,   fk,     fkhat,  r,      ssq1,
     &                wmach,  work,   y,      yhatm)
c
c ##DETERMINE THE CONSTRAINTS THAT WILL BE PLACED ON THE QUADRATIC
c ##PROGRAM
c
          call FNDCNS(fail,   ichpnt, iconrc, id1,    id2,    ierrmd,
     &                inspnt, mclin,  mctotl, mdat,   mdim,   mfail,
     &                mrecl,  msglev, nclin,  ndim,   nlog,   npdat,
     &                opchr,  p,      bl,     bspc,   bu,     c,
     &                dspc,   e,      enorm,  fk,     wmach,  y,
     &                work,   cfiles)
c
c ##FIND LEAST SQUARES SOLUTION IF DESIRED
c
          call MINERR(fail,   inform, istate, iterqp, iwork,  kx,
     &                lenw,   mclin,  mdim,   msglev, nclin,  ndim,
     &                nlog,   nrankr, opchr,  ainp,   b,      bl,
     &                bu,     c,      clamda, cvec,   f,      obj,
     &                ssq1,   work)
c
c ##COMPUTE THE STANDARD DEVIATION OF THE ERROR
c
          call TGTNRM(fail,   mbrent, mintg,  msglev, nlog,   opchr,
     &                p,      obj,    sqemin, tnrmsq, tol,    wmach,
     &                yhatm)
c
c ##COMPUTE THE INITIAL SMOOTHING PARAMETER
c
          call FRSTSP(fail,   mdat,   mdim,   msglev, ndim,   nlog,
     &                opchr,  p,      drvhst, dspc,   fkhat,  sp(1),
     &                tnrmsq, tol,    wmach,  work,   yhatm)
c
c ##GUESS INITIAL SOLUTION
c
          ismth = iabs(opchr(1))
          if (itnsmp .gt. 0) then
            call LSOPTN('cold start')
            call LSOPTN('problem type = qp4')
            call DLOAD(ndim, 1.d+0, f, 1)
            weak = .false.
            done = .false.
          else
            itersv = 0
            done   = .true.
            call DCOPY(ndim, f, 1, f2, 1)
            if (opchr(8) .gt. 0) then
              sphst(itersv)  = 0.d+0
            else
              if ((ismth .eq. 5) .or. (ismth .eq. 7))
     &          drvhst(itersv) = drvhst(1)
              sphst(itersv)  = sp(1)
            end if
            call NRMDIF(fail,   iter,   mdat,   mdim,   msglev, ndim,
     &                  nlog,   opchr,  p,      bl,     drvhst, dspc,
     &                  f,      fkhat,  obj,    objvec, solspc,
     &                  sphst(itersv),  tnrmsq, wmach,  work,   yhatm,
     &                  yhatr,  smooth)
          end if
c
c ##FIND THE BEST SMOOTHING PARAMETER
c
          iter = 0
20        if ((fail(1) .le. 0) .and. (.not. done)) then
            iter = iter + 1
            if ((ismth .eq. 4) .and. (iter .eq. 2))
     &        call LSOPTN('Warm Start')
c
c ###DETERMINE THE SMOOTHING PARAMETER FOR THE NEXT ITERATE
c
            if (ismth .eq. 4) then
              call FSPD  (fail,   inform, iter,   kfiles, mdim,
     &                    msglev, ndim,   nlog,   opchr,  p,
     &                    a,      dsd,    dspc,   sp,     spdsda,
     &                    spinfb, tnrmsq, tol,    wmach,  bnd)
            else if ((ismth .eq. 1) .or. (ismth .eq. 3))then
              call FSPC  (fail,   ispc,   iter,   opchr,  drvhst, drvx,
     &                    drvy,   sp(1),  sphst,  tol,    weak)
            end if
            sphst(iter) = sp(1)
c
c ##PREPARE THE QP INPUTS, AND SOLVE THE RESULTING QP
c
            call CPHSWS(fail,   iter,   kx,     kxr,    mdim,   mdat,
     &                  msglev, ndim,   nlog,   nranka, nrankr, opchr,
     &                  p,      a,      ainp,   dspc,   r,      sp(1),
     &                  wmach,  work)
            call LSSOL (fail,   nranka, ndim,   nclin,  mclin,  mdim,
     &                  c,      bl,     bu,     cvec,   istate, kx,
     &                  f,      ainp,   b,      inform, iterqp, obj,
     &                  clamda, iwork,  mdim,   work,   lenw,   ssq1 )
c
c ##REPORT ON THE RESULTS OF THE QP SOLVER
c
            call QPREP (fail,   inform, istate, iterqp, kx,     msglev,
     &                  nclin,  ndim,   nlog,   clamda, sp(1),  weak)
c
c ##COMPUTE THE CONSTRAINED VLAM OR RLAM
c
            if ((ismth .eq. 1) .or. (ismth .eq. 3))
     &        call CMINH (fail,   ispc,   istate, iter,   itnsmp, iwork,
     &                    mdat,   mdim,   msglev, nclin,  ndim,   nlog,
     &                    opchr,  p,      bl,     bu,     drvhst, dspc,
     &                    fkhat,  sp(1),  sphst,  tnrmsq, tol,    wmach,
     &                    work,   yhatm)
c
c ###COMPUTE THE FIDELITY OF THE SOLUTION TO THE DATA
c
            call NRMDIF(fail,   iter,   mdat,   mdim,   msglev, ndim,
     &                  nlog,   opchr,  p,      bl,     drvhst, dspc,
     &                  f,      fkhat,  obj,    objvec, solspc, sp(1),
     &                  tnrmsq, wmach,  work,   yhatm,  yhatr,  smooth)
c
c ###DECIDE IF THE SEARCH FOR THE "CORRECT" SMOOTHING PARAMETER SHOULD
c ###BE STOPPED
c
            if (ismth .eq. 4) then
              call DONED (fail,   inform, iter,   itersv, itnsmp,
     &                    mdim,   miter,  msglev, ndim,   nlog,
     &                    drvhst, dsd,    f,      f2,     sp,
     &                    spdsda, sphst,  spinfb, sqemin, tnrmsq,
     &                    tol,    wmach,  done,   bnd)
            else if ((ismth .eq. 1) .or. (ismth .eq. 3)) then
              call DONEC (fail,   ispc,   iter,   itersv, itnsmp,
     &                    msglev, ndim,   nlog,   opchr,  drvhst,
     &                    drvx,   drvy,   f,      f2,     sphst,
     &                    tol,    done,   weak)
            else
              itersv = 1
              done   = .true.
              call DCOPY (ndim,   f,      1,      f2,     1)
            end if
            go to 20
          end if
c
c ##COMPUTE PROPERTIES OF THE SOLUTION THAT WILL BE WRITTEN TO THE
c ##OUTPUT FILE
c
          call PREOUT(fail,   ichpnt, inspnt, istate, itersv, mdat,
     &                mdim,   mfail,  miter,  mrecl,  msglev, ndim,
     &                nlog,   npdat,  opchr,  p,      diam,   drvhst,
     &                dspc,   enorm,  f2,     fk,     fkhat,  objvec,
     &                solspc, sphst,  volx,   wmach,  work,   y,
     &                yhatm,  yhatr,  ypctdf, yr,     cfiles)
c
c ##REPORT THE SOLUTION AND ITS PROPERTIES
c
          call OUTPUT(fail,   ichpnt, iconrc, inspnt, itersv, mctotl,
     &                mdat,   mdim,   merr,   mfail,  minsiv, minst,
     &                miter,  mpar,   msglev, nch,    ndim,   ninst,
     &                nlog,   npar,   npdat,  opchr,  bl,     bspc,
     &                bu,     diam,   drvhst, dspc,   einv,   f2,
     &                objvec, par,    solspc, sphst,  tol,    tpar,
     &                volx,   wmach,  y,      ypctdf, yr,     header,
     &                inme,   parnm,  smooth)
c *
c *-end key .ne. information
c *
          end if
c
c ##EXECUTE ANY COMMANDS THE USERS HAS PLACED IN OUTSPC
c
          if (opchr(15) .ne. 0)
     &    call OUTSPC(fail,   ichpnt, iconrc, inspnt, itersv, kerr,
     &                kfiles, mctotl, mdat,   mdim,   merr,   mfail,
     &                minsiv, minst,  mintg,  miter,  mpar,   msglev,
     &                nch,    ndim,   ninst,  nlog,   npar,   npdat,
     &                opchr,  p,      bl,     bspc,   bu,     diam,
     &                drvhst, dspc,   e,      einv,   enorm,  f,
     &                f2,     fk,     objvec, par,    solspc, sphst,
     &                tol,    tpar,   volx,   wmach,  work,   y,
     &                yhatm,  yhatr,  yinf,   ypctdf, yr,     fname,
     &                header, inme,   inpnm,  parnm,  tdstnm, tfnme)
c
c ##CHECK IF WE'VE HAD TOO MANY INVERSION FAILURES
c
          if (fail(1) .gt. 0) then
            call UINDEX(fail(mfail - 1))
            nfail    = nfail + 1
            if (nfail .ge. nterm) call ERROR(27, nscr, msgscr,
     &      nfail, 0.d00, ' ', fail(mfail))
          else
            fail(mfail - 1) = -abs(fail(mfail - 1))
            nfail    = 0
          end if
c
c ##RESET THE FAIL VECTOR
c
        do 100 ifail = 1, mfail - 2
          fail(ifail) = -abs(fail(ifail))
100     continue
c
c ##OR ECHO THE KERNEL INFORMATION FILE
c
        else if (key .eq. 'INSTRUMENT') then
          call RDINST(line,   inme,   npar,   par,    parnm,  minst,
     &                ntoken, tokens, maxtok, njob,   iinst,  ierr)
          newlin = .false.
          if (ierr .eq. 1)
     &      call ERROR(17, nscr, msgscr, 0, 0.d+0, line, fail(mfail))
        else if (key .eq. 'POINTS') then
          if (ndigit .eq. 1) then
            nplot = abs(IDSAFE(dvalue(1)))
            if (nplot .le. 1) nplot = 2
          else
            call WARNG(2, nscr, msgscr, 0, 0.d+0, ' ')
          end if
        else if (key .eq. 'SMOOTHING') then
          call FNDKEY(smooth, msmth, tokens, ntoken, tie,    loctie,
     &                sorted)
          if (loctie .gt. 0) then
            opchr(1) = loctie
          else
            call ERROR(74, nscr, msgscr, 0, 0.d+0, tie, fail(mfail))
          end if
        else if (key .eq. 'OUTSPC') then
          call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie,
     &                sorted)
          opchr(15) = 0
          if (tie .eq. 'YES') opchr(15) = 1
        else if (key .eq. 'INPSPC') then
          call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie,
     &                sorted)
          opchr(14) = 0
          if (tie .eq. 'YES') opchr(14) = 1
        else if (key .eq. 'WEAK') then
          call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie,
     &                sorted)
          if (tie .eq. 'NO') then
            opchr(13) = 0
          else
            opchr(13) = 1
          end if
        else if (key .eq. 'MESSAGE') then
          call FNDKEY(fsuffx, msuffx, tokens, ntoken, tie,    loctie,
     &                .false.)
          if (loctie .gt. 0) then
            msglev(loctie) = IDSAFE(dvalue(1))
          else
            call WARNG(6, nscr, msgscr, 0, 0.d+0, tie)
          end if
        else if (key .eq. 'LEAST') then
          call FNDKEY(ties, maxtie, tokens, ntoken, tie,    loctie,
     &                sorted)
          if (tie .eq. 'NO') then
            opchr(8) = 0
          else
            opchr(8) = 1
          end if
        else if (key .eq. 'ECHO') then
          call KERECH(fail,   mfail,  mpar,   mrecl,  msgscr, nscr,
     &                work,   cfiles, tokens)
c
c ##OR GENERATE A DATA FILE CONTAINING KERNEL FUNCTION VALUES
c
        else if (key .eq. 'CALCULATE') then
          jcha = IDSAFE(dvalue(1))
          call KERPLT(iinst,  jcha,   msgscr, nch,    ncha,   nplot,
     &                nscr,   dspc,   work,   inme)
c
c ##OR (DROP/RECOVER) THE KERNEL FUNCTION MATRIX
c
        else if ((key .eq. 'DROP') .or. (key .eq. 'RECOVER')) then
          if (ntoken .eq. 1) then
            drfile = cfiles(8)
          else
            call OPUPPR(tokens(2), 'L')
            drfile = tokens(2)
          end if
          call RECDRP(fail,   ichpnt, id1,    id2,    inspnt, mdat,
     &                minsiv, mpar,   msgscr, ndim,   ninst,  npdat,
     &                nscr,   opchr,  diam,   dspc,   fk,     sp(1),
     &                tpar,   drfile, key)
c
c ##OR SIMULATE DATA
c
        else if (key .eq. 'SIMULATE') then
          smfile = tokens(2)
          call OPUPPR(smfile, 'L')
c
c        ##READ micron.sim TO FIND OUT ABOUT THE TEST FUNCTIONS
c        ##BEFORE THE DATA IS GENERATED
c
          call SIMINI(fail,   maxtok, mfail,   mfun,   msgscr, mtfpar,
     &                nscr,   ntfpar, ttfpar,  cfiles, tfnme,  tfpnm,
     &                tokens)
          call SIMDAT(fail,   ichpnt, inspnt, isumfn, maxfev, maxtok,
     &                mdat,   mfail,  mfun,   minfev, minsiv, minst,
     &                mintg,  mpar,   mrecl,  msgscr, msum,   nch,
     &                ninst,  npar,   npdat,  nscr,   nseed,  ntfpar,
     &                decdiv, dpar,   dspc,   e,      errmin, errpct,
     &                par,    tfpar,  tol,    tpar,   ttfpar, wmach,
     &                work,   y,      cfiles, inme,   parnm,  smfile,
     &                tdstnm, tfnme,  tfpnm,  tokens)
c
c ##CREATE THE KERNEL INFORMATION FILE
c
        else if (key .eq. 'CREATE') then
          make = .false.
          call KERINF(fail,   maxfev, maxtok, mfail,  minfev, minst,
     &                mintg,  mpar,   mrecl,  msgscr, nch,    npar,
     &                nscr,   decdiv, dspc,   par,    tol,    wmach,
     &                work,   cfiles, inme,   parnm,  tokens)
c
c ##OR UPDATE THE KERNEL INFORMATION FILE
c
        else if (key .eq. 'UPDATE') then
          make = .false.
          call UPDATE(fail,   maxfev, mfail,  minfev, minst,  mintg,
     &                mpar,   mrecl,  msgscr, nch,    npar,   nscr,
     &                decdiv, dspc,   par,    tol,    wmach,  work,
     &                cfiles, inme,   parnm)
c
c ##ELSE THE KEYWORD WAS NOT RECOGNIZED
c
        else
          call ERROR(3, nscr, msgscr, 0, 0.d+0, key, fail(mfail))
        end if
        if (key .eq. 'SKIP') skip = skip - 1
        if ((skip .le. 0) .and. (newlin)) then
          ntoken = maxtok
          call RDSTR (njob,   ierr,   line,   ntoken, tokens)
        end if
        go to 10
      end if
      if (fail(mfail) .gt. 0) then
        call SCRNM(19, nscr, msgscr, 0, 0.d+0, ' ')
      else
        call SCRNM(18, nscr, msgscr, 0, 0.d+0, ' ')
      end if
      stop
*
*___end of MICRON___*
*
      end
      subroutine   FSPC(fail,   ispc,   iter,   opchr,  drvhst, drvx,
     &                  drvy,   smp,    sphst,  tol,    weak)
      integer           fail,   ispc,   iter,   opchr
      double precision  drvhst, drvx,   drvy,   smp,    sphst,  tol
      logical           weak
      dimension         fail(*), ispc(*), opchr(*), drvhst(0:*),
     &                  drvx(*), drvy(*), sphst(0:*), tol(*)
************************************************************************
* FSPC chooses the next smoothing parameter to use to minimize
* V(lambda) or R(lambda).
* We make sure lambda is large enough to prevent a weak minimum.
*
* LAST MODIFIED:12 DEC 1988
************************************************************************
      double precision  sclini, sclmul, smpsmn, smpwmx, spscl,  zero
      parameter        (sclini = 5.0d+0)
      parameter        (sclmul = 1.5d+0)
      parameter        (zero   = 0.0d+0)
      save              smpsmn, smpwmx, spscl
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      if (iter    .eq. 1) then
        smpwmx  = -1.d+0
        smpsmn  = -1.d+0
        ispc(1) = 0
        ispc(2) = 0
        return
      else if (weak) then
        if (smp .gt. smpwmx) smpwmx = smp
      else
        if ((smpsmn .lt. zero) .or. (smp .lt. smpsmn)) smpsmn = smp
      end if
      if (ispc(1) .eq. 0) then
        smp = sphst(iter)
      else if (ispc(1) .eq. 1) then
        if (weak .and. (opchr(13) .eq. 0) .and. (smpsmn .gt. zero)) then
          smp = dsqrt( smpwmx * smpsmn )
        else
          spscl = sclini
          smp   = dexp(drvx(2)) * spscl
        end if
      else if (ispc(2) .eq. 0) then
        spscl = spscl * sclmul
        if (drvy(2) .lt. drvy(3)) then
          smp = dexp(drvx(2)) / spscl
        else
          smp = dexp(drvx(3)) * spscl
        end if
      else if (ispc(2) .ge. 1)  then
        call DBMCOR(drvx,   drvy,   tol(2))
        smp = dexp(drvx(1))
      end if
c
c ##AVOID WEAK MINIMUMS
c
      if ((opchr(13) .eq. 0) .and. (smp .le. smpwmx)) then
        if (smpsmn .gt. zero) then
          smp = dsqrt(smpsmn * smpwmx)
        else
          smp = smp * spscl
        end if
        drvx(1) = dlog(smp)
      end if
      return
*
*___End of FSPC___*
*
      end
      subroutine RECDRP(fail,   ichpnt, id1,    id2,    inspnt, mdat,
     &                  minsiv, mpar,   msgscr, ndim,   ninst,  npdat,
     &                  nscr,   opchr,  diam,   dspc,   fk,     smp,
     &                  tpar,   drfile, job)
      integer           fail,   ichpnt, id1,    id2,    inspnt, mdat,
     &                  minsiv, mpar,   msgscr, ndim,   ninst,  npdat,
     &                  nscr,   opchr
      double precision  diam,   dspc,   fk,     smp,    tpar
      character*(*)     drfile, job
      dimension         fail(*), ichpnt(*), inspnt(*),  opchr(*),
     &                  diam(*), dspc(*),   fk(mdat,*), tpar(*)
************************************************************************
* RECDRP drops or recovers the kernel function matrix fk, par, and
* the pointer matrices in the unformatted data file micron.r.
*
* LAST MODIFIED: 6 FEB 1989
************************************************************************
      integer           idat,   idim,   ioerr,  ipar,   joerr,  mdpar,
     &                  ndrp
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
      call FUNIT(ndrp)
      mdpar = minsiv * mpar
      if (job .eq. 'DROP') then
        call SCRNM(49, nscr, msgscr, 0, 0.d+0, drfile)
        open(ndrp, file=drfile, form = 'unformatted',
     &       status='unknown')
        write(ndrp) id1, id2, ndim,  ninst,  npdat,  smp
        do 100 idat = 1, npdat
          write(ndrp) ichpnt(idat), inspnt(idat)
          write(ndrp) (fk(idat,idim),idim = 1, ndim)
100     continue
        write(ndrp) (diam(idim), idim = 1, ndim)
        write(ndrp) dspc(5), dspc(6)
        write(ndrp) (tpar(ipar), ipar = 1, mdpar)
      else
c
c ##ELSE RECOVER
c
        call SCRNM(49, nscr, msgscr, 1, 0.d+0, drfile)
        open(ndrp, file=drfile, status='old', form = 'unformatted',
     &       iostat = ioerr, err = 1000)
        read(ndrp, iostat = joerr, err = 1000) id1, id2, ndim,  ninst,
     &                                         npdat, smp
        do 200 idat = 1, npdat
          read(ndrp, iostat = joerr, err = 1000)
     &         ichpnt(idat), inspnt(idat)
          read(ndrp, iostat = joerr, err = 1000)
     &         (fk(idat,idim), idim = 1, ndim)
200     continue
        read(ndrp, iostat = joerr, err = 1000)
     &       (diam(idim), idim = 1, ndim)
        read(ndrp, iostat = joerr, err = 1000)
     &       dspc(5), dspc(6)
        read(ndrp, iostat = joerr, err = 1000)
     &       (tpar(ipar), ipar = 1, mdpar)
        dspc(1) = diam(1)
        dspc(2) = diam(ndim)
      end if
c
c ##CHECK FOR ERRORS
c
1000  if (ioerr .ne. 0) then
        call ERROR(97, nscr, msgscr, ioerr, 0.d+0, 'open', fail(2))
      else
        close(ndrp)
        if (joerr .ne. 0)
     &  call ERROR(97, nscr, msgscr, joerr, 0.d+0, 'read', fail(2))
      end if
      return
*
*___end of RECDRP___*
*
      end
      subroutine FILES (fail,   kfiles, mfail,  msglev, nextra, ninp,
     &                  nlog,   skip,   fname,  inpnm,  key)
      integer           fail,   kfiles, mfail,  msglev, nextra, ninp,
     &                  nlog,   skip
      character*(*)     fname,  inpnm,  key
      dimension         fail(*), msglev(*), nlog(*)
************************************************************************
*
* FILES opens  all the output files used by MICRON.  Each file name
* has the form "root"."extension".  All the files have the same root
* name which is listed as an entry in the micron.jbs data file.
* The extension names
* are (.inp .ech .out .log .plt).  The file is not opened if the message
* level of that file is less than or equal to 0;  the input file is not
* opened if skip is positive.
*
* INPUTS:
* msgech  - integer parameter that controls the amount of output
*           to the ".ech" file
* skip    - a integer variable that is positive if the input file should
*           not be opened
*
* OUTPUTS:
* nech    - the unit number for the ".ech" file
*
* LAST REVISION:13 JAN 1988
************************************************************************
      integer           iclose, ioerr,  iopen,  iopens, iscr,   itemp,
     &                  joerr,  nopen
      double precision  xfile
      character*12      filenm
      character*4       suffix
      logical           oplog
      parameter        (nopen = 4)
      dimension         suffix(nopen)
      data suffix      /'.ech', '.log', '.out', '.plt'/
*-----------------------------------------------------------------------
      if (fail(1) .gt. 0) return
c
c ##THE FOLLOWING STATEMNET HELPS US CLEAN UP IF AN ERROR OCCURRS
c
      call ILOAD(nopen, 0, nlog, 1)
      ninp   = 0
      iopens = 1
c
      itemp  = nlog(4)
      if (key .eq. 'INFORMATION') then
        iscr      = 11
        msglev(4) = 0
      else
        iscr      = 10
      end if
      if (skip .gt. 0) then
        call SCRNM(9, nlog(5), msglev(5), 0, 0.d+0, ' ')
      else
        call SCRNM(iscr, nlog(5), msglev(5), 0, 0.d00, inpnm)
        call FUNIT(ninp)
        open(ninp, iostat=ioerr, err=1000, file=inpnm, status='old')
        rewind(ninp)
      end if
c
      do 100 iopen = 1, nopen
        iopens = iopen
        if (msglev(iopen) .gt. 0) then
          call MAKEFN(fname, suffix(iopen), filenm)
          call FUNIT(nlog(iopen))
          open(nlog(iopen), iostat=joerr, err=1000, file=filenm,
     &         status='unknown')
        end if
100   continue
c
1000  xfile = dble(kfiles + nextra)
      if (joerr .ne. 0) call ERROR(8,  nlog(5), msglev(5),
     &  joerr, xfile, filenm, fail(3))
      if (ioerr .ne. 0)
     &  call ERROR(28, nlog(5), msglev(5), ioerr, xfile, inpnm, fail(3))
      if (fail(3) .gt. 0) then
        call UINDEX(fail(1))
c
c ##IF AN ERROR OCCURRED, CLOSE ALL THE FILES WHICH HAVE BEEN
c ##OPENED
c
        if (ninp .ne. 0) inquire(ninp, opened = oplog)
        if (oplog) close(ninp)
        do 200 iclose = 1, iopens - 1
          if (nlog(iclose) .ne. 0) inquire(nlog(iclose), opened = oplog)
          if (oplog) close(nlog(iclose))
200     continue
        call ILOAD(nopen, 0, nlog, 1)
      end if
      msglev(4) = itemp
      return
*
*___end of FILES___*
*
      end
      subroutine HELLO (fail,   mfail,  msglev, nlog)
      integer           fail,   mfail,  msglev, nlog
      dimension         fail(*), msglev(*), nlog(*)
************************************************************************
*
* HELLO shows the MICRON logo if micron.dat has been read successfully
* and the message level is large enough.
*
* LAST MODIFIED: 7 DEC 1987
************************************************************************
      integer           itime,  jtime,  ntitle
      character*1       t,d,m,i,c,r,o,n,e,s,a
      character*8       q
      parameter        (t='*',d='-',m='M',c='C',r='R',i='I',o='O',n='N')
      parameter        (e='.',s='''',a=',',q='*-*-*-*-')
c
      if (msglev(5) .le. 0) return
      if (fail(mfail) .le. 0) then
        ntitle = nlog(5)
        if (msglev(5) .lt. 100) then
          write(ntitle, 12)
        else
          write(ntitle,2) a,q,q,q,q,q,q,q,t,a
          write(ntitle,3) a,t,t,a
          write(ntitle,4) t,t
          write(ntitle,4) t,t
          write(ntitle,5) t,a,c,d,c,r,d,r,d,r,a,a,o,d,o,a,t
          write(ntitle,6) t,m,m,m,m,i,c,s,r,r,o,o,n,n,a,n,t
          write(ntitle,7) t,i,c,r,d,r,d,r,s,o,o,n,n,a,n,t
          write(ntitle,8) t,i,c,a,r,r,o,a,a,o,n,n,n,t
          write(ntitle,9) t,m,m,i,i,i,s,c,d,c,r,r,o,d,o,n,n,t
          write(ntitle,4) t,t
          write(ntitle,4) t,t
          write(ntitle,10) t,t
          write(ntitle,11) t,t
          write(ntitle,4) t,t
          write(ntitle,4) t,t
          write(ntitle,3) s,t,t,s
          write(ntitle,2) s,q,q,q,q,q,q,q,t,s
          do 100 itime=1,1000
            do 110 jtime=1,100
110         continue
100       continue
          write(ntitle,1)
        end if
      else
       write(ntitle, 13)
      end if
1     format('  ',//////////////////////////////////)
2     format(t9,a1,7a8,2a1)
3     format(t8,2a1,t67,2a1)
4     format(t8,a1,t68,a1)
5     format(t8,a1,3x,'M',6x,'M',3x,'III',4x,4a1,4x,6a1,5x,5a1,
     &       4x,'N',5x,'N',3x,a1)
6     format(t8,a1,3x,2a1,4x,2a1,3(4x,a1),3x,a1,4x,a1,4x,a1,5x,a1,3x,
     &       3a1,3x,a1,3x,a1)
7     format(t8,a1,3x,'M M  M M',2(4x,a1),8x,6a1,4x,a1,5x,a1,3x,a1,2x,
     &       2a1,1x,a1,3x,a1)
8     format(t8,a1,3x,'M  MM  M',3(4x,a1),3x,a1,4x,a1,4x,2a1,3x,2a1,3x,
     &       a1,4x,2a1,3x,a1)
9     format(t8,a1,3x,a1,6x,a1,3x,3a1,4x,4a1,4x,a1,5x,a1,5x,3a1,5x,a1,
     &       5x,a1,3x,a1)
10    format(t8,a1,t22,'Multi-instrument Inversion using',t68,a1)
11    format(t8,a1,t25,'Constrained RegularizatiON',t68,a1)
12    format(//,' MICRON is working.',//)
13    format(//)
      return
*
*___end of HELLO___*
*
      end
