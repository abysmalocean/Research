/*     SANDIA MATHEMATICAL PROGRAM LIBRARY */
/*     APPLIED MATHEMATICS DIVISION 2613 */
/*     SANDIA LABORATORIES */
/*     ALBUQUERQUE, NEW MEXICO  87185 */
/*     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978 */
/*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    ISSUED BY SANDIA LABORATORIES                     * */
/*  *                   A PRIME CONTRACTOR TO THE                       * */
/*  *                UNITED STATES DEPARTMENT OF ENERGY                 * */
/*  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * * */
/*  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   * */
/*  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      * */
/*  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    * */
/*  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  * */
/*  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      * */
/*  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     * */
/*  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      * */
/*  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          * */
/*  * OWNED RIGHTS.                                                     * */
/*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     * */
/*  * PART IS SAND77-1441.                                              * */
/*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*   WRITTEN BY L. F. SHAMPINE AND M. K. GORDON */

/*   ABSTRACT */

/*   SUBROUTINE  ODE  INTEGRATES A SYSTEM OF  NEQN  FIRST ORDER */
/*   ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM */
/*             DY(I)/DT = F(T,Y(1),Y(2),...,Y(NEQN)) */
/*             Y(I) GIVEN AT  T . */
/*   THE SUBROUTINE INTEGRATES FROM  T  TO  TOUT .  ON RETURN THE */
/*   PARAMETERS IN THE CALL LIST ARE SET FOR CONTINUING THE INTEGRATION. */
/*   THE USER HAS ONLY TO DEFINE A NEW VALUE  TOUT  AND CALL  ODE  AGAIN. */

/*   THE DIFFERENTIAL EQUATIONS ARE ACTUALLY SOLVED BY A SUITE OF CODES */
/*   DE ,  STEP1 , AND  INTRP .  ODE  ALLOCATES VIRTUAL STORAGE IN THE */
/*   ARRAYS  WORK  AND  IWORK  AND CALLS  DE .  DE  IS A SUPERVISOR WHICH */
/*   DIRECTS THE SOLUTION.  IT CALLS ON THE ROUTINES  STEP1  AND  INTRP */
/*   TO ADVANCE THE INTEGRATION AND TO INTERPOLATE AT OUTPUT POINTS. */
/*   STEP1  USES A MODIFIED DIVIDED DIFFERENCE FORM OF THE ADAMS PECE */
/*   FORMULAS AND LOCAL EXTRAPOLATION.  IT ADJUSTS THE ORDER AND STEP */
/*   SIZE TO CONTROL THE LOCAL ERROR PER UNIT STEP IN A GENERALIZED */
/*   SENSE.  NORMALLY EACH CALL TO  STEP1  ADVANCES THE SOLUTION ONE STEP */
/*   IN THE DIRECTION OF  TOUT .  FOR REASONS OF EFFICIENCY  DE */
/*   INTEGRATES BEYOND  TOUT  INTERNALLY, THOUGH NEVER BEYOND */
/*   T+10*(TOUT-T), AND CALLS  INTRP  TO INTERPOLATE THE SOLUTION AT */
/*   TOUT .  AN OPTION IS PROVIDED TO STOP THE INTEGRATION AT  TOUT  BUT */
/*   IT SHOULD BE USED ONLY IF IT IS IMPOSSIBLE TO CONTINUE THE */
/*   INTEGRATION BEYOND  TOUT . */

/*   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, */
/*   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL */
/*   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON. */
/*   FOR USAGE DETAILS, ALSO SEE SLA-73-1060. */

/*   THE PARAMETERS REPRESENT -- */
/*      F -- SUBROUTINE F(T,Y,YP) TO EVALUATE DERIVATIVES YP(I)=DY(I)/DT */
/*      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED */
/*      Y(*) -- SOLUTION VECTOR AT  T */
/*      T -- INDEPENDENT VARIABLE */
/*      TOUT -- POINT AT WHICH SOLUTION IS DESIRED */
/*      RELERR,ABSERR -- RELATIVE AND ABSOLUTE ERROR TOLERANCES FOR LOCAL */
/*           ERROR TEST.  AT EACH STEP THE CODE REQUIRES */
/*             ABS(LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR */
/*           FOR EACH COMPONENT OF THE LOCAL ERROR AND SOLUTION VECTORS */
/*      IFLAG -- INDICATES STATUS OF INTEGRATION */
/*      WORK(*),IWORK(*) -- ARRAYS TO HOLD INFORMATION INTERNAL TO CODE */
/*           WHICH IS NECESSARY FOR SUBSEQUENT CALLS */

/*   FIRST CALL TO ODE -- */

/*   THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR THE ARRAYS */
/*   IN THE CALL LIST, */
/*      Y(NEQN), WORK(100+21*NEQN), IWORK(5), */
/*   DECLARE  F  IN AN EXTERNAL STATEMENT, SUPPLY THE SUBROUTINE */
/*   F(T,Y,YP) TO EVALUATE */
/*      DY(I)/DT = YP(I) = F(T,Y(1),Y(2),...,Y(NEQN)) */
/*   AND INITIALIZE THE PARAMETERS -- */
/*      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED */
/*      Y(*) -- VECTOR OF INITIAL CONDITIONS */
/*      T -- STARTING POINT OF INTEGRATION */
/*      TOUT -- POINT AT WHICH SOLUTION IS DESIRED */
/*      RELERR,ABSERR -- RELATIVE AND ABSOLUTE LOCAL ERROR TOLERANCES */
/*      IFLAG -- +1,-1.  INDICATOR TO INITIALIZE THE CODE.  NORMAL INPUT */
/*           IS +1.  THE USER SHOULD SET IFLAG=-1 ONLY IF IT IS */
/*           IMPOSSIBLE TO CONTINUE THE INTEGRATION BEYOND  TOUT . */
/*   ALL PARAMETERS EXCEPT  F ,  NEQN  AND  TOUT  MAY BE ALTERED BY THE */
/*   CODE ON OUTPUT SO MUST BE VARIABLES IN THE CALLING PROGRAM. */

/*   OUTPUT FROM  ODE  -- */

/*      NEQN -- UNCHANGED */
/*      Y(*) -- SOLUTION AT  T */
/*      T -- LAST POINT REACHED IN INTEGRATION.  NORMAL RETURN HAS */
/*           T = TOUT . */
/*      TOUT -- UNCHANGED */
/*      RELERR,ABSERR -- NORMAL RETURN HAS TOLERANCES UNCHANGED.  IFLAG=3 */
/*           SIGNALS TOLERANCES INCREASED */
/*      IFLAG = 2 -- NORMAL RETURN.  INTEGRATION REACHED  TOUT */
/*            = 3 -- INTEGRATION DID NOT REACH  TOUT  BECAUSE ERROR */
/*                   TOLERANCES TOO SMALL.  RELERR ,  ABSERR  INCREASED */
/*                   APPROPRIATELY FOR CONTINUING */
/*            = 4 -- INTEGRATION DID NOT REACH  TOUT  BECAUSE MORE THAN */
/*                   500 STEPS NEEDED */
/*            = 5 -- INTEGRATION DID NOT REACH  TOUT  BECAUSE EQUATIONS */
/*                   APPEAR TO BE STIFF */
/*            = 6 -- INTEGRATION DID NOT REACH  TOUT  BECAUSE SOLUTION */
/*                   VANISHED MAKING PURE RELATIVE ERROR IMPOSSIBLE. */
/*                   MUST USE NON-ZERO  ABSERR  TO CONTINUE. */
/*            = 7 -- INVALID INPUT PARAMETERS (FATAL ERROR) */
/*           THE VALUE OF  IFLAG  IS RETURNED NEGATIVE WHEN THE INPUT */
/*           VALUE IS NEGATIVE AND THE INTEGRATION DOES NOT REACH  TOUT , */
/*           I.E., -3, -4, -5, -6. */
/*      WORK(*),IWORK(*) -- INFORMATION GENERALLY OF NO INTEREST TO THE */
/*           USER BUT NECESSARY FOR SUBSEQUENT CALLS. */

/*   SUBSEQUENT CALLS TO  ODE -- */

/*   SUBROUTINE  ODE  RETURNS WITH ALL INFORMATION NEEDED TO CONTINUE */
/*   THE INTEGRATION.  IF THE INTEGRATION REACHED  TOUT , THE USER NEED */
/*   ONLY DEFINE A NEW  TOUT  AND CALL AGAIN.  IF THE INTEGRATION DID NOT */
/*   REACH  TOUT  AND THE USER WANTS TO CONTINUE, HE JUST CALLS AGAIN. */
/*   IN THE CASE  IFLAG=6 , THE USER MUST ALSO ALTER THE ERROR CRITERION. */
/*   THE OUTPUT VALUE OF  IFLAG  IS THE APPROPRIATE INPUT VALUE FOR */
/*   SUBSEQUENT CALLS.  THE ONLY SITUATION IN WHICH IT SHOULD BE ALTERED */
/*   IS TO STOP THE INTEGRATION INTERNALLY AT THE NEW  TOUT , I.E., */
/*   CHANGE OUTPUT  IFLAG=2  TO INPUT  IFLAG=-2 .  ERROR TOLERANCES MAY */
/*   BE CHANGED BY THE USER BEFORE CONTINUING.  ALL OTHER PARAMETERS MUST */
/*   REMAIN UNCHANGED. */

int ode(U_fp f, int *neqn, double *y, double *t,double *tout, double *relerr, double *abserr,int *iflag, double *work, int *iwork)
{
    /* Initialized data */
    static int ialpha = 1;
    static int ibeta = 13;
    static int isig = 25;
    static int iv = 38;
    static int iw = 50;
    static int ig = 62;
    static int iphase = 75;
    static int ipsi = 76;
    static int ix = 88;
    static int ih = 89;
    static int ihold = 90;
    static int istart = 91;
    static int itold = 92;
    static int idelsn = 93;

    extern int de(U_fp, int *, double *, double *, double *, double *, double *, int *, 
	double *, double *, double *, double *, 
	double *, double *, double *, double *, 
	double *, double *, double *, double *, logical *,
	double *, double *, double *, double *, logical *, 
	double *, double *, int *, logical *, int *, int *, int *);
    static int ip, iyp, iwt, iyy, iphi;
    static logical nornd, start, phase1;
    static int iypout;

    /* Parameter adjustments */
    --y;
    --work;
    --iwork;

    /* Function Body */
    iyy = 100;
    iwt = iyy + *neqn;
    ip = iwt + *neqn;
    iyp = ip + *neqn;
    iypout = iyp + *neqn;
    iphi = iypout + *neqn;
    if (abs(*iflag) < 2 || abs(*iflag) > 6) 
    {
//	goto L1;
//L1:
	    de_((U_fp)f, neqn, &y[1], t, tout, relerr, abserr, iflag, &work[iyy], 
		&work[iwt], &work[ip], &work[iyp], &work[iypout], &work[iphi], 
		&work[ialpha], &work[ibeta], &work[isig], &work[iv], &work[iw], 
		&work[ig], &phase1, &work[ipsi], &work[ix], &work[ih], &work[ihold], 
		&start, &work[itold], &work[idelsn], &iwork[1], &nornd, &iwork[3], &iwork[4], &iwork[5]);
	    work[istart] = -1.0;
	    if (start) 
  	    {
		work[istart] = 1.0;
	    }	
	    work[iphase] = -1.0;
	    if (phase1) 
	    {
		work[iphase] = 1.0;
	    }
	    iwork[2] = -1;
	    if (nornd) 
	    {
		iwork[2] = 1;
	    }
    }
    start = work[istart] > 0.0;
    phase1 = work[iphase] > 0.0;
    nornd = iwork[2] != -1;
    return 0;
} /* ode_ */

