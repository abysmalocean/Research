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

/*   SUBROUTINE  STEP1  IS NORMALLY USED INDIRECTLY THROUGH SUBROUTINE */
/*   ODE .  BECAUSE  ODE  SUFFICES FOR MOST PROBLEMS AND IS MUCH EASIER */
/*   TO USE, USING IT SHOULD BE CONSIDERED BEFORE USING  STEP1  ALONE. */

/*   SUBROUTINE STEP1 INTEGRATES A SYSTEM OF  NEQN  FIRST ORDER ORDINARY */
/*   DIFFERENTIAL EQUATIONS ONE STEP, NORMALLY FROM X TO X+H, USING A */
/*   MODIFIED DIVIDED DIFFERENCE FORM OF THE ADAMS PECE FORMULAS.  LOCAL */
/*   EXTRAPOLATION IS USED TO IMPROVE ABSOLUTE STABILITY AND ACCURACY. */
/*   THE CODE ADJUSTS ITS ORDER AND STEP SIZE TO CONTROL THE LOCAL ERROR */
/*   PER UNIT STEP IN A GENERALIZED SENSE.  SPECIAL DEVICES ARE INCLUDED */
/*   TO CONTROL ROUNDOFF ERROR AND TO DETECT WHEN THE USER IS REQUESTING */
/*   TOO MUCH ACCURACY. */

/*   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, */
/*   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL */
/*   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON. */
/*   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING */
/*   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, */
/*   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060. */


/*   THE PARAMETERS REPRESENT -- */
/*      F -- SUBROUTINE TO EVALUATE DERIVATIVES */
/*      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED */
/*      Y(*) -- SOLUTION VECTOR AT X */
/*      X -- INDEPENDENT VARIABLE */
/*      H -- APPROPRIATE STEP SIZE FOR NEXT STEP.  NORMALLY DETERMINED BY */
/*           CODE */
/*      EPS -- LOCAL ERROR TOLERANCE */
/*      WT(*) -- VECTOR OF WEIGHTS FOR ERROR CRITERION */
/*      START -- LOGICAL VARIABLE SET .TRUE. FOR FIRST STEP,  .FALSE. */
/*           OTHERWISE */
/*      HOLD -- STEP SIZE USED FOR LAST SUCCESSFUL STEP */
/*      K -- APPROPRIATE ORDER FOR NEXT STEP (DETERMINED BY CODE) */
/*      KOLD -- ORDER USED FOR LAST SUCCESSFUL STEP */
/*      CRASH -- LOGICAL VARIABLE SET .TRUE. WHEN NO STEP CAN BE TAKEN, */
/*           .FALSE. OTHERWISE. */
/*      YP(*) -- DERIVATIVE OF SOLUTION VECTOR AT  X  AFTER SUCCESSFUL */
/*           STEP */
/*   THE ARRAYS  PHI, PSI  ARE REQUIRED FOR THE INTERPOLATION SUBROUTINE */
/*   INTRP .  THE ARRAY  P  IS INTERNAL TO THE CODE.  THE REMAINING NINE */
/*   VARIABLES AND ARRAYS ARE INCLUDED IN THE CALL LIST ONLY TO ELIMINATE */
/*   LOCAL RETENTION OF VARIABLES BETWEEN CALLS. */

/*   INPUT TO STEP1 */

/*      FIRST CALL -- */

/*   THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR ALL ARRAYS */
/*   IN THE CALL LIST, NAMELY */

/*     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12), */
/*    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13) */
/*                              --                --    **NOTE** */

/*   THE USER MUST ALSO DECLARE  START ,  CRASH ,  PHASE1  AND  NORND */
/*   LOGICAL VARIABLES AND  F  AN EXTERNAL SUBROUTINE, SUPPLY THE */
/*   SUBROUTINE  F(X,Y,YP)  TO EVALUATE */
/*      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN)) */
/*   AND INITIALIZE ONLY THE FOLLOWING PARAMETERS. */
/*      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED */
/*      Y(*) -- VECTOR OF INITIAL VALUES OF DEPENDENT VARIABLES */
/*      X -- INITIAL VALUE OF THE INDEPENDENT VARIABLE */
/*      H -- NOMINAL STEP SIZE INDICATING DIRECTION OF INTEGRATION */
/*           AND MAXIMUM SIZE OF STEP.  MUST BE VARIABLE */
/*      EPS -- LOCAL ERROR TOLERANCE PER STEP.  MUST BE VARIABLE */
/*      WT(*) -- VECTOR OF NON-ZERO WEIGHTS FOR ERROR CRITERION */
/*      START -- .TRUE. */

/*   STEP1  REQUIRES THAT THE L2 NORM OF THE VECTOR WITH COMPONENTS */
/*   LOCAL ERROR(L)/WT(L)  BE LESS THAN  EPS  FOR A SUCCESSFUL STEP.  THE */
/*   ARRAY  WT  ALLOWS THE USER TO SPECIFY AN ERROR TEST APPROPRIATE */
/*   FOR HIS PROBLEM.  FOR EXAMPLE, */
/*      WT(L) = 1.0  SPECIFIES ABSOLUTE ERROR, */
/*            = ABS(Y(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF THE */
/*                 L-TH COMPONENT OF THE SOLUTION, */
/*            = ABS(YP(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF */
/*                 THE L-TH COMPONENT OF THE DERIVATIVE, */
/*            = AMAX1(WT(L),ABS(Y(L)))  ERROR RELATIVE TO THE LARGEST */
/*                 MAGNITUDE OF L-TH COMPONENT OBTAINED SO FAR, */
/*            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  SPECIFIES A MIXED */
/*                 RELATIVE-ABSOLUTE TEST WHERE  RELERR  IS RELATIVE */
/*                 ERROR,  ABSERR  IS ABSOLUTE ERROR AND  EPS = */
/*                 AMAX1(RELERR,ABSERR) . */

/*      SUBSEQUENT CALLS -- */

/*   SUBROUTINE  STEP1  IS DESIGNED SO THAT ALL INFORMATION NEEDED TO */
/*   CONTINUE THE INTEGRATION, INCLUDING THE STEP SIZE  H  AND THE ORDER */
/*   K , IS RETURNED WITH EACH STEP.  WITH THE EXCEPTION OF THE STEP */
/*   SIZE, THE ERROR TOLERANCE, AND THE WEIGHTS, NONE OF THE PARAMETERS */
/*   SHOULD BE ALTERED.  THE ARRAY  WT  MUST BE UPDATED AFTER EACH STEP */
/*   TO MAINTAIN RELATIVE ERROR TESTS LIKE THOSE ABOVE.  NORMALLY THE */
/*   INTEGRATION IS CONTINUED JUST BEYOND THE DESIRED ENDPOINT AND THE */
/*   SOLUTION INTERPOLATED THERE WITH SUBROUTINE  INTRP .  IF IT IS */
/*   IMPOSSIBLE TO INTEGRATE BEYOND THE ENDPOINT, THE STEP SIZE MAY BE */
/*   REDUCED TO HIT THE ENDPOINT SINCE THE CODE WILL NOT TAKE A STEP */
/*   LARGER THAN THE  H  INPUT.  CHANGING THE DIRECTION OF INTEGRATION, */
/*   I.E., THE SIGN OF  H , REQUIRES THE USER SET  START = .TRUE. BEFORE */
/*   CALLING  STEP1  AGAIN.  THIS IS THE ONLY SITUATION IN WHICH  START */
/*   SHOULD BE ALTERED. */

/*   OUTPUT FROM STEP1 */

/*      SUCCESSFUL STEP -- */

/*   THE SUBROUTINE RETURNS AFTER EACH SUCCESSFUL STEP WITH  START  AND */
/*   CRASH  SET .FALSE. .  X  REPRESENTS THE INDEPENDENT VARIABLE */
/*   ADVANCED ONE STEP OF LENGTH  HOLD  FROM ITS VALUE ON INPUT AND  Y */
/*   THE SOLUTION VECTOR AT THE NEW VALUE OF  X .  ALL OTHER PARAMETERS */
/*   REPRESENT INFORMATION CORRESPONDING TO THE NEW  X  NEEDED TO */
/*   CONTINUE THE INTEGRATION. */

/*      UNSUCCESSFUL STEP -- */

/*   WHEN THE ERROR TOLERANCE IS TOO SMALL FOR THE MACHINE PRECISION, */
/*   THE SUBROUTINE RETURNS WITHOUT TAKING A STEP AND  CRASH = .TRUE. . */
/*   AN APPROPRIATE STEP SIZE AND ERROR TOLERANCE FOR CONTINUING ARE */
/*   ESTIMATED AND ALL OTHER INFORMATION IS RESTORED AS UPON INPUT */
/*   BEFORE RETURNING.  TO CONTINUE WITH THE LARGER TOLERANCE, THE USER */
/*   JUST CALLS THE CODE AGAIN.  A RESTART IS NEITHER REQUIRED NOR */
/*   DESIRABLE. */

/* *********************************************************************** */
/* *  THE ONLY MACHINE DEPENDENT CONSTANTS ARE BASED ON THE MACHINE UNIT * */
/* *  ROUNDOFF ERROR  U  WHICH IS THE SMALLEST POSITIVE NUMBER SUCH THAT * */
/* *  1.0+U .GT. 1.0  .  THE USER MUST CALCULATE  U  AND INSERT          * */
/* *  TWOU=2.0*U  AND  FOURU=4.0*U  IN THE DATA STATEMENT BEFORE CALLING * */
/* *  THE CODE.  THE ROUTINE  MACHIN  CALCULATES  U .                    * */


int step1(S_fp f, int *neqn, double *y, double *x, 
	double *h__, double *eps, double *wt, logical *start, 
	double *hold, int *k, int *kold, logical *crash, 
	double *phi, double *p, double *yp, double *psi, 
	double *alpha, double *beta, double *sig, double *v, 
	double *w, double *g, logical *phase1, int *ns, logical *nornd)
{
    /* Initialized data */
    static double twou = 1.54e-8;
    static double fouru = 3.08e-8;
    static double two[13] = { 2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,2048.,4096.,8192. };
    static double gstr[13] = { .5,.0833,.0417,.0264,.0188,.0143,.0114,.00936,.00789,.00679,.00592,.00524,.00468 };

    /* System generated locals */
    int phi_dim1, phi_offset, i1, i2;
    double d1, d2, d3;

    /* Builtin functions */
    double d_sign(double *, double *), sqrt(double), pow_dd(double *, double *);

    /* Local variables */
    static int i__, j, l;
    static double r__;
    static int iq, im1, km1, km2, ip1, kp1, kp2;
    static double erk, err, tau, rho, sum;
    static int nsm2, nsp1, nsp2;
    static double absh, hnew;
    static int knew;
    static double xold, erkm1, erkm2, erkp1, temp1, temp2, temp3, temp4, 
	    temp5, temp6, p5eps;
    static int ifail;
    static double reali, round;
    static int limit1, limit2;
    static double realns;


    /* Parameter adjustments */
    --yp;
    --p;
    phi_dim1 = *neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --y;
    --psi;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;

    /* Function Body */
/* *********************************************************************** */


/*       ***     BEGIN BLOCK 0     *** */
/*   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE */
/*   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A */
/*   STARTING STEP SIZE. */
/*                   *** */

/*   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE */

    *crash = TRUE_;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L5;
    }
    d1 = fouru * abs(*x);
    *h__ = d_sign(&d1, h__);
    return 0;
L5:
    p5eps = *eps * .5;

/*   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE */

    round = 0.;
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
/* L10: */
/* Computing 2nd power */
	d1 = y[l] / wt[l];
	round += d1 * d1;
    }
    round = twou * sqrt(round);
    if (p5eps >= round) {
	goto L15;
    }
    *eps = round * 2. * (fouru + 1.);
    return 0;
L15:
    *crash = FALSE_;
    g[1] = 1.;
    g[2] = .5;
    sig[1] = 1.;
    if (! (*start)) {
	goto L99;
    }

/*   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP */

    (*f)(x, &y[1], &yp[1]);
    sum = 0.;
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	phi[l + phi_dim1] = yp[l];
	phi[l + (phi_dim1 << 1)] = 0.;
/* L20: */
/* Computing 2nd power */
	d1 = yp[l] / wt[l];
	sum += d1 * d1;
    }
    sum = sqrt(sum);
    absh = abs(*h__);
    if (*eps < sum * 16. * *h__ * *h__) {
	absh = sqrt(*eps / sum) * .25;
    }
/* Computing MAX */
    d2 = absh, d3 = fouru * abs(*x);
    d1 = max(d2,d3);
    *h__ = d_sign(&d1, h__);
    *hold = 0.;
    *k = 1;
    *kold = 0;
    *start = FALSE_;
    *phase1 = TRUE_;
    *nornd = TRUE_;
    if (p5eps > round * 100.) {
	goto L99;
    }
    *nornd = FALSE_;
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
/* L25: */
	phi[l + phi_dim1 * 15] = 0.;
    }
L99:
    ifail = 0;
/*       ***     END BLOCK 0     *** */

/*       ***     BEGIN BLOCK 1     *** */
/*   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING */
/*   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED. */
/*                   *** */

L100:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    km2 = *k - 2;

/*   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT */
/*   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE */

    if (*h__ != *hold) {
	*ns = 0;
    }
    if (*ns <= *kold) {
	++(*ns);
    }
    nsp1 = *ns + 1;
    if (*k < *ns) {
	goto L199;
    }

/*   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH */
/*   ARE CHANGED */

    beta[*ns] = 1.;
    realns = (double) (*ns);
    alpha[*ns] = 1. / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.;
    if (*k < nsp1) {
	goto L110;
    }
    i1 = *k;
    for (i__ = nsp1; i__ <= i1; ++i__) {
	im1 = i__ - 1;
	temp2 = psi[im1];
	psi[im1] = temp1;
	beta[i__] = beta[im1] * psi[im1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	reali = (double) i__;
/* L105: */
	sig[i__ + 1] = reali * alpha[i__] * sig[i__];
    }
L110:
    psi[*k] = temp1;

/*   COMPUTE COEFFICIENTS G(*) */

/*   INITIALIZE V(*) AND SET W(*). */

    if (*ns > 1) {
	goto L120;
    }
    i1 = *k;
    for (iq = 1; iq <= i1; ++iq) {
	temp3 = (double) (iq * (iq + 1));
	v[iq] = 1. / temp3;
/* L115: */
	w[iq] = v[iq];
    }
    goto L140;

/*   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*) */

L120:
    if (*k <= *kold) {
	goto L130;
    }
    temp4 = (double) (*k * kp1);
    v[*k] = 1. / temp4;
    nsm2 = *ns - 2;
    if (nsm2 < 1) {
	goto L130;
    }
    i1 = nsm2;
    for (j = 1; j <= i1; ++j) {
	i__ = *k - j;
/* L125: */
	v[i__] -= alpha[j + 1] * v[i__ + 1];
    }

/*   UPDATE V(*) AND SET W(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i1 = limit1;
    for (iq = 1; iq <= i1; ++iq) {
	v[iq] -= temp5 * v[iq + 1];
/* L135: */
	w[iq] = v[iq];
    }
    g[nsp1] = w[1];

/*   COMPUTE THE G(*) IN THE WORK VECTOR W(*) */

L140:
    nsp2 = *ns + 2;
    if (kp1 < nsp2) {
	goto L199;
    }
    i1 = kp1;
    for (i__ = nsp2; i__ <= i1; ++i__) {
	limit2 = kp2 - i__;
	temp6 = alpha[i__ - 1];
	i2 = limit2;
	for (iq = 1; iq <= i2; ++iq) {
/* L145: */
	    w[iq] -= temp6 * w[iq + 1];
	}
/* L150: */
	g[i__] = w[1];
    }
L199:
/*       ***     END BLOCK 1     *** */

/*       ***     BEGIN BLOCK 2     *** */
/*   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED */
/*   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K, */
/*   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED. */
/*                   *** */

/*   CHANGE PHI TO PHI STAR */

    if (*k < nsp1) {
	goto L215;
    }
    i1 = *k;
    for (i__ = nsp1; i__ <= i1; ++i__) {
	temp1 = beta[i__];
	i2 = *neqn;
	for (l = 1; l <= i2; ++l) {
/* L205: */
	    phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
	}
/* L210: */
    }

/*   PREDICT SOLUTION AND DIFFERENCES */

L215:
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
	phi[l + kp1 * phi_dim1] = 0.;
/* L220: */
	p[l] = 0.;
    }
    i1 = *k;
    for (j = 1; j <= i1; ++j) {
	i__ = kp1 - j;
	ip1 = i__ + 1;
	temp2 = g[i__];
	i2 = *neqn;
	for (l = 1; l <= i2; ++l) {
	    p[l] += temp2 * phi[l + i__ * phi_dim1];
/* L225: */
	    phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
	}
/* L230: */
    }
    if (*nornd) {
	goto L240;
    }
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
	p[l] = y[l] + tau;
/* L235: */
	phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
/* L245: */
	p[l] = y[l] + *h__ * p[l];
    }
L250:
    xold = *x;
    *x += *h__;
    absh = abs(*h__);
    (*f)(x, &p[1], &yp[1]);

/*   ESTIMATE ERRORS AT ORDERS K,K-1,K-2 */

    erkm2 = 0.;
    erkm1 = 0.;
    erk = 0.;
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	temp3 = 1. / wt[l];
	temp4 = yp[l] - phi[l + phi_dim1];
	if (km2 < 0) {
	    goto L265;
	} else if (km2 == 0) {
	    goto L260;
	} else {
	    goto L255;
	}
L255:
/* Computing 2nd power */
	d1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
	erkm2 += d1 * d1;
L260:
/* Computing 2nd power */
	d1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
	erkm1 += d1 * d1;
L265:
/* Computing 2nd power */
	d1 = temp4 * temp3;
	erk += d1 * d1;
    }
    if (km2 < 0) {
	goto L280;
    } else if (km2 == 0) {
	goto L275;
    } else {
	goto L270;
    }
L270:
    erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
L275:
    erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
L280:
    temp5 = absh * sqrt(erk);
    err = temp5 * (g[*k] - g[kp1]);
    erk = temp5 * sig[kp1] * gstr[*k - 1];
    knew = *k;

/*   TEST IF ORDER SHOULD BE LOWERED */

    if (km2 < 0) {
	goto L299;
    } else if (km2 == 0) {
	goto L290;
    } else {
	goto L285;
    }
L285:
    if (max(erkm1,erkm2) <= erk) {
	knew = km1;
    }
    goto L299;
L290:
    if (erkm1 <= erk * .5) {
	knew = km1;
    }

/*   TEST IF STEP SUCCESSFUL */

L299:
    if (err <= *eps) {
	goto L400;
    }
/*       ***     END BLOCK 2     *** */

/*       ***     BEGIN BLOCK 3     *** */
/*   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) . */
/*   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE */
/*   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR */
/*   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE */
/*   PRECISION. */
/*                   *** */

/*   RESTORE X, PHI(*,*) AND PSI(*) */

    *phase1 = FALSE_;
    *x = xold;
    i1 = *k;
    for (i__ = 1; i__ <= i1; ++i__) {
	temp1 = 1. / beta[i__];
	ip1 = i__ + 1;
	i2 = *neqn;
	for (l = 1; l <= i2; ++l) {
/* L305: */
	    phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[
		    l + ip1 * phi_dim1]);
	}
/* L310: */
    }
    if (*k < 2) {
	goto L320;
    }
    i1 = *k;
    for (i__ = 2; i__ <= i1; ++i__) {
/* L315: */
	psi[i__ - 1] = psi[i__] - *h__;
    }

/*   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP */
/*   SIZE */

L320:
    ++ifail;
    temp2 = .5;
    if ((i1 = ifail - 3) < 0) {
	goto L335;
    } else if (i1 == 0) {
	goto L330;
    } else {
	goto L325;
    }
L325:
    if (p5eps < erk * .25) {
	temp2 = sqrt(p5eps / erk);
    }
L330:
    knew = 1;
L335:
    *h__ = temp2 * *h__;
    *k = knew;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L340;
    }
    *crash = TRUE_;
    d1 = fouru * abs(*x);
    *h__ = d_sign(&d1, h__);
    *eps += *eps;
    return 0;
L340:
    goto L100;
/*       ***     END BLOCK 3     *** */

/*       ***     BEGIN BLOCK 4     *** */
/*   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE */
/*   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE */
/*   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP. */
/*                   *** */
L400:
    *kold = *k;
    *hold = *h__;

/*   CORRECT AND EVALUATE */

    temp1 = *h__ * g[kp1];
    if (*nornd) {
	goto L410;
    }
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
	y[l] = p[l] + rho;
/* L405: */
	phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
    }
    goto L420;
L410:
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
/* L415: */
	y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
    }
L420:
    (*f)(x, &y[1], &yp[1]);

/*   UPDATE DIFFERENCES FOR NEXT STEP */

    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
	phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
/* L425: */
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * 
		phi_dim1];
    }
    i1 = *k;
    for (i__ = 1; i__ <= i1; ++i__) {
	i2 = *neqn;
	for (l = 1; l <= i2; ++l) {
/* L430: */
	    phi[l + i__ * phi_dim1] += phi[l + kp1 * phi_dim1];
	}
/* L435: */
    }

/*   ESTIMATE ERROR AT ORDER K+1 UNLESS: */
/*     IN FIRST PHASE WHEN ALWAYS RAISE ORDER, */
/*     ALREADY DECIDED TO LOWER ORDER, */
/*     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE */

    erkp1 = 0.;
    if (knew == km1 || *k == 12) {
	*phase1 = FALSE_;
    }
    if (*phase1) {
	goto L450;
    }
    if (knew == km1) {
	goto L455;
    }
    if (kp1 > *ns) {
	goto L460;
    }
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) {
/* L440: */
/* Computing 2nd power */
	d1 = phi[l + kp2 * phi_dim1] / wt[l];
	erkp1 += d1 * d1;
    }
    erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

/*   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER */
/*   FOR NEXT STEP */

    if (*k > 1) {
	goto L445;
    }
    if (erkp1 >= erk * .5) {
	goto L460;
    }
    goto L450;
L445:
    if (erkm1 <= min(erk,erkp1)) {
	goto L455;
    }
    if (erkp1 >= erk || *k == 12) {
	goto L460;
    }

/*   HERE ERKP1 .LT. ERK .LT. AMAX1(ERKM1,ERKM2) ELSE ORDER WOULD HAVE */
/*   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED */

/*   RAISE ORDER */

L450:
    *k = kp1;
    erk = erkp1;
    goto L460;

/*   LOWER ORDER */

L455:
    *k = km1;
    erk = erkm1;

/*   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP */

L460:
    hnew = *h__ + *h__;
    if (*phase1) {
	goto L465;
    }
    if (p5eps >= erk * two[*k]) {
	goto L465;
    }
    hnew = *h__;
    if (p5eps >= erk) {
	goto L465;
    }
    temp2 = (double) (*k + 1);
    d1 = p5eps / erk;
    d2 = 1. / temp2;
    r__ = pow_dd(&d1, &d2);
/* Computing MAX */
    d1 = .5, d2 = min(.9,r__);
    hnew = absh * max(d1,d2);
/* Computing MAX */
    d2 = hnew, d3 = fouru * abs(*x);
    d1 = max(d2,d3);
    hnew = d_sign(&d1, h__);
L465:
    *h__ = hnew;
    return 0;
/*       ***     END BLOCK 4     *** */
} /* step1_ */

