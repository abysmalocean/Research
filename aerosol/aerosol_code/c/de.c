 /*   ODE  MERELY ALLOCATES STORAGE FOR  DE  TO RELIEVE THE USER OF THE */
/*   INCONVENIENCE OF A LONG CALL LIST.  CONSEQUENTLY  DE  IS USED AS */
/*   DESCRIBED IN THE COMMENTS FOR  ODE . */

/*   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, */
/*   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL */
/*   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON. */


/* *********************************************************************** */
/* *  THE ONLY MACHINE DEPENDENT CONSTANT IS BASED ON THE MACHINE UNIT   * */
/* *  ROUNDOFF ERROR  U  WHICH IS THE SMALLEST POSITIVE NUMBER SUCH THAT * */
/* *  1.0+U .GT. 1.0 .  U  MUST BE CALCULATED AND  FOURU=4.0*U  INSERTED * */
/* *  IN THE FOLLOWING DATA STATEMENT BEFORE USING  DE .  THE ROUTINE    * */
/* *  MACHIN  CALCULATES  U .  FOURU  AND  TWOU=2.0*U  MUST ALSO BE      * */
/* *  INSERTED IN SUBROUTINE  STEP1 BEFORE CALLING  DE .                 * */
      

int de(S_fp f, int *neqn, double *y, double *t, 
	double *tout, double *relerr, double *abserr, int *
	iflag, double *yy, double *wt, double *p, double *yp, 
	double *ypout, double *phi, double *alpha, double *
	beta, double *sig, double *v, double *w, double *g, 
	logical *phase1, double *psi, double *x, double *h__, 
	double *hold, logical *start, double *told, double *
	delsgn, int *ns, logical *nornd, int *k, int *kold, 
	int *isnold)
{
    /* Initialized data */

    static double fouru = 3.08e-8;
    static int maxnum = 500;

    /* System generated locals */
    int phi_dim1, phi_offset, i__1;
    double d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    int i_sign(int *, int *);
    double d_sign(double *, double *);

    /* Local variables */
    static int l;
    static double del, eps;
    static int isn, kle4;
    static double tend;
    extern /* Subroutine */ int step1_(S_fp, int *, double *, 
	    double *, double *, double *, double *, logical *,
	     double *, int *, int *, logical *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, logical *, int *, logical *);
    static logical crash, stiff;
    extern int intrp(double *, double *, double *, double *, double *, int *, int *, double *, double *);
    static double absdel, abseps;
    extern int errchk(int *, int *);
    static double releps;
    static int nostep;

   /* Parameter adjustments */
    phi_dim1 = *neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --ypout;
    --yp;
    --p;
    --wt;
    --yy;
    --y;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;
    --psi;

    /* Function Body */
/* *********************************************************************** */

/*   THE CONSTANT  MAXNUM  IS THE MAXIMUM NUMBER OF STEPS ALLOWED IN ONE */
/*   CALL TO  DE .  THE USER MAY CHANGE THIS LIMIT BY ALTERING THE */
/*   FOLLOWING STATEMENT */

/*            ***            ***            *** */
/*   TEST FOR IMPROPER PARAMETERS */

    if (*neqn < 1) 
    {
	goto L10;
    }
    if (*t == *tout) 
    {
	goto L10;
    }
    if (*relerr < (float)0. || *abserr < 0.) 
    {
	goto L10;
    }
    eps = max(*relerr,*abserr);
    if (eps <= 0.) 
    {
	goto L10;
    }
    if (*iflag == 0) 
    {
	goto L10;
    }
    isn = i_sign(&c__1, iflag);
    *iflag = abs(*iflag);
    if (*iflag == 1) 
    {
	goto L20;
    }
    if (*t != *told) 
    {
	goto L10;
    }
    if (*iflag >= 2 && *iflag <= 5) 
    {
	goto L20;
    }
    if (*iflag == 6 && *abserr > 0.) 
    {
	goto L20;
    }
    if (*iflag != 6 && *iflag != 7) 
    {
	goto L10;
    }
    if (*iflag == 6) 
    {
	errchk(&c_n78, (int*)"IN ODE, IFLAG = 6 AND ABSERR NON-POSITIVE ON INPUT -- EXECUTION TERMINATED.");
    }
    if (*iflag == 7) 
    {
	errchk(&c_n54, (int*)"IN ODE, IFLAG = 7 ON INPUT -- EXECUTION TERMINATED.", (ftnlen)54);
    }
    return 0;
L10:
    *iflag = 7;
    errchk(&c__36, (int*)"IN ODE, INVALID INPUT PARAMETERS.", (ftnlen)36);
    return 0;

/*   ON EACH CALL SET INTERVAL OF INTEGRATION AND COUNTER FOR NUMBER OF */
/*   STEPS.  ADJUST INPUT ERROR TOLERANCES TO DEFINE WEIGHT VECTOR FOR */
/*   SUBROUTINE  STEP1 */

L20:
    del = *tout - *t;
    absdel = abs(del);
    tend = *t + del * 10.;
    if (isn < 0) 
    {
	tend = *tout;
    }
    nostep = 0;
    kle4 = 0;
    stiff = FALSE_;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1) 
    {
	goto L30;
    }
    if (*isnold < 0) 
    {
	goto L30;
    }
    if (*delsgn * del > 0.) 
    {
	goto L50;
    }

/*   ON START AND RESTART ALSO SET WORK VARIABLES X AND YY(*), STORE THE */
/*   DIRECTION OF INTEGRATION AND INITIALIZE THE STEP SIZE */

L30:
    *start = TRUE_;
    *x = *t;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
/* L40: */
	yy[l] = y[l];
    }
    *delsgn = d_sign(&c_b174, &del);
/* Computing MAX */
    d__3 = fouru * abs(*x), d__4 = (d__1 = *tout - *x, abs(d__1));
    d__2 = max(d__3,d__4);
    d__5 = *tout - *x;
    *h__ = d_sign(&d__2, &d__5);

/*   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN */

L50:
    if ((d__1 = *x - *t, abs(d__1)) < absdel) 
    {
	goto L60;
    }
    intrp_(x, &yy[1], tout, &y[1], &ypout[1], neqn, kold, &phi[phi_offset], &psi[1]);
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   IF CANNOT GO PAST OUTPUT POINT AND SUFFICIENTLY CLOSE, */
/*   EXTRAPOLATE AND RETURN */

L60:
    if (isn > 0 || (d__1 = *tout - *x, abs(d__1)) >= fouru * abs(*x)) 
    {
	goto L80;
    }
    *h__ = *tout - *x;
    (*f)(x, &yy[1], &yp[1]);
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
/* L70: */
	y[l] = yy[l] + *h__ * yp[l];
    }
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   TEST FOR TOO MANY STEPS */

L80:
    if (nostep < maxnum) 
    {
	goto L100;
    }
    *iflag = isn << 2;
    if (stiff) 
    {
	*iflag = isn * 5;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
/* L90: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP */

L100:
/* Computing MIN */
    d__3 = abs(*h__), d__4 = (d__1 = tend - *x, abs(d__1));
    d__2 = min(d__3,d__4);
    *h__ = d_sign(&d__2, h__);
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
	wt[l] = releps * (d__1 = yy[l], abs(d__1)) + abseps;
	if (wt[l] <= 0.) 
        {
	    goto L140;
	}
/* L110: */
    }
    step1_((S_fp)f, neqn, &yy[1], x, h__, &eps, &wt[1], start, hold, k, kold,
	 &crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], 
	&beta[1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd);

/*   TEST FOR TOLERANCES TOO SMALL */

    if (! crash) 
    {
	goto L130;
    }
    *iflag = isn * 3;
    *relerr = eps * releps;
    *abserr = eps * abseps;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
/* L120: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   AUGMENT COUNTER ON NUMBER OF STEPS AND TEST FOR STIFFNESS */

L130:
    ++nostep;
    ++kle4;
    if (*kold > 4) 
    {
	kle4 = 0;
    }
    if (kle4 >= 50) 
    {
	stiff = TRUE_;
    }
    goto L50;

/*   RELATIVE ERROR CRITERION INAPPROPRIATE */

L140:
    *iflag = isn * 6;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) 
    {
/* L150: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    errchk(&c__97, (int*)"IN ODE, PURE RELATIVE ERROR IMPOSSIBLE WHEN SOLUTION VANISHES. USE A NON-ZERO VALUE OF ABSERR", (ftnlen)97);
    return 0;
}

