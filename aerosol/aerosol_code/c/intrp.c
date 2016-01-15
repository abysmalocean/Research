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


/*   THE METHODS IN SUBROUTINE  STEP1  APPROXIMATE THE SOLUTION NEAR  X */
/*   BY A POLYNOMIAL.  SUBROUTINE  INTRP  APPROXIMATES THE SOLUTION AT */
/*   XOUT  BY EVALUATING THE POLYNOMIAL THERE.  INFORMATION DEFINING THIS */
/*   POLYNOMIAL IS PASSED FROM  STEP1  SO  INTRP  CANNOT BE USED ALONE. */

/*   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, */
/*   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL */
/*   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON. */
/*   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING */
/*   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, */
/*   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060. */

/*   INPUT TO INTRP -- */

/*   THE USER PROVIDES STORAGE IN THE CALLING PROGRAM FOR THE ARRAYS IN */
/*   THE CALL LIST */
/*      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),PSI(12) */
/*   AND DEFINES */
/*      XOUT -- POINT AT WHICH SOLUTION IS DESIRED. */
/*   THE REMAINING PARAMETERS ARE DEFINED IN  STEP1  AND PASSED TO */
/*   INTRP  FROM THAT SUBROUTINE */

/*   OUTPUT FROM  INTRP -- */

/*      YOUT(*) -- SOLUTION AT  XOUT */
/*      YPOUT(*) -- DERIVATIVE OF SOLUTION AT  XOUT */
/*   THE REMAINING PARAMETERS ARE RETURNED UNALTERED FROM THEIR INPUT */
/*   VALUES.  INTEGRATION WITH  STEP1  MAY BE CONTINUED. */


int intrp(double *x, double *y, double *xout, double *yout, double *ypout, int *neqn, int *kold,double *phi, double *psi)
{
    /* Initialized data */
    static double g[13] = { 1.0 };
    static double rho[13] = { 1.0 };
    /* System generated locals */
    int phi_dim1, phi_offset, i1, i2;
    /* Local variables */
    static int i, j, l;
    static double w[13], hi;
    static int ki, jm1;
    static double eta;
    static int kip1;
    static double term, temp1, temp2, temp3, gamma;
    static int limit1;
    static double psijm1;

    /* Parameter adjustments */
    phi_dim1 = *neqn;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --ypout;
    --yout;
    --y;
    --psi;

    /* Function Body */
    hi = *xout - *x;
    ki = *kold + 1;
    kip1 = ki + 1;

    /*INITIALIZE W(*) FOR COMPUTING G(*)*/
    i1 = ki;
    for (i = 1; i <= i1; ++i) 
    {
	temp1 = (double) i;
	w[i - 1] = 1.0 / temp1;
    }
    term = 0.0;

    /*COMPUTE G(*)*/
    i1 = ki;
    for (j = 2; j <= i1; ++j) 
    {
	jm1 = j - 1;
	psijm1 = psi[jm1];
	gamma = (hi + term) / psijm1;
	eta = hi / psijm1;
	limit1 = kip1 - j;
	i2 = limit1;
	for (i = 1; i <= i2; ++i) 
	{
	    w[i-1] = gamma*w[i-1]-eta*w[i];
	}
	g[j-1] = w[0];
	rho[j-1] = gamma*rho[jm1-1];
	term = psijm1;
    }

    /*INTERPOLATE*/
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) 
    {
	ypout[l] = 0.0;
	yout[l] = 0.0;
    }
    i1 = ki;
    for (j = 1; j <= i1; ++j) 
    {
	i = kip1 - j;
	temp2 = g[i - 1];
	temp3 = rho[i - 1];
	i2 = *neqn;
	for (l = 1; l <= i2; ++l) 
	{
	    yout[l] += temp2 * phi[l + i * phi_dim1];
	    ypout[l] += temp3 * phi[l + i * phi_dim1];
	}
    }
    i1 = *neqn;
    for (l = 1; l <= i1; ++l) 
    {
	yout[l] = y[l] + hi * yout[l];
    }
    return 0;
} 
