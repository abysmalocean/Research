/* ----------------------------------------------------------------------- */
/* *    THIS FUNCTION COMPUTES THE VALUES OF THE M(WX,T) DISTRIBUTION    * */
/* *    BY INTERPOLATING FROM THE VALUES OF A AT THE GRID POINTS.        * */
/* *    NOTE THAT A(I) IS EQUAL TO THE 'M' DISTRIBUTION FUNCTION AT W(I) * */
/* *                                                                     * */
/* *            VARIABLES                                                * */
/* *                                                                     * */
/* *    M(WX,TIME)=N(X,TIME)*DERIVATIVE OF X WITH RESPECT TO WX          * */
/* *    NPTS=NUMBER OF POINTS IN WX ARRAY  (INPUT)                       * */
/* *    NVARM1=NUMBER OF GRID POINTS MINUS ONE  (INPUT)                  * */
/* *    WX=LN(X/XA)/LN(XB/XA)   (INPUT ARRAY)                            * */
/* *    X(I)=PARTICLE SIZE AT THE I-TH GRID POINT                        * */
/* *    XA=SMALLEST PARTICLE SIZE IN CONTINUOUS REGIME                   * */
/* *    XB=LARGEST PARTICLE SIZE IN CONTINUOUS REGIME                    * */
/* *    ZDIST(I)=ARRAY OF M(WX,TIME)   (OUTPUT)                          * */
/* *                                                                     * */
/* *    INTERPOLATION METHOD                                             * */
/* *       ISPLIN=0  CUBIC SPLINES                                       * */
/* *                  REF. GERALD,C.F.,'APPLIED NUMERICAL ANALYSIS',     * */
/* *                       PP.474-488, ADDISON-WESLEY 1978               * */
/* *                       ALSO CARNAHAN,B., LUTHER,H.A., AND WILKES,J.O.* */
/* *                       'APPLIED NUMERICAL METHODS',P.63,WILEY 1969   * */
/* *              1  LINEAR SPLINES                                      * */
/* *              2  LOGARITHMIC SPLINES (POWER LAW)                     * */
/* ----------------------------------------------------------------------- */
/*        DETERMINE BOUNDING GRID POINTS OF WX ASSUMING THAT THE POINTS */
/*        WX(I),I=1,NPTS ARE IN ASCENDING ORDER */

/* Subroutine */ 
int distw(int *isplin, int *nvarm1, int *npts, double *wx, double *zdist)
{
    /* System generated locals */
    int index1, index2;
    double d1, d2, d3, d4;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static int i, j, jstart;

    /* Parameter adjustments */
    --zdist;
    --wx;

    /* Function Body */
    jstart = 1;
    index1 = *npts;
    for (i = 1; i <= index1; ++i) 
    {
	index2 = *nvarm1;
	for (j = jstart; j <= index2; ++j) 
	{
	    if (wx[i] <= aersl2_1.w[j]) 
	    {
//		goto L20;
		//L20:
		if (*isplin == 1) 
		{
	    		//goto L1;
			/*LINEAR SPLINES */
			//L1:
			zdist[i] = aersl2_1.a[j]+(wx[i]-aersl2_1.w[j])/(aersl2_1.w[j-1]-aersl2_1.w[j])*(aersl2_1.a[j-1]-aersl2_1.a[j]);
			//goto L3;
			jstart = j;
		}	
		if (*isplin == 2) 
		{
			//goto L2;
			/*LOGARITHMIC SPLINES */
			//L2:
			zdist[i] = log(aersl2_1.a[j]);
			zdist[i] += (wx[i]-aersl2_1.w[j])/(aersl2_1.w[j-1]-aersl2_1.w[j])*(log(aersl2_1.a[j-1])-zdist[i]);
			zdist[i] = exp(zdist[i]);
		}

		/*CUBIC SPLINES */
		if (j == 1) 
		{
	    		//goto L30;
			//L30:
			/* Computing 3rd power */
			d1 = wx[i]-aersl2_1.w[j-1];
			/* Computing 2nd power */
			d2 = aersl2_1.w[j]-aersl2_1.w[j-1];
			zdist[i] = (aersl7_1.p[j-1]*(d1*(d1*d1))+(aersl2_1.a[j]*6.0-d2*d2*aersl7_1.p[j-1])*(wx[i]-aersl2_1.w[j-1])+aersl2_1.a[j-1]*6.0*(aersl2_1.w[j]-wx[i]))/((aersl2_1.w[j]-aersl2_1.w[j-1])*6.0);
			//goto L3;
			jstart = j;
		}
		if (j == *nvarm1) 
		{
	    		//goto L40;
			//L40:
			/* Computing 3rd power */
			d1 = aersl2_1.w[j]-wx[i];
			/*Computing 2nd power */
			d2 = aersl2_1.w[j]-aersl2_1.w[j-1];
			zdist[i] = (aersl7_1.p[j-2]*(d1*(d1*d1))+aersl2_1.a[j]*6.0*(wx[i]-aersl2_1.w[j-1])+(aersl2_1.a[j-1]*6.0-d2*d2*aersl7_1.p[j-2])*(aersl2_1.w[j]-wx[i]))/((aersl2_1.w[j]-aersl2_1.w[j-1])*6.0);
			//goto L3;
			jstart = j;
		}
		/* Computing 3rd power */
		d1 = aersl2_1.w[j]-wx[i];
		/* Computing 3rd power */
		d2 = wx[i]-aersl2_1.w[j-1];
		/* Computing 2nd power */
		d3 = aersl2_1.w[j]-aersl2_1.w[j-1];
		/* Computing 2nd power */
		d4 = aersl2_1.w[j]-aersl2_1.w[j-1];
		zdist[i] = (aersl7_1.p[j-2]*(d1*(d1*d1))+aersl7_1.p[j-1]*(d2*(d2*d2))+(aersl2_1.a[j]*6.0-d3*d3*aersl7_1.p[j-1])*(wx[i]-aersl2_1.w[j-1])+(aersl2_1.a[j-1]*6.0-d4*d4*aersl7_1.p[j-2])*(aersl2_1.w[j]-wx[i]))/((aersl2_1.w[j]-aersl2_1.w[j-1])*6.0);
		//goto L3;
		jstart = j;
	    }
	}
	j = *nvarm1;
	//L3:
	//jstart = j;
    }
    return 0;
} /* distw_ */
