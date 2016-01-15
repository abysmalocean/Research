/* ----------------------------------------------------------------------- 		 */
/* *    THIS FUNCTION CALCULATES THE DERIVATIVE OF THE 'M' DISTRIBUTION  		 */
/* *    FUNCTION W.R.T. W, FOR THE NVAR GRID POINTS                      		 */
/* *      REF. SEE REFERENCES TO ROUTINE DISTW                           		 */
/* *                                                                     		 */
/* *              VARIABLES                                              		 */
/* *    A(I)='M' DISTRIBUTION AT W(I)                                    		 */
/* *    DERIV=OUTPUT ARRAY OF DERIVATIVES AT GRID POINTS                 		 */
/* *    NVAR=NUMBER OF GRID POINTS  (INPUT)                              		 */
/* *    W(I)=LOGARITHMICALLY TRANSFORMED PARTICLE SIZE AT THE I-TH GRID POINT 		 */
/* *    X(I)=PARTICLE SIZE AT THE I-TH GRID POINT                        		 */
/*      P array contains the second derivatives at the grid points and passed to distw by the common block aersl7 */     

int derivt(int *isplin, int *nvar, double *deriv)
{
    /* System generated locals */
    int index1, i, nvarm1, nvarm2;
    /* Parameter adjustments */
    --deriv;
    /* Function Body */
    nvarm1 = *nvar - 1;
    nvarm2 = *nvar - 2;
    deriv[1]=(aersl2_1.a[1]-aersl2_1.a[0])/(aersl2_1.w[c__1]-aersl2_1.w[c__1-1])-(aersl2_1.w[c__1]-aersl2_1.w[c__1-1])*aersl7_1.p[0]/6.0;
    index1 = nvarm2;
    for (i = 2; i <= nvarm2; ++i) 
    {
	deriv[i]=(aersl2_1.a[i]-aersl2_1.a[i-1])/(aersl2_1.w[i]-aersl2_1.w[i-1])-(aersl2_1.w[i]-aersl2_1.w[i-1])*(aersl7_1.p[i-2]*2.0+aersl7_1.p[i-1])/6.0;
    }
    deriv[nvarm1]=(aersl2_1.a[*nvar-1]-aersl2_1.a[nvarm1-1])/(aersl2_1.w[nvarm1]-aersl2_1.w[nvarm1-1])-(aersl2_1.w[nvarm1]
	     - aersl2_1.w[nvarm1 - 1]) * aersl7_1.p[nvarm2 - 1]/3.0;
    deriv[*nvar] = deriv[nvarm1]+(aersl2_1.w[nvarm1]-aersl2_1.w[nvarm1-1])*0.5*aersl7_1.p[nvarm2-1];
    return 0;
}

