/* ----------------------------------------------------------------------- */
/* *    THIS ROUTINE IS USED FOR DETERMINING THE SECOND DERIVATIVE OF    * */
/* *    THE DISTRIBTUION AT THE GRID POINTS BY SOLVING A SYSTEM OF LINEAR* */
/* *    ALGEBRAIC EQUATIONS WITH A TRI-DIAGONAL COEFFICIENT MATRIX.      * */
/* *    THE ROUTINE IS REQUIRED FOR A CUBIC SPLINE INTERPOLATION.        * */
/* *    REFERENCE; CARNAHAN,B., LUTHER,H.A. AND WILKES,J.O. 'APPLIED     * */
/* *    NUMERICAL METHODS,' JOHN WILEY AND SONS, PAGE 446 (1969)         * */
/* *                                                                     * */
/* *               VARIABLES                                             * */
/* *                                                                     * */
/* *    DIAG(I)=THE I-TH DIAGONAL ELEMENT OF THE MATRIX                  * */
/* *    NEWSET=FLAG TO INDICATE A NEW SET OF GRID POINTS                 * */
/* *           0=SAME SET OF GRID POINTS                                 * */
/* *           1=NEW SET OF GRID POINTS                                  * */
/* *    P(I)=THE SECOND DERIVATIVE AT THE (I+1)-TH GRID POINT            * */
/* *    SUB(I)=THE (I-1)-TH SUB DIAGONAL ELEMENT OF THE MATRIX           * */
/* *    SUPER(I)=THE I-TH SUPER-DIAGONAL ELEMENT OF THE MATRIX           * */
/* ----------------------------------------------------------------------- */

/* Subroutine */ 
int setup(int *isplin, int *newset, int *nvar)
{
    /* System generated locals */
    int i1, i2;
    double d1;

    /* Local variables */
    static double g[38];
    static int i, j;
    static double sub[38], diag[38], super[38];
    static int nvarm2, nvarm3;

    /*COMPUTE TRI-DIAGONAL MATRIX */

    if (*newset == 0) 
    {
	goto L6;
    }
    nvarm2 = *nvar - 2;
    nvarm3 = nvarm2 - 1;
    i1 = nvarm2;
    for (i = 1; i <= i1; ++i) 
    {
    /* Computing 2nd power */
	d1 = aersl2_1.w[i] - aersl2_1.w[i - 1];
	sub[i - 1] = d1 * d1 / 6.0;
	i2 = i + 1;
	super[i - 1] = (aersl2_1.w[i] - aersl2_1.w[i - 1]) * (aersl2_1.w[i2] - aersl2_1.w[i2 - 1]) / 6.0;
    /* L1: */
	diag[i - 1] = (sub[i - 1] + super[i - 1]) * 2.;
    }
    i1 = nvarm2;
    for (i = 2; i <= i1; ++i) 
    {
    /* L2: */
	diag[i - 1] -= sub[i - 1] * super[i - 2] / diag[i - 2];
    }

    /*SOLVE THE LINEAR SYSTEM WITH THE TRI-DIAGONAL MATRIX */

L6:
    i1 = varm2;
    for (i = 1; i <= i1; ++i) 
    {
/*L5:*/
	i2 = i + 1;
	aersl7_1.p[i-1]=(aersl2_1.w[i]-aersl2_1.w[i-1])*(aersl2_1.a[i+1]-aersl2_1.a[i])/(aersl2_1.w[i2]-aersl2_1.w[i2-1])-aersl2_1.a[i]+aersl2_1.a[i-1];
    }
    g[0] = aersl7_1.p[0] / diag[0];
    i2 = nvarm2;
    for (i = 2; i <= i2; ++i) 
    {
/* L3: */
	g[i - 1] = (aersl7_1.p[i - 1] - sub[i - 1] * g[i - 2]) / diag[i - 1];
    }
    aersl7_1.p[nvarm2 - 1] = g[nvarm2 - 1];
    i2 = nvarm3;
    for (j = 1; j <= i2; ++j) 
    {
	i = nvarm2 - j;
/* L4: */
	aersl7_1.p[i - 1] = g[i - 1] - super[i - 1] * aersl7_1.p[i] / diag[i - 1];
    }
    return 0;
} /* setup_ */
