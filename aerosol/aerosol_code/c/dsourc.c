/* ----------------------------------------------------------------------- */
/* *    THIS IS AN EXAMPLE OF A SOURCE AND/OR REMOVAL SUBROUTINE IN THE  * */
/* *    THE DISCRETE REGIME                                              * */
/* *                                                                     * */
/* *                    VARIABLES                                        * */
/* *                                                                     * */
/* *    CLSCON(I)=I-TH CLUSTER CONCENTRATION  (INPUT IN 1/CC)            * */
/* *    I=CLUSTER NUMBER  (INPUT)                                        * */
/* *    REMOVE(J,I)=THE J-TH REMOVAL TERM FOR THE I-TH CLUSTER (OUTPUT   * */
/* *      ARRAY IN 1/SEC/CC)                                             * */
/* *    SOURCE(J,I)=THE J-TH SOURCE TERM FOR THE I-TH CLUSTER  (OUTPUT   * */
/* *      ARRAY IN 1/SEC/CC)                                             * */
/* *    TIME=TIME  (INPUT IN SEC)                                        * */
/* ----------------------------------------------------------------------- */

/* Subroutine */ 
int dsourc_(int *npts, double *time, double *clscon, double *source, double *remove)
{
    /* System generated locals */
    int p_count;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static int i;

    /* Parameter adjustments */
    remove -= 4;
    source -= 4;
    --clscon;

    /* Function Body */
    p_count = *npts;
    for (i = 1; i <= p_count; ++i) 
    {
	remove[i * 3 + 1] = clscon[i] * - 0.01;
	remove[i * 3 + 2] = exp(-((doublereal) i)) * -.1 * clscon[i] * 	exp(*time * - 0.01);
	source[i * 3 + 1] = 0.0;
	source[i * 3 + 2] = exp(-((doublereal) i)) * 1e3;
	if ((doublereal) i > 1.) 
	{
		//goto L1;
	}
	source[i * 3 + 1] = 1e7;
//L1:
//	;
    }
    return 0;
} /* dsourc_ */

