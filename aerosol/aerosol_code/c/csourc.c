/* ----------------------------------------------------------------------- */
/* *    THIS IS AN EXAMPLE OF A SOURCE AND/OR REMOVAL SUBROUTINE IN THE  * */
/* *    CONTINUOUS REGIME                                                * */
/* *                                                                     * */
/* *                  VARIABLES                                          * */
/* *                                                                     * */
/* *    DISTV(I)=THE DISTRIBUTION FUNCTION OF PARTICLE VOLUME OF THE I-TH* */
/* *      PARTICLE SIZE  (INPUT ARRAY IN  1/CC/CC)                       * */
/* *    GLSIG=NATURAL LOGARITHM OF THE GEOMETRIC STANDARD DEVIATION      * */
/* *    GMEAN=GEOMETRIC MEAN VOLUME  (CC)                                * */
/* *    REMOVE(J,I)=THE J-TH REMOVAL TERM  FOR THE I-TH PARTICLE SIZE    * */
/* *                (OUTPUT ARRAY IN 1/CC/SEC/CC)                        * */
/* *    SOURCE(J,I)=THE J-TH SOURCE TERM  FOR THE I-TH PARTICLE SIZE     * */
/* *                (OUTPUT ARRAY IN 1/CC/SEC/CC)                        * */
/* *    TIME=TIME (INPUT IN SEC)                                         * */
/* *    TOT=NUMBER PER CC                                                * */
/* *    V(I)=I-TH PARTICLE VOLUME  (INPUT ARRAY IN CC)                   * */
/* ----------------------------------------------------------------------- */

/* Subroutine */ 
int csourc_(int *npts, double *time, double *v,	double *distv, double *source, double *remove)
{
    /* System generated locals */
    int p_count;
    double d1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static int i;

    /* Parameter adjustments */
    remove -= 4;
    source -= 4;
    --distv;
    --v;

    /* Function Body */
    p_count = *npts;
    for (i = 1; i <= p_count; ++i) 
    {
/* Computing 2nd power */
	d1 = log(v[i] / c_b565) / c_b566;
	source[i * 3 + 1] = .3989 * c_b564 * exp(-.5 * (d1 * d1)) / (v[i] * c_b566);
/* Computing 2nd power */
	d1 = log(v[i] / c_b568) / c_b569;
	source[i * 3 + 2] = .3989 * c_b564 * exp(-.5 * (d1 * d1)) / (v[i] * c_b569);
	remove[i * 3 + 1] = v[i] * -5e6 * exp(*time * -.005) * distv[i];
/* L1: */
	remove[i * 3 + 2] = v[i] * -6e4 * exp(*time * -.01) * distv[i];
    }
    return 0;
} /* csourc_ */

