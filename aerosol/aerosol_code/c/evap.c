/* ----------------------------------------------------------------------- */
/* *     THIS IS AN EXAMPLE OF A EVAPORATION FUNCTION SUBROUTINE         * */
/* *    REFS. SKINNER,L.M. AND SAMBLES,J.R., 'THE KELVIN EQUATION - A    * */
/* *    REVIEW', J. AEROSOL SCIENCE, VOL. 3, PP. 199-210 (1972) AND      * */
/* *    FRIEDLANDER,S.K., 'SMOKE, DUST AND HAZE', P.229, JOHN WILEY AND  * */
/* *    SONS, NEW YORK (1977)                                            * */
/* *                                                                     * */
/* *                    VARIABLES                                        * */
/* *                                                                     * */
/* *    BOLTZ=BOLTZMANN'S CONSTANT  (ERG/DEG K)                          * */
/* *    COEF(I)=EVAPORATION COEFFICIENT  (OUTPUT ARRAY IN 1/SEC)         * */
/* *    D=PARTICLE DIAMETER  (CM)                                        * */
/* *    EQNUM=EQUILIBRIUM NUMBER CONCENTRATION  (1/CC)                   * */
/* *    NPTS=NUMBER OF POINTS  (INPUT)                                   * */
/* *    surten=SURFACE TENSION  (DYNES/CM)                               * */
/* *    TEMPK=TEMPERATURE  (DEG K)                                       * */
/* *    TIME=TIME  (INPUT IN SEC)                                        * */
/* *    V(I)=PARTICLE VOLUME  (INPUT ARRAY IN CC)                        * */
/* *    VMOLEC=MOLECULAR VOLUME  (CC)                                    * */
/* ----------------------------------------------------------------------- */

/* Subroutine */ 
int evap(int *npts, double *time, double *v, double *coef)
{
    /* Initialized data */

    static double surten = 70.;
    static double vmolec = 6.54e-23;
    static double boltz = 1.38e-16;
    static double tempk = 298.;
    static double eqnum = 1e6;

    /* System generated locals */
    int p_count;
    double d1;

    /* Builtin functions */
    double pow_dd(double *, double *), exp(double);

    /* Local variables */
    static double d;
    static int i;
    extern /* Subroutine */ int beta_(int *, double *, double *, double *, double *);

    /* Parameter adjustments */
    --coef;
    --v;

    /* Function Body */
    p_count = *npts;
    for (i = 1; i <= p_count; ++i) 
    {
	coef[i] = 0.o;
	if ((doublereal) aersl1_1.iparm[0] == 1.0) 
	{
	    //goto L1;
	}
	d1 = v[i] * 1.9099;
	d = pow_dd(&d1, &c_b574);
	d1 = v[i] - vmolec;
	beta_(&c__1, time, &vmolec, &d1, &coef[i]);
	coef[i] = coef[i] * eqnum * exp(surten * 4. * vmolec / (d * boltz * tempk));
//L1:
//	;
    }
    return 0;
} /* evap_ */

