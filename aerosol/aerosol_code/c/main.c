/* main.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct 
{
	double a[60];
	double x[60];
	double w[40];
} aersl2_;

#define aersl2_1 aersl2_

struct aersl5_1_ 
{
	int isave; 
	int icondn;
	int igrid; 
	int iprnt;
	int iptim; 
	int iread; 
	int isetsw;
	int isetup;
	int ispswt;
	int isplin;
	int isstat;
	int iswtch;
	int issflg;
	int itimdp;
	int ncremv;
	int ncsorc;
	int ndcls[20];
	int ndiscr;
	int ndremv;
	int ndsorc;
	int nmulp1;
	int nmxvar;
	int nprnt;
	int nptbet;
	int nquadp;
	int nvar;
	int nvart;
	int nvarp1;
	int nvarp2;
	int nvarm1;
	int maxevp;
	int minw;
};

#define aersl5_1 (*(struct aersl5_1_ *) &aersl5_)

union {
    struct {
	doublereal tout[16], dumtim, xa, xb, xboxa, algxba, relerr, abserr;
    } _1;
    struct {
	doublereal tout[16], time, xa, xb, xboxa, algxba, relerr, abserr;
    } _2;
    struct {
	doublereal tout[16], dumtim, xa, xb, xboxa, relerr, abserr;
    } _3;
} aersl6_;

#define aersl6_1 (aersl6_._1)
#define aersl6_2 (aersl6_._2)
#define aersl6_3 (aersl6_._3)

struct 
{
	int iparm[10];
} aersl1_;

#define aersl1_1 aersl1_

union {
    struct {
	doublereal z__[24], factor[24], zmass[24];
    } _1;
    struct {
	doublereal z__[24], factor[24], zvolm[24];
    } _2;
} aersl3_;

#define aersl3_1 (aersl3_._1)
#define aersl3_2 (aersl3_._2)

struct aersl8_1_ {
    doublereal wquad[12], wtquad[12];
};

#define aersl8_1 (*(struct aersl8_1_ *) &aersl8_)

struct {
    doublereal cndint, evpint;
    integer n1tilda;
} cnevint_;

#define cnevint_1 cnevint_

struct {
    doublereal besr[800]	/* was [20][40] */, betadc[800]	/* was [20][
	    40] */, betdsr[210], evppls[40], evpstr[60], fact1a[960]	/* 
	    was [24][40] */, fact1b[960]	/* was [24][40] */, fact2[
	    1440]	/* was [24][60] */, facte1[24], facte3[4];
    integer maxgdc[40];
    doublereal wd[800]	/* was [20][40] */, wmz1[960]	/* was [24][40] */, 
	    wpls[40], z1[960]	/* was [24][40] */, z2[960]	/* was [24][
	    40] */, ze[4];
} aersl4_;

#define aersl4_1 aersl4_

struct {
    doublereal p[38];
} aersl7_;

#define aersl7_1 aersl7_

/* Initialized data */

struct {
    integer fill_1[3];
    integer e_2;
    integer fill_3[1];
    integer e_4;
    integer fill_5[37];
    integer e_6;
    integer fill_7[7];
    } aersl5_ = { {0}, 2, {0}, 1, {0}, 24 };

struct 
{
    doublereal e_1[24];
} aersl8_= {.064056892862605626085, .191118867473616309159, 
	    .315042679696163374387, .433793507626045138487, 
	    .545421471388839535658, .648093651936975569252, 
	    .740124191578554364244, .820001985973902921954, 
	    .886415527004401034213, .938274552002732758524, 
	    .974728555971309498198, .99518721999702136018, 
	    .127938195346752156974, .125837456346828296121, 
	    .121670472927803391204, .115505668053725601353, 
	    .107444270115965634783, .09761865210411388827, 
	    .086190161531953275917, .073346481411080305734, 
	    .059298584915436780746, .044277438817419806169, 
	    .028531388628933663181, .012341229799987199547 };


/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n2 = -2;
static integer c_n1 = -1;
static doublereal c_b174 = 1.;
static doublereal c_b181 = .33333333333333;
static integer c__4 = 4;
static doublereal c_b564 = 100.;
static doublereal c_b565 = 4.08e-19;
static doublereal c_b566 = 1.22;
static doublereal c_b568 = 5e-18;
static doublereal c_b569 = 1.8;
static doublereal c_b574 = .33333333333;
static integer c__59 = 59;
static integer c_n78 = -78;
static integer c_n54 = -54;
static integer c__36 = 36;
static integer c__97 = 97;

/* * */
/* *    NOTE: IN THE WRITEUP THE I-TH CLUSTER MEANS A PARTICLE IN THE */
/* *          DISCRETE REGIME CONTAINING I MONOMERS.  THE I-TH MULTIPLET */
/* *          IS THE I-TH GRID POINT, WHICH IS SPECIFIED BY THE NUMBER OF */
/* *          MONOMERS CONTAINED IN THE PARTICLE. */
/* * */
/* *    R.J.=RIGHT JUSTIFIED INTEGER */
/* *    GDE=GENERAL DYNAMIC EQUATION */
/* *    DC GDE=DISCRETE-CONTINUOUS GENERAL DYNAMIC EQUATION */
/* * */
/* *                  FIRST CARD OF INPUT */
/* * */
/* *    IN COLUMNS 1-2 THE NUMBER OF CASES TO BE COMPUTED    R.J. */
/* *    IF A BLANK CARD IS USED ONE CASE WILL BE ASSUMED. */
/* * */
/* * */
/* *    FOR EACH CASE THE FOLLOWING SET OF CARDS MUST FOLLOW.  UNLESS */
/* *    SPECIFIED AS RIGHT JUSTIFIED (R.J.) INTEGER, ALL INPUT DATA IS TO */
/* *    BE IN REAL FORMAT. */
/* * */
/* *                  FIRST CARD/CASE */
/* * */
/* *    COLUMNS */
/* *     1-40    TITLE TO BE PRINTED AS A HEADING FOR THIS CASE */
/* *    41-50    SMALLEST PARTICLE SIZE  (IF THE DC GDE IS USED THIS IS */
/*                THE MONOMER SIZE. THE PROGRAM AUTOMATICALLY DETERMINES */
/* *               THE SMALLEST SIZE IN THE CONTINUOUS REGIME) */
/* *    51-60    LARGEST PARTICLE SIZE */
/* *    61-70    RELATIVE INTEGRATION ERROR */
/* *    71-80    ABSOLUTE INTEGRATION ERROR */
/* * */
/* *                  SECOND CARD/CASE */
/* * */
/* *     1- 2    NUMBER OF GRID POINTS   (10-40)             R.J. */
/* *     3- 4    MAXIMUM NUMBER OF GRID POINTS  (10-40)      R.J. */
/* *     5- 6    NUMBER OF DISCRETE SIZES (0,2-20)           R.J. */
/* *               (IF SET TO ZERO THE CONTINUOUS GDE IS ASSUMED) */
/* *     7- 8    NUMBER OF DISCRETE MULTIPLETS IN THE CONTINUOUS REGIME */
/* *               (0,2-20)                                  R.J. */
/* *               (MUST BE SET TO ZERO FOR THE CONTINUOUS GDE) */
/* *     9-10    NUMBER OF OUTPUTS (1-16)                    R.J. */
/* *    11-12    THE OUTPUT TIME NUMBER WHEN THE PROGRAM IS TO SWITCH */
/* *             INTERPOLATION METHODS.  FOR EXAMPLE, IF 5 OUTPUT TIMES */
/* *             ARE GIVEN AS 0.0,60.,300.,600., AND 1200. SECONDS, THEN */
/* *             A 2 WOULD RESULT IN THE PROGRAM SWITHCING FROM THE */
/* *             INITIALLY SPECIFIED INTERPOLATION METHOD TO THE OTHER */
/* *             ALTERNATIVE IMMEDIATELY AFTER 60. SECONDS OF SIMULATION. */
/* *             A NUMBER LESS THAN OR EQUAL TO 1 OR GREATER THAN THE */
/* *             NUMBER OF OUTPUTS WOULD RESULT IN NO CHANGE OF */
/* *             INTERPOLATING METHODS.                       R.J. */
/* *    13-14    THE OUTPUT TIME NUMBER WHEN THE PROGRAM IS TO */
/* *             SWITCH TO THE STEADY STATE APPROXIMATION FOR THE */
/* *             DISCRETE REGIME.  IF SET TO -1 THE PROGRAM WILL PRINT */
/* *             THE STEADY STATE CONCENTRATIONS IN ADDITION TO THE */
/* *             TRANSIENT CONCENTRATIONS.  IF SET TO -2 THE PROGRAM */
/* *             WILL SWITCH TO THE STEADY STATE APPROXIMATION WHEN */
/* *             THE SPECIFIED ERROR TOLERANCES ARE MET.  (SET TO */
/* *             ZERO FOR THE CONTINUOUS GDE)                 R.J. */
/* *    15-20            *** BLANK *** */
/* *    21-23    GRID POINT SELECTION */
/* *               0=PROGRAM WILL SELECT LOGARITHMICALLY SPACED GRID */
/* *                 IN PARTICLE SIZE */
/* *               FOR VALUES GREATER THAN OR EQUAL TO ONE, THIS IS THE */
/* *                 NUMBER OF GRID POINTS TO BE READ IN ASCENDING ORDER. */
/* *                 FOR NEGATIVE VALUES, READ IN DESCENDING ORDER THE */
/* *                 ABSOLUTE VALUE NUMBER OF GRID POINTS. THE PROGRAM */
/* *                 WILL SELECT THE REMAINING GRID POINTS AS */
/* *                 LOGARITHMICALLY SPACED IN PARTICLE SIZE.  R.J. */
/* *    24-30            *** BLANK *** */
/* *       31    NUMBER OF INTERPOLATING POINTS BETWEEN GRID POINTS */
/* *       32    NUMBER OF CONTINUOUS SOURCE TERMS  (0-3) */
/* *       33    NUMBER OF CONTINUOUS REMOVAL TERMS (0-3) */
/* *       34    NUMBER OF DISCRETE SOURCE TERMS    (0-3) */
/* *       35    NUMBER OF DISCRETE REMOVAL TERMS   (0-3) */
/* *       36    0=NO CONDENSATION, 1=WITH CONDENSATION */
/* *       37    INITIAL SPLINE INTERPOLATION METHOD NUMBER */
/* *               0=CUBIC SPLINES */
/* *               1=LINEAR SPLINES */
/* *               2=LOGARITHMIC SPLINES (POWER LAW FIT) */
/* *       38    CALL SETUP BEFORE INTERPOLATING  (0=NO/1=YES) */
/* *       39    SPLINE INTERPOLATION METHOD NUMBER FOR SWITHCED METHOD */
/* *       40    CALL SETUP FOR SWITCHED METHOD  (0=NO/1=YES) */
/* *       41    TIME DEPENDENT COAGULATION AND EVAPORATION */
/* *               COEFFICIENTS (0=NO/1=YES) */
/* *    42-50            *** BLANK *** */
/* *    51-60    USER DEFINED PARAMETER SWITCHES */
/* * */
/* *                  THIRD CARD(S)/CASE */
/* * */
/* *    IN FORMAT 8E10.4 THE TIMES FOR WHICH OUTPUT IS DESIRED IN */
/* *    ASCENDING ORDER. IF THE FIRST OUTPUT TIME IS GREATER */
/* *    THAN ZERO THE SEVENTH CARD(S)/CASE MUST BE INCLUDED.  THE NUMBER */
/* *    OF OUTPUT TIMES IS SPECIFIED IN COLUMNS 9-10 OF THE SECOND */
/* *    CARD/CASE */
/* * */
/* *                  FOURTH CARD(S)/CASE */
/* * */
/* *    IN FORMAT 16I5 THE DISCRETE MULTIPLET NUMBERS IN ASCENDING */
/* *    ORDER.  THE FIRST AND LAST GRID POINTS ARE AUTOMATICALLY FIXED AT */
/* *    THE SMALLEST AND LARGEST PARTICLE SIZE RESPECTIVELY, IN THE */
/* *    CONTINUOUS REGIME.  THEREFORE THE FIRST MULTIPLET NUMBER GIVEN */
/* *    MUST BE GREATER THAN OR EQUAL TO THE NUMBER OF DISCRETE SIZES PLUS */
/* *    TWO. THE MULTIPLET NUMBERS MUST RANGE FROM THE NUMBER OF DISCRETE */
/* *    SIZES PLUS TWO TO TWICE THE NUMBER OF DISCRETE SIZES. THE LAST */
/* *    MULTIPLET NUMBER MUST BE EQUAL TO TWICE THE NUMBER OF DISCRETE */
/* *    SIZES. THE NUMBER OF MULTIPLETS GIVEN SHOULD CORRESPOND TO */
/* *    COLUMNS 7-8 IN THE SECOND CARD/CASE */
/* *    (IF NO DISCRETE SIZES ARE USED DO NOT PUT IN ANY CARD) */
/* * */
/* *                  FIFTH CARD(S)/CASE */
/* * */
/* *    IN FORMAT 8E10.4 THE INITIAL NUMBER CONCENTRATION IN THE */
/* *    DISCRETE REGIME IN ASCENDING PARTICLE SIZE ORDER. IF THE */
/* *    STEADY STATE APPROXIMATION IS REQUESTED FOR THE FIRST */
/* *    OUTPUT TIME (I.E. COLUMNS 13-14 IN THE SECOND CARD/CASE */
/* *    IS 01), THEN THESE VALUES WILL BE USED AS THE INITIAL */
/* *    GUESS FOR DETERMINING THE STEADY STATE CONCENTRATIONS */
/* *    (IF NO DISCRETE SIZES ARE USED DO NOT PUT IN ANY CARD) */
/* * */
/* *                  SIXTH CARD(S)/CASE */
/* * */
/* *    IN FORMAT 8E10.4 THE PARTICLE SIZES OF THE GRID POINTS IN */
/* *    ASCENDING OR DESCENDING ORDER AS SPECIFIED IN CARD COLUMNS 21-23 */
/* *    IN THE SECOND CARD/CASE. */
/* *    THE NUMBER OF GRID POINTS GIVEN MUST BE LESS THAN THE NUMBER OF */
/* *    GRID POINTS MINUS THE NUMBER OF MULTIPLETS GIVEN MINUS TWO. */
/* *    FOR THE CONTINUOUS GDE THE FIRST GRID POINT IS AUTOMATICALLY SET */
/* *    AT THE SMALLEST PARTICLE SIZE.  THEREFORE GRID POINTS ARE READ */
/* *    IN STARTING WITH THE SECOND ACTUAL GRID POINT. FOR THE DISCRETE - */
/* *    CONTINUOUS GDE GRID POINTS ARE READ IN STARTING AFTER THE LAST */
/* *    MULTIPLET SIZE AND SHOULD THEREFORE HAVE A SIZE GREATER THAN */
/* *    THE MONOMER SIZE TIMES TWICE THE NUMBER OF DISCRETE SIZES. IF */
/* *    GRID POINTS ARE READ IN DECSENDING ORDER THE FIRST GRID POINT */
/* *    SIZE READ IN MUST BE LESS THAN THE LARGEST SIZE SPECIFIED IN */
/* *    COLUMNS 51-60 OF THE FIRST CARD/CASE. */
/* *    (NOT TO BE USED IF LOGARITHMICALLY SPACED GRID POINTS ARE */
/* *     REQUESTED, I.E. COLUMNS 21-23 IN THE SECOND CARD/CASE IS ZERO) */
/* * */
/* *                    SEVENTH CARD(S)/CASE */
/* * */
/* *    IN FORMAT 8E10.4 THE DISTRIBUTION VALUES FOR INITIAL TIMES */
/* *    GREATER THAN ZERO.  EVEN IF GRID POINTS ARE READ IN */
/* *    DESCENDING ORDER, THE DISTRIBUTION VALUES ARE READ IS */
/* *    ASCENDING ORDER. THIS CARD(S) IS REQUIRED ONLY IF THE FIRST */
/* *    OUTPUT TIME IS GREATER THAN ZERO */
/* * */
/* *    TO TRANSFER PARAMETER VALUES TO THE ROUTINES */
/* *    INSERT THE FOLLOWING: COMMON/AERSL1/IPARM(10) */
/* *    IF YOU ENCOUNTER ANY TROUBLE USING THE PROGRAM CONTACT */
/* *                      FRED GELBARD */
/* *                      DEPT. OF CHEMICAL ENGINEERING */
/* *                      CALTECH 208-41 */
/* *                      PASADENA,C.A.  91125 */
/* *                      TEL. 213-7956811 EXT 1671 OR EXT 1115 */
/* *---------------------------------------------------------------------- */
/* * */
/* *          VARIABLES USED BY THE PROGRAM */
/* * */
/* *    A(I)=FOR I=1 TO NVAR IT IS THE VALUE OF THE CONTINUOUS */
/* *      DISTRIBUTION AT THE I-TH GRID POINT.  FOR I=NVAR+1 TO NVART */
/* *      IT IS THE NUMBER OF PARTICLES AT THE DISCRETE SIZE I-NVAR */
/* *    ABSERR=ABSOLUTE INTEGRATION ERROR (SEE ODE CODE WRITE UP) */
/* *    ADIF(I)=THIS ARRAY IS USED TO TRANSFER THE A ARRAY TO THE ODE */
/* *      PACKAGE. */
/* *    ASAVE(I)=SAVED A ARRAY FROM LAST OUTPUT TIME */
/* *    ALGXBA=LN(XB/XA) */
/* *    ATRY(I)=NEXT ITERATION FOR THE CONCENTRATION OF THE I-TH CLUSTER, */
/* *      USING THE STEADY STATE APPROXIMATION */
/* *    BESR(I,J)=COAGULATION COEFFICIENT TIMES THE JACOBIAN FOR THE I-TH */
/* *      DISCRETE CLUSTER WITH A PARTICLE IN THE CONTINUOUS REGIME SUCH */
/* *      THAT A PARTICLE AT THE J-TH GRID POINT IS FORMED */
/* *    BETADC(I,J)=THE COAGULATION COEFFICIENT OF THE I-TH CLUSTER WITH */
/* *      THE J-TH GRID POINT */
/* *    BETDSR(I)=STORAGE ARRAY FOR CLUSTER-CLUSTER COAGULATION */
/* *      COEFFICIENTS */
/* *    BOTTOM=LOWER LIMIT OF INTEGRATION */
/* *    DADT(I)=THE DERIVATIVE OF ADIF(I) WITH RESPECT TO TIME */
/* *    DIAG(I)=DIAGONAL ELEMENTS OF TRI-DIAGONAL MATRIX */
/* *    DISTD=NUMBER DISTRIBUTION (NUMBER/CC) */
/* *    DISTS=SURFACE AREA DISTRIBUTION  (SQUARE CM/CC) */
/* *    DISTV=VOLUME DISTRIBUTION  (CC/CC) */
/* *    DMEAN=MEAN DIAMETER  (CM) */
/* *    EVAPR=TOTAL MONOMER EVAPORATION RATE FROM THE CONTINUOUS REGIME */
/* *    EVPSTR(I)=EVAPORATION COEFFICIENT AT THE I-TH GRID POINT. FOR */
/* *      I=NVARP2 TO NVART, THIS IS THE EVAPORATION COEFFICIENT OF THE */
/* *      I-TH CLUSTER */
/* *    FACT1A(I,J)=THE I-TH QUADRATURE FACTOR FOR THE FIRST TERM OF THE */
/* *      FIRST COAGULATION INTEGRAL FOR THE J-TH GRID POINT */
/* *    FACT1B(I,J)=THE I-TH QUADRATURE FACTOR FOR THE SECOND TERM OF THE */
/* *      FIRST COAGULATION INTEGRAL FOR THE J-TH GRID POINT */
/* *    FACT2(I,J)=THE I-TH QUADRATURE FACTOR FOR THE SECOND COAGULATION */
/* *      INTEGRAL FOR THE J-TH GRID POINT. FOR J=NVARP1 TO NVART THIS IS */
/* *      THE QUADRATURE FACTOR FOR THE SCAVENGING OF THE (J-NVAR)-TH */
/* *      CLUSTER BY PARTICLES IN THE CONTINUOUS REGIME */
/* *    FACTE1(I)=THE I-TH QUADRATURE FACTOR FOR THE GAIN OF MONOMER DUE */
/* *      TO EVAPORATION FROM THE CONTINUOUS REGIME */
/* *    FACTE3(I)=THE I-TH QUADRATURE FACTOR FOR EVAPORATION FROM THE */
/* *      CONTINUOUS REGIME TO FORM THE LARGEST DISCRETE PARTICLE */
/* *    FACTOR(I)=THE I-TH QUADRATURE WEIGHT FACTOR FOR INTEGRATION OVER */
/* *      THE ENTIRE CONTINUOUS DOMAIN */
/* *    HEADNG=CHARACTER STRING FOR HEADING */
/* *    ICASE=THE CASE NUMBER BEING SOLVED */
/* *    ICONDN=  0=NO CONDNESATION/1=WITH CONDENSATION */
/* *    IFLAG=FLAG TO ODE, 1 INDICATES THIS IS THE FIRST CALL */
/* *    IGRID=NUMBER OF GRID POINTS GIVEN */
/* *    IPARM(I)=THE I-TH USER SUPPLIED PARAMETER */
/* *    IPRNT=PRINTER LOGICAL UNIT NUMBER */
/* *    IREAD=CARD READER LOGICAL UNIT NUMBER */
/* *    ISETSW=OUTPUT NUMBER TO SWITCH INTERPOLATION METHOD */
/* *    ISETUP=FLAG TO CALL SETUP BEFORE INTERPOLATING (0=NO/1=YES) */
/* *    ISPLIN=FLAG TO INDICATE  0=CUBIC SPLINE INTERPOLATION */
/* *                             1=LINEAR SPLINE INTERPOLATION */
/* *                             2=LOGARITHMIC INTERPOLATION */
/* *    ISPSWT=FLAG TO CALL SETUP FOR SWITCHED INTERPOLATION METHOD */
/* *      (0=NO/1=YES) */
/* *    ISSFLG=STEADY STATE FLAG, 0=DO NOT USE THE STEADY STATE */
/* *      APPROXIMATION, 1=USE THE STEADY STATE APPROXIMATION FOR THE */
/* *      DISCRETE REGIME */
/* *    ISSTAT=OUTPUT NUMBER AT WHICH THE STEADY STATE APPROXIMATION IS */
/* *      TO BE USED.  IF ZERO DO NOT USE THE S.S. APPROXIMATION.  IF */
/* *      -1 PRINT S.S. BUT DO NOT USE.  IF -2 SWITCH TO S.S. WHEN */
/* *      SPECIFIED ERROR TOLERANCES HAVE BEEN MET */
/* *    ISTOP=ERROR FLAG, 1=ERROR,0=NO ERROR,-1=CORRECTABLE ERROR,CODE */
/* *            WILL REVERT TO PREVIOUS OUTPUT TIME,-2=PRINT ONLY DISCRETE */
/* *            REGIME AS DETERMINED FROM THE STEADY STATE APPROXIMATION */
/* *    ISWTCH=THE OUTPUT NUMBER FOR WHICH THE PROGRAM IS TO */
/* *      SWITCH INTERPOLATION METHODS */
/* *    ITIMDP=FLAG FOR TIME DEPENDENT COAGULATION AND EVAPORATION */
/* *      COEFFICIENTS (0=TIME INDEPENDENT/1=TIME DEPENDENT) */
/* *    JPRTIM=THE PRINTING DO-LOOP COUNTER */
/* *    MAXEVP=THE LAST GRID POINT IN THE CONTINUOUS REGIME FOR WHICH */
/* *      EVAPORATION IS STILL INCLUDED */
/* *    MAXGDC(I)=MAXIMUM DISCRETE CLUSTER SIZE THAT CAN COAGULATE WITH A */
/* *      PARTICLE IN THE CONTINUOUS REGIME TO FORM A PARTICLE AT THE I-TH */
/* *      GRID POINT */
/* *    MAXSZ=THE LARGEST DISCRETE CLUSTER THAT CAN COAGULATE WITH ANOTHER */
/* *      DISCRETE CLUSTER TO FROM THE CLUSTER OF INTEREST */
/* *    MINW=SMALLEST GRID POINT FOR WHICH THE FIRST COAGULATION INTEGRAL */
/* *      IS TO BE COMPUTED */
/* *    MXCASE=MAXIMUM NUMBER OF CASES TO SOLVE. WILL DEFAULT TO ONE IF */
/* *      LESS THAN OR EQUAL TO ZERO */
/* *    NCREMV=NUMBER OF CONTINUOUS REMOVALS */
/* *    NCSORC=NUMBER OF CONTINUOUS SOURCES */
/* *    NDCLS(I)=THE DISCRETE SIZE NUMBER OF THE I-TH GRID POINT GIVEN IN */
/* *      THE CONTINUOUS REGIME */
/* *    NDISCR=NUMBER OF DISCRETE SIZES */
/* *    NDREMV=NUMBER OF DISCRETE REMOVALS */
/* *    NDSORC=NUMBER OF DISCRETE SOURCES */
/* *    NMULGV=NUMBER OF DISCRETE MULTIPLETS GIVEN IN THE CONTINUOUS */
/* *      REGIME */
/* *    NMULP1=NUMBER OF MULTIPLETS GIVEN PLUS ONE */
/* *    NMULP2=NUMBER OF MULTIPLETS GIVEN PLUS TWO */
/* *    NMXVAR=MAXIMUM NUMBER OF GRID POINTS */
/* *    NPRNT=NUMBER OF OUTPUTS */
/* *    NPTBET=NUMBER OF INTERPOLATED POINTS BETWEEN GRID POINTS */
/* *    NQUADP=NUMBER OF QUADRATURE POINTS (MUST BE AN EVEN NUMBER) */
/* *    NQUAD2=HALF THE NUMBER OF QUADRATURE POINTS */
/* *    NTROUB=NUMBER OF POINTS AT WHICH DISTRIBUTION IS NEGATIVE */
/* *    NVAR=NUMBER OF VARIABLES IN THE CONTINUOUS REGIME */
/* *          (EQUAL TO THE NUMBER OF GRID POINTS) */
/* *    NVARP1=NVAR+1=NUMBER OF VARIABLES PLUS ONE */
/* *    NVARP2=NVAR+2=NUMBER OF VARIABLES PLUS TWO */
/* *    NVART=TOTAL NUMBER OF VARIABLES=NVAR+NDISCR */
/* *    P(I)=STORAGE ARRAY FOR THE SECOND DERIVATIVE OF M(W,T) WITH */
/* *      RESPECT TO W AT THE I+1 GRID POINT.  AT THE FIRST AND LAST GRID */
/* *      POINT THE SECOND DERIVATIVE IS ASSUMED TO BE ZERO. */
/* *    RELERR=RELATIVE ERROR (SEE ODE CODE WRITE UP) */
/* *    SCAVNG(I)=SCAVENGING RATE OF THE I-TH CLUSTER */
/* *    TOTALQ=ARRAY CONTAINING THE FIRST MOMENTS OF AEROSOL ADDED */
/* *      TO AND REMOVED FROM THE CONTINUOUS REGIME, AND ADDED TO */
/* *      AND REMOVED FROM THE DISCRETE REGIME FOR EACH MECHANISM, */
/* *      RESPECTIVELY */
/* *    TOTNUM=TOTAL NUMBER OF PARTICLES/CC */
/* *    TOTSUR=TOTAL SURFACE AREA  (SQUARE CM/CC) */
/* *    TOTVOL=TOTAL VOLUME OF PARTICLES  (CC/CC) */
/* *    TOUT(I)=THE I-TH OUTPUT TIME */
/* *    SUB(I)=SUB-DIAGONAL ELEMENTS OF TRI-DIAGONAL MATRIX */
/* *    SUPER(I)=SUPER-DIAGONAL ELEMENTS OF TRI-DIAGONAL MATRIX */
/* *    TXA=LN(2)/LN(XB/XA)= THE W SIZE OF A PARTICLE PRODUCED WHEN */
/* *      TWO OF THE SMALLEST PATICLES OF SIZE W=0 COAGULATE */
/* *    UPPER=UPPER LIMIT OF INTEGRATION */
/* *    W(I)=LN(X/XA)/LN(XB/XA) OF THE I-TH GRID POINT */
/* *    WBET(I)=THE I-TH W VALUE INTERPOLATED BETWEEN TWO GRID POINTS */
/* *    WD(I,J)=THE W OF A PARTICLE IN THE CONTINUOUS REGIME, SUCH THAT */
/* *      WHEN IT COAGUALTES WITH THE I-TH DISCRETE CLUSTER A PARTICLE OF */
/* *      SIZE W(J) IS FORMED */
/* *    WMZ1(I,J)=THE W SIZE OF A PARTICLE WHICH WHEN IT COAGULATES WITH A */
/* *      PARTICLE AT THE I-TH QUADRATURE POINT, A PARTICLE AT THE J-TH */
/* *      GRID POINT IS FORMED */
/* *    WPLS(I)=THE W SIZE OF THE I-TH GRID POINT PLUS A MONOMER */
/* *    X(I)=PARTICLE SIZE AT THE I-TH GRID POINT */
/* *    XA=SIZE OF SMALLEST PARTICLE */
/* *    XB=SIZE OF LARGEST PARTICLE */
/* *    XBET=PARTICLE SIZE OF INTERPOLATED POINTS BETWEEN TWO GRID POINTS */
/* *    XBOXA=XB/XA */
/* *    XTROUB(I)=PARTICLE SIZE ARRAY AT WHICH DISTRIBUTION IS NEGATIVE */
/* *    Z(I)=THE W SIZE OF THE I-TH QUADRATURE POINT OVER THE ENTIRE */
/* *      CONTINUOUS REGIME */
/* *    ZMASS(I)=THE PARTICLE SIZE AT THE I-TH QUADRATURE POINT */
/* *    Z1(I,J)=THE I-TH QUADRATURE POINT FOR THE FIRST COAGULATION */
/* *      INTEGRAL FOR THE J-TH GRID POINT */
/* *    Z2(I,J)=THE I-TH QUADRATURE POINT FOR THE SECOND COAGUALTION */
/* *      INTEGRAL FOR THE J-TH GRID POINT */
/* *---------------------------------------------------------------------- */
/*     PROGRAM AEROSL */
/* ----------------------------------------------------------------------- */
/* *    THIS IS THE MAIN PROGRAM FOR THE AEROSOL CODE                    * */

/* ----------------------------------------------------------------------- */
int MAIN__()
{
    /* Format strings */
    static char fmt_30[] = "(\002 ENTER INPUT FILE: \002)";
    static char fmt_31[] = "(a40)";
    static char fmt_32[] = "(\002 ENTER OUTPUT FILE: \002)";
    static char fmt_1[] = "(i2)";
    static char fmt_2[] = "(////52x,\002NUMBER OF CASES =\002,i3)";
    static char fmt_20[] = "(//47x,\002NO DATA INPUT ERRORS FOUND\002)";
    static char fmt_15[] = "(\0020*** ERROR TOLERANCES TOO SMALL, WILL BE IN\
CREASED, TIME=\002,e12.4)";
    static char fmt_12[] = "(\0020*** SOLUTION VANISHES, WILL USE ABSOLUTE E\
RROR OF 12**-30, TIME=\002,e10.4)";
    static char fmt_4[] = "(\0020*** SOMETHING IS WRONG, MORE THAN 500 STEPS\
 NEEDED TO REACH TIME=\002,e10.4,\002, RECHECK INPUT OR CONTACT FRED GELBAR\
D\002)";
    static char fmt_11[] = "(\0020*** EQUATIONS SEEM TO BE STIFF TIME REACHE\
D= \002,e12.4)";
    static char fmt_19[] = "(\0020*** INVALID INPUT TO ODE PACKAGE\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    olist o__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(), f_open(olist *);
    int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;
    extern int ode_(U_fp,integer *,doublereal *,doublereal *,doublereal *,doublereal *,doublereal *,integer *,doublereal *,integer *);
    static doublereal time, work[1360];
    extern int check_(integer *);
    static integer icase;
    static doublereal asave[20];
    static integer iwork[5], istop;
    extern int input_(), setup_(integer *, integer *, 
	    integer *), coagcf_(integer *), chktbl_(integer *);
    extern int diffun_(...);
    static integer maxcas, iodepk;
    extern int steady_(doublereal *);
    static integer jprtim, jstart;
    static char filnam1[40], filnam2[40];
    extern int output_(integer *);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_30, 0 };
    static cilist io___2 = { 0, 5, 0, fmt_31, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_32, 0 };
    static cilist io___5 = { 0, 5, 0, fmt_31, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_2, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_15, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_12, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_4, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_19, 0 };



    s_wsfe(&io___1);
    e_wsfe();
    s_rsfe(&io___2);
    do_fio(&c__1, filnam1, (ftnlen)40);
    e_rsfe();
    s_wsfe(&io___4);
    e_wsfe();
    s_rsfe(&io___5);
    do_fio(&c__1, filnam2, (ftnlen)40);
    e_rsfe();
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = 40;
    o__1.ofnm = filnam1;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = 40;
    o__1.ofnm = filnam2;
    o__1.orl = 0;
    o__1.osta = "NEW";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    io___7.ciunit = aersl5_1.iread;
    s_rsfe(&io___7);
    do_fio(&c__1, (char *)&maxcas, (ftnlen)sizeof(integer));
    e_rsfe();
    if (maxcas <= 0) {
	maxcas = 1;
    }
    io___9.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___9);
    do_fio(&c__1, (char *)&maxcas, (ftnlen)sizeof(integer));
    e_wsfe();

/*        IN THE FOLLOWING DO LOOP TERMINATING WITH STATEMENT NUMBER 10 */
/*        COMPUTE EACH CASE */

    i__1 = maxcas;
    for (icase = 1; icase <= i__1; ++icase) {

/*        READ AND CHECK INPUT DATA. IF AN INPUT ERROR IS FOUND (I.E. */
/*        ISTOP=1), SKIP THE CASE. IF NO ERRORS ARE FOUND (I.E. ISTOP=0) */
/*        CALL THE ROUTINE SETUP IF NEEDED FOR INTERPOLATION AND CONTINUE */

	input_();
	iodepk = 0;
	time = aersl6_1.tout[0];
	aersl6_1.dumtim = time;
	istop = 0;
	aersl5_1.iptim = 1;
	if (aersl5_1.isetup == 1) {
	    setup_(&aersl5_1.isplin, &c__1, &aersl5_1.nvar);
	}
	check_(&istop);
	if (istop == 1) {
	    goto L10;
	}
	io___14.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___14);
	e_wsfe();
	if (aersl5_1.isetup == 1) {
	    setup_(&aersl5_1.isplin, &c__1, &aersl5_1.nvar);
	}

/*        CHECK FOR NEGATIVE VALUES OF THE DISTRIBUTION.  IF A NEGATIVE */
/*        VALUE WAS FOUND AND GRID POINTS WERE ADDED (I.E. ISTOP=-1), */
/*        RECHECK THE INTERPOLATION SCHEME.  IF A NEGATIVE VALUE WAS */
/*        FOUND AND NO GRID POINTS WERE ADDED (I.E. ISTOP=1) SKIP THE */
/*        CASE */

L3:
	chktbl_(&istop);
	if (istop == -1) {
	    goto L3;
	}
	if (istop == 1) {
	    goto L10;
	}

/*        NOW THAT THE INPUT DATA AND DISTRIBUTION HAVE BEEN CHECKED, */
/*        OUTPUT THE INITIAL CONDITIONS AND INITIALIZE OUTPUT COUNTERS */
/*        (I.E. IPTIM AND JPRTIM), AND THE FLAG FOR THE O.D.E. PACKAGE */
/*        (I.E. IODEPK).  TWO VARIABLES ARE USED FOR THE OUTPUT COUNTERS */
/*        SO AS NOT TO PASS THE DO-LOOP COUNTER IN THE COMMON BLOCK */
/*        /AERSL5/, TO AVOID THE POSSIBLITY OF THE USER CHANGING */
/*        THE COUNTER IN HIS ROUTINE */

	coagcf_(&c__0);
	output_(&c__0);
	if (aersl5_1.isstat == 1) {
	    aersl5_1.issflg = 1;
	}
	if (aersl5_1.issflg == 1) {
	    steady_(&aersl6_1.dumtim);
	}
	if (aersl5_1.issflg == 1) {
	    output_(&c_n2);
	}
	if (aersl5_1.isstat >= 0 || aersl5_1.isstat == -2) {
	    goto L22;
	}
	i__2 = aersl5_1.ndiscr;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L26: */
	    asave[i__ - 1] = aersl2_1.a[aersl5_1.nvar + i__ - 1];
	}
	steady_(&aersl6_1.dumtim);
	output_(&c_n2);
	i__2 = aersl5_1.ndiscr;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L27: */
	    aersl2_1.a[aersl5_1.nvar + i__ - 1] = asave[i__ - 1];
	}
L22:
	if (aersl5_1.nprnt == 1) {
	    goto L10;
	}

	jprtim = 2;
	aersl5_1.iptim = jprtim;
L6:
	jstart = jprtim;
	iodepk = 1;

/*        IN THE FOLLOWING DO LOOP TERMINATING WITH STATEMENT NUMBER 5 */
/*        COMPUTE THE SOLUTION OF THE CASE FOR NPRNT OUTPUT TIMES */

	i__2 = aersl5_1.nprnt;
	for (jprtim = jstart; jprtim <= i__2; ++jprtim) {
	    aersl5_1.iptim = jprtim;
L17:
	    i__3 = aersl5_1.nvar + aersl5_1.ndiscr * (1 - aersl5_1.issflg);
	    ode_((U_fp)diffun_, &i__3, aersl2_1.a, &time, &aersl6_1.tout[
		    jprtim - 1], &aersl6_1.relerr, &aersl6_1.abserr, &iodepk, 
		    work, iwork);
	    aersl6_1.dumtim = time;
	    switch (iodepk) {
		case 1:  goto L13;
		case 2:  goto L13;
		case 3:  goto L9;
		case 4:  goto L7;
		case 5:  goto L8;
		case 6:  goto L14;
		case 7:  goto L18;
	    }
L9:
	    io___21.ciunit = aersl5_1.iprnt;
	    s_wsfe(&io___21);
	    do_fio(&c__1, (char *)&time, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    goto L17;
L14:
	    io___22.ciunit = aersl5_1.iprnt;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&time, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    aersl6_1.abserr = 1e-30;
	    goto L17;
L13:
	    if (aersl5_1.isetup == 1) {
		setup_(&aersl5_1.isplin, &c__0, &aersl5_1.nvar);
	    }
	    if (aersl5_1.issflg == 1 && aersl5_1.itimdp == 1) {
		coagcf_(&c__0);
	    }
	    if (aersl5_1.issflg == 1) {
		steady_(&aersl6_1.dumtim);
	    }

/*        IF GRID POINTS HAVE BEEN ADDED (I.E. ISTOP=-1), RESTART */
/*        CALCULATION AND TRY AGAIN TO REACH TOUT(JPRTIM). IF THE */
/*        DISTRIBUTION IS NEGATIVE AND NO GRID POINTS HAVE BEEN ADDED */
/*        (I.E. ISTOP=1), SKIP TO THE NEXT CASE. */

	    chktbl_(&istop);
	    if (istop == 1) {
		goto L10;
	    }
	    if (istop == 0) {
		goto L21;
	    }
	    time = aersl6_1.tout[jprtim - 2];
	    aersl6_1.dumtim = time;
	    i__3 = istop + aersl5_1.itimdp;
	    coagcf_(&i__3);
	    if (aersl5_1.issflg == 1) {
		steady_(&aersl6_1.dumtim);
	    }
	    output_(&istop);
	    if (jprtim == aersl5_1.isstat + 1) {
		aersl5_1.issflg = 0;
	    }
	    goto L6;
L21:
	    output_(&istop);

/*        CHECK IF THE SPLINE INTERPOLATION METHOD IS TO BE CHANGED OR */
/*        IF THE STEADY STATE APPROXIMATION IS TO BE USED FOR THE */
/*        DISCRETE REGIME */

	    if (aersl5_1.isstat >= 0 || aersl5_1.issflg == 1) {
		goto L23;
	    }
	    i__3 = aersl5_1.ndiscr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L24: */
		asave[i__ - 1] = aersl2_1.a[aersl5_1.nvar + i__ - 1];
	    }
	    if (aersl5_1.itimdp == 1) {
		coagcf_(&c__0);
	    }
	    steady_(&aersl6_1.dumtim);
	    output_(&c_n2);
	    if (aersl5_1.isstat == -1) {
		goto L28;
	    }
	    i__3 = aersl5_1.nvart;
	    for (i__ = aersl5_1.nvarp1; i__ <= i__3; ++i__) {
		if ((d__1 = aersl2_1.a[i__ - 1] - asave[i__ - aersl5_1.nvar - 
			1], abs(d__1)) / aersl6_1.relerr > aersl2_1.a[i__ - 1]
			 + aersl6_1.abserr) {
		    goto L28;
		}
/* L29: */
	    }
	    aersl5_1.issflg = 1;
	    iodepk = 1;
	    goto L23;
L28:
	    i__3 = aersl5_1.ndiscr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L25: */
		aersl2_1.a[aersl5_1.nvar + i__ - 1] = asave[i__ - 1];
	    }

L23:
	    if (jprtim != aersl5_1.iswtch) {
		goto L16;
	    }
	    aersl5_1.isplin = aersl5_1.ispswt;
	    aersl5_1.isetup = aersl5_1.isetsw;
	    output_(&c_n1);

L16:
	    if (jprtim != aersl5_1.isstat) {
		goto L5;
	    }
	    iodepk = 1;
	    aersl5_1.issflg = 1;
	    if (aersl5_1.itimdp == 1) {
		coagcf_(&c__0);
	    }
	    steady_(&aersl6_1.dumtim);
	    output_(&c_n2);
L5:
	    ;
	}
	goto L10;

/*        PRINT ERROR MESSAGE AND PROCEED TO NEXT CASE */

L7:
	io___23.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___23);
	do_fio(&c__1, (char *)&time, (ftnlen)sizeof(doublereal));
	e_wsfe();
	goto L10;
L8:
	io___24.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___24);
	e_wsfe();
	goto L10;
L18:
	io___25.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___25);
	e_wsfe();
L10:
	;
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

int input_()
{
    /* Format strings */
    static char fmt_1[] = "(10a4,4e10.4/7i2,6x,i3,7x,11i1,9x,10i1)";
    static char fmt_3[] = "(\0021\002,43x,10a4/\0020SMALLEST PARTICLE VOLUM\
E\002,e33.4,11x,\002LARGEST PARTICLE VOLUME\002,e28.4/\0020RELATIVE INTEGRAT\
ION ERROR\002,e29.4,11x,\002ABSOLUTE INTEGRATION ERROR\002,e23.4/\0020TOTAL \
NUMBER OF GRID POINTS\002,i26,13x,\002MAXIMUM NUMBER OF GRID POINTS\002,i18,/\
\0020NUMBER OF DISCRETE SIZES\002,i29,13x,\002NUMBER OF MULTIPLETS GIVEN\002\
,i21,/\0020NUMBER OF OUTPUT TIMES\002,i31,13x,\002SWITCH INTERPOLATION AT OU\
TPUT NUMBER\002,i10/\0020STEADY STATE AT OUTPUT NUMBER\002,i24,13x,\002NUMBE\
R OF GRID POINTS GIVEN\002,i20/\0020NUMBER OF INTERPOLATING POINTS\002,i23,1\
3x,\002NUMBER OF CONTINUOUS SOURCE TERMS\002,i14)";
    static char fmt_4[] = "(\0020NUMBER OF CONTINUOUS REMOVAL TERMS\002,i19,\
13x,\002NUMBER OF DISCRETE SOURCE TERMS\002,i16/\0020NUMBER OF DISCRETE REMO\
VAL TERMS\002,i21,13x,\002CONDENSATION (0=NO/1=YES)\002,i22/\0020INITIAL SPL\
INE INTERPOLATION METHOD NUMBER\002,i11,13x,\002CALL SETUP BEFORE INTERPOLAT\
ING (0=NO/1=YES)\002,i3/\0020SPLINE NUMBER FOR SWITCHED METHOD\002,i20,13x\
,\002CALL SETUP FOR SWITCHED METHOD (0=NO/1=YES)\002,i4/\0020TIME DEPENDENT \
COAG./EVAP. COEFFICIENT (0=NO/1=YES)\002,i2,13x,\002USER SUPPLIED PARAMETER \
SWITCHES\002,5x,10i2)";
    static char fmt_5[] = "(8e10.4)";
    static char fmt_23[] = "(/55x,\002OUTPUT TIMES\002/\0020NO.\002,7x,\002T\
IME\002,6x,\002NO.\002,7x,\002TIME\002,6x,\002NO.\002,7x,\002TIME\002,6x,\
\002NO.\002,7x,\002TIME\002,6x,\002NO.\002,7x,\002TIME\002,6x,\002N0.\002,7x,\
\002TIME\002/6(i3,e14.4,3x))";
    static char fmt_28[] = "(16i5)";
    static char fmt_40[] = "(/47x,\002GRID POINTS AT THE FOLLOWING MULTIPL\
ETS\002,(/12i8))";
    static char fmt_12[] = "(//55x,\002GRID POINTS\002/\0020NO.\002,7x,\002S\
IZE\002,6x,\002NO.\002,7x,\002SIZE\002,6x,\002NO.\002,7x,\002SIZE\002,6x,\
\002NO.\002,7x,\002SIZE\002,6x,\002NO.\002,7x,\002SIZE\002,6x,\002N0.\002,7x,\
\002SIZE\002/6(i3,e14.4,3x))";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(), 
	    s_wsfe(cilist *), e_wsfe();
    double log(doublereal), pow_dd(doublereal *, doublereal *), d_sign(
	    doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, ii;
    extern /* Subroutine */ int xintl_(integer *, doublereal *, doublereal *);
    static integer igdpm1, igdpm2, nquad2, nvgdm1, nmulp2;
    static doublereal headng[10];
    static integer nmulgv;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 0, 0, fmt_1, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_3, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_4, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_5, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_23, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_28, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_5, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_5, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_5, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_12, 0 };


/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE READS IN THE INPUT DATA, PRINTS THE INPUT DATA,  * */
/* *    DETERMINES INITIAL GRID POINT LOCATIONS, AND THE VALUE OF THE    * */
/* *    INITIAL DISTRIBUTION AT THE GRID POINTS                          * */
/* ----------------------------------------------------------------------- */

/*        READ AND ECHO PRINT INPUT DATA */

    io___26.ciunit = aersl5_1.iread;
    s_rsfe(&io___26);
    for (i__ = 1; i__ <= 10; ++i__) {
	do_fio(&c__1, (char *)&headng[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    do_fio(&c__1, (char *)&aersl6_2.xa, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.xb, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.relerr, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.abserr, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl5_1.nvar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nmxvar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndiscr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nmulgv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nprnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.iswtch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isstat, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.igrid, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nptbet, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ncsorc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ncremv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndsorc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndremv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.icondn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isplin, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isetup, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ispswt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isetsw, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.itimdp, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 10; ++i__) {
	do_fio(&c__1, (char *)&aersl1_1.iparm[i__ - 1], (ftnlen)sizeof(
		integer));
    }
    e_rsfe();

/*             FORCED DEFAULT VALUES */

    aersl5_1.issflg = 0;
    if (aersl5_1.ndiscr > 0) {
	goto L2;
    }
    nmulgv = 0;
    aersl5_1.ndsorc = 0;
    aersl5_1.ndremv = 0;
    aersl5_1.isstat = 0;

L2:
    if (aersl6_2.relerr <= 0.) {
	aersl6_2.relerr = .001;
    }
    if (aersl5_1.nvar < 10) {
	aersl5_1.nvar = 10;
    }
    if (aersl5_1.nmxvar < aersl5_1.nvar) {
	aersl5_1.nmxvar = aersl5_1.nvar;
    }
    if (aersl5_1.nprnt <= 0) {
	aersl5_1.nprnt = 1;
    }

    io___30.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___30);
    for (i__ = 1; i__ <= 10; ++i__) {
	do_fio(&c__1, (char *)&headng[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    do_fio(&c__1, (char *)&aersl6_2.xa, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.xb, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.relerr, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl6_2.abserr, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&aersl5_1.nvar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nmxvar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndiscr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nmulgv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nprnt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.iswtch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isstat, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.igrid, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.nptbet, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ncsorc, (ftnlen)sizeof(integer));
    e_wsfe();
    io___31.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___31);
    do_fio(&c__1, (char *)&aersl5_1.ncremv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndsorc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ndremv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.icondn, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isplin, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isetup, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.ispswt, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.isetsw, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&aersl5_1.itimdp, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 10; ++i__) {
	do_fio(&c__1, (char *)&aersl1_1.iparm[i__ - 1], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();

    io___32.ciunit = aersl5_1.iread;
    s_rsfe(&io___32);
    i__1 = aersl5_1.nprnt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&aersl6_2.tout[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_rsfe();
    io___33.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___33);
    i__1 = aersl5_1.nprnt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aersl6_2.tout[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wsfe();
    aersl5_1.nvarm1 = aersl5_1.nvar - 1;
    aersl5_1.nvarp1 = aersl5_1.nvar + 1;
    aersl5_1.nvarp2 = aersl5_1.nvar + 2;
    aersl5_1.nvart = aersl5_1.nvar + aersl5_1.ndiscr;
    aersl5_1.nmulp1 = nmulgv + 1;
    nmulp2 = nmulgv + 2;
    igdpm1 = abs(aersl5_1.igrid) + aersl5_1.nmulp1;
    igdpm2 = igdpm1 + 1;

/*        IF THE DISCRETE CONTINUOUS GDE IS USED READ IN THE MULTIPLET */
/*        GRID POINTS. FOR THE DISCRETE CONTINUOUS GDE, XA CORRESPONDS */
/*        TO THE MONOMER SIZE ON INPUT. AFTER COMPUTING THE CLUSTER */
/*        SIZES, XA IS THE SMALLEST SIZE IN THE CONTINUOUS REGIME */

    if (aersl5_1.ndiscr == 0) {
	goto L30;
    }
    aersl5_1.ndcls[0] = aersl5_1.ndiscr + 1;
    io___37.ciunit = aersl5_1.iread;
    s_rsfe(&io___37);
    i__1 = aersl5_1.nmulp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&aersl5_1.ndcls[i__ - 1], (ftnlen)sizeof(
		integer));
    }
    e_rsfe();
    aersl5_1.ndcls[aersl5_1.nmulp1 - 1] = aersl5_1.ndiscr << 1;
    io___38.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___38);
    i__1 = aersl5_1.nmulp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&aersl5_1.ndcls[i__ - 1], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();
    io___39.ciunit = aersl5_1.iread;
    s_rsfe(&io___39);
    i__1 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&aersl2_1.a[i__ - 1], (ftnlen)sizeof(doublereal)
		);
    }
    e_rsfe();
    i__1 = aersl5_1.ndiscr;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	aersl2_1.x[aersl5_1.nvar + i__ - 1] = aersl6_2.xa * (doublereal) i__;
    }
    i__1 = aersl5_1.nmulp1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L29: */
	aersl2_1.x[i__ - 1] = aersl6_2.xa * (doublereal) aersl5_1.ndcls[i__ - 
		1];
    }
    aersl6_2.xa = aersl2_1.x[0];

L30:
    aersl6_2.xboxa = aersl6_2.xb / aersl6_2.xa;
    aersl6_2.algxba = log(aersl6_2.xboxa);

/*                       ********************* */
/*        THE REMAINDER OF THE ROUTINE DETERMINES GRID POINT LOCATIONS, */
/*        THE VALUE OF THE DISTRIBUTION AT THE GRID POINTS, AND THE */
/*        QUADRATURE POINTS AND WEIGHT FACTORS OVER THE ENTIRE CONTINUOUS */
/*        DOMAIN. */
/*                       ********************* */

    aersl2_1.x[0] = aersl6_2.xa;
    aersl2_1.x[aersl5_1.nvar - 1] = aersl6_2.xb;
    aersl2_1.w[0] = 0.;
    aersl2_1.w[aersl5_1.nvar - 1] = 1.;

/*        IF USER SUPPLIES GRID POINTS, READ THEM IN AND COMPUTE THE */
/*        REMAINING GRID POINTS AS LOGARITHMICALLY SPACED IN X.  IF NO */
/*        GRID POINTS ARE SUPPLIED BY THE USER, ALL GRID POINTS WILL */
/*        BE LOGARITHMICALLY SPACED IN X.  GRID POINTS ARE READ IN FROM */
/*        SMALLEST TO LARGEST UNLESS IGRID IS NEGATIVE, THEN THE */
/*        ABSOLUTE VALUE OF IGRID POINTS ARE READ IN FROM LARGEST TO */
/*        SMALLEST. */

    if (nmulgv == 0) {
	goto L15;
    }
    i__1 = aersl5_1.nmulp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L21: */
	aersl2_1.w[i__ - 1] = log(aersl2_1.x[i__ - 1] / aersl6_2.xa) / 
		aersl6_2.algxba;
    }
L15:
    if (aersl5_1.igrid != 0) {
	goto L16;
    }

/*        HERE IF NO GRID POINTS ARE SUPPLIED BY USER */

    i__1 = aersl5_1.nvarm1;
    for (i__ = nmulp2; i__ <= i__1; ++i__) {
/* L14: */
	aersl2_1.w[i__ - 1] = (doublereal) (i__ - aersl5_1.nmulp1) / (
		doublereal) (aersl5_1.nvar - aersl5_1.nmulp1) * (1. - 
		aersl2_1.w[aersl5_1.nmulp1 - 1]) + aersl2_1.w[aersl5_1.nmulp1 
		- 1];
    }
    goto L7;

/*        HERE IF USER SUPPLIES GRID POINT LOCATIONS */

L16:
    if (aersl5_1.igrid > 0) {
	io___40.ciunit = aersl5_1.iread;
	s_rsfe(&io___40);
	i__1 = igdpm1;
	for (i__ = nmulp2; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&aersl2_1.x[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsfe();
    }
    if (aersl5_1.igrid < 0) {
	io___41.ciunit = aersl5_1.iread;
	s_rsfe(&io___41);
	i__1 = igdpm1;
	for (i__ = nmulp2; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&aersl2_1.x[aersl5_1.nvarm1 + nmulp2 - i__ 
		    - 1], (ftnlen)sizeof(doublereal));
	}
	e_rsfe();
    }
    i__1 = igdpm1;
    for (ii = nmulp2; ii <= i__1; ++ii) {
	i__ = ii;
	if (aersl5_1.igrid < 0) {
	    i__ = aersl5_1.nvarm1 + nmulp2 - ii;
	}
/* L9: */
	aersl2_1.w[i__ - 1] = log(aersl2_1.x[i__ - 1] / aersl6_2.xa) / 
		aersl6_2.algxba;
    }

/*        IF ALL THE GRID POINTS HAVE NOT BEEN SPECIFIED BY THE USER, */
/*        CALCULATE THE REMAINING GRID POINTS */

    if (igdpm1 >= aersl5_1.nvarm1) {
	goto L19;
    }
    if ((doublereal) aersl5_1.igrid < 0.) {
	goto L8;
    }

/*        HERE IF GRID POINTS ARE READ IN, IN ASCENDING ORDER */

    i__1 = aersl5_1.nvarm1;
    for (i__ = igdpm2; i__ <= i__1; ++i__) {
/* L17: */
	aersl2_1.w[i__ - 1] = (doublereal) (i__ - igdpm1) / (doublereal) (
		aersl5_1.nvar - igdpm1) * (1. - aersl2_1.w[igdpm1 - 1]) + 
		aersl2_1.w[igdpm1 - 1];
    }
    goto L7;

/*        HERE IF GRID POINTS ARE READ IN, IN DESCENDING ORDER */

L8:
    nvgdm1 = aersl5_1.nvar + aersl5_1.igrid - 1;
    i__1 = nvgdm1;
    for (i__ = nmulp2; i__ <= i__1; ++i__) {
/* L10: */
	aersl2_1.w[i__ - 1] = (doublereal) (i__ - aersl5_1.nmulp1) / (
		doublereal) (aersl5_1.nvar + aersl5_1.igrid - aersl5_1.nmulp1)
		 * (aersl2_1.w[aersl5_1.nvar + aersl5_1.igrid - 1] - 
		aersl2_1.w[aersl5_1.nmulp1 - 1]) + aersl2_1.w[aersl5_1.nmulp1 
		- 1];
    }

/*        NOW THAT ALL THE W VALUES OF THE GRID POINTS HAVE BEEN */
/*        DETERMINED CALCULATE THE PARTICLE SIZES (I.E. X(I)) AT THE */
/*        GRID POINTS WHICH WERE NOT SPECIFIED BY THE USER */

L7:
    i__1 = aersl5_1.nvarm1;
    for (ii = igdpm2; ii <= i__1; ++ii) {
	i__ = ii;
	if (aersl5_1.igrid < 0) {
	    i__ = aersl5_1.nmulp1 + ii - igdpm1;
	}
/* L34: */
	aersl2_1.x[i__ - 1] = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &
		aersl2_1.w[i__ - 1]);
    }

L19:
    io___44.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___44);
    i__1 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aersl2_1.x[i__ - 1], (ftnlen)sizeof(doublereal)
		);
    }
    e_wsfe();

/*        DETERMINE INITIAL VALUES OF THE DISTRIBUTION TO INITIALIZE A(I) */
/*        IF THE INITIAL TIME IS GREATER THAN ZERO, THIS MEANS THAT THE */
/*        USER IS RESTARTING A PROBLEM AND THE INITIAL VALUES OF A(I) */
/*        ARE TO BE READ IN */

/*     IF(TOUT(1).GT.0.D0)GO TO 13 */
    xintl_(&aersl5_1.nvar, aersl2_1.x, aersl2_1.a);
    i__1 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L6: */
	aersl2_1.a[i__ - 1] = aersl2_1.a[i__ - 1] * aersl6_2.algxba * 
		aersl2_1.x[i__ - 1];
    }
    goto L11;
/*  13 READ(IREAD,5)(A(I),I=1,NVAR) */
/*     DO 22 I=1,NVAR */
/*  22 A(I)=A(I)*ALGXBA*X(I) */

/*        COMPUTE THE QUADRATURE POINTS AND WEIGHT FACTORS FOR INTEGRATIO */
/*        OVER THE ENTIRE CONTINUOUS DOMAIN TO OBTAIN THE FIRST MOMENT */

L11:
    nquad2 = aersl5_1.nquadp / 2;
    i__1 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	i__2 = nquad2 - i__ + 1, i__3 = i__ - nquad2;
	j = max(i__2,i__3);
	aersl3_1.factor[i__ - 1] = aersl8_1.wtquad[j - 1] * .5;
	d__1 = (doublereal) (i__ - nquad2) - .5;
	aersl3_1.z__[i__ - 1] = (d_sign(&c_b174, &d__1) * aersl8_1.wquad[j - 
		1] + 1.) * .5;
/* L32: */
	aersl3_1.zmass[i__ - 1] = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &
		aersl3_1.z__[i__ - 1]);
    }
    return 0;
} 

/* Subroutine */ int output_(integer *istop)
{
    /* Format strings */
    static char fmt_6[] = "(\0021\002,49x,\002TIME =\002,e10.4,\002        \
 \002//22x,\002PARTICLE\002,10x,\002PARTICLE\002,6x,\002DN/ALOG(DP)\002,10x\
,\002DS/ALOG(DP)\002,3x,\002DV/ALOG(DP)\002/23x,\002VOLUME\002,11x,\002DIAME\
TER\002,8x,\002DISTRIBUTION\002,6x,\002DISTRIBUTION\002,6x,\002DISTRIBUTIO\
N\002/19x,\002V (CC)\002,11x,\002D (CM)\002,13x,\002N(LOG(D))\002,10x,\002S(\
LOG(D))\002,8x,\002V(LOG(D))\002/13x,\002NO.\002,9x,\002  \002,16x,\002  \
\002,14x,\002    \002,13x,\002        \002,11x,\002     \002/)";
    static char fmt_2[] = "(10x,i5,5x,e10.4,4e18.4)";
    static char fmt_7[] = "(/50x,\002INTERPOLATED DISTRIBUTION\002/)";
    static char fmt_13[] = "(/20x,\002MEAN VOLUME\002,6x,\002MEAN DIAMETE\
R\002,5x,\002TOTAL NUMBER\002,4x,\002TOTAL SURFACE AREA\002,3x,\002TOTAL VOL\
UME\002/24x,\002  \002,16x,\002  \002,16x,\002    \002,12x,\002        \002,\
11x,\002     \002/e30.4,4e18.4)";
    static char fmt_14[] = "(//56x,\002VOLUME ADDED\002//51x,\002MECHANIS\
M\002,5x,\002(CC/CC)\002/(i56,e18.4))";
    static char fmt_16[] = "(//54x,\002STEADY STATE\002)";
    static char fmt_9[] = "(//53x,\002DISCRETE REGIME\002//22x,\002PARTICL\
E\002,10x,\002PARTICLE\002,11x,\002NUMBER\002,9x,\002SURFACE AREA\002,8x,\
\002VOLUME\002/23x,\002VOLUME\002,11x,\002DIAMETER\002,7x,\002CONCENTRATIO\
N\002,6x,\002CONCENTRATION\002,4x,\002CONCENTRATION\002/25x,\002V\002,17x\
,\002D\002,17x,\002N\002,18x,\002S\002,16x,\002V\002/13x,\002NO.\002,9x,\002\
  \002,16x,\002  \002,14x,\002    \002,13x,\002        \002,11x,\002     \
\002/)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal d__;
    static integer i__, j;
    static doublereal wbet, xbet, dmean, vmean, distd;
    extern /* Subroutine */ int total_(doublereal *);
    static doublereal dists, distv;
    extern /* Subroutine */ int distw_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal zdist[24];
    static integer nspnr;
    static doublereal totalq[12];
    static integer ndspnr;
    static doublereal totnum, totvol, totsur;

    /* Fortran I/O blocks */
    static cilist io___47 = { 0, 0, 0, fmt_6, 0 };
    static cilist io___53 = { 0, 0, 0, fmt_2, 0 };
    static cilist io___54 = { 0, 0, 0, fmt_7, 0 };
    static cilist io___59 = { 0, 0, 0, fmt_2, 0 };
    static cilist io___65 = { 0, 0, 0, fmt_13, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_14, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_16, 0 };
    static cilist io___70 = { 0, 0, 0, fmt_9, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_2, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_13, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_14, 0 };


/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE PRINTS THE DISTRIBUTION AND SEVERAL MOMENTS OF   * */
/* *    THE DISTRIBUTION AT EACH OUTPUT TIME                             * */
/* ----------------------------------------------------------------------- */
    if (*istop == -2) {
	goto L15;
    }

/*        OUTPUT DISTRIBUTION AT THE GRID POINTS */

    io___47.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___47);
    do_fio(&c__1, (char *)&aersl6_2.time, (ftnlen)sizeof(doublereal));
    e_wsfe();
    i__1 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = aersl2_1.x[i__ - 1] * 1.90986;
	d__ = pow_dd(&d__1, &c_b181);
	distd = aersl2_1.a[i__ - 1] * 6.90776 / aersl6_2.algxba;
	dists = d__ * 3.14159 * d__ * distd;
	distv = aersl2_1.x[i__ - 1] * distd;
/* L1: */
	io___53.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aersl2_1.x[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&distd, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&dists, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&distv, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*        COMPUTE INTERPOLATION BETWEEN GRID POINTS */

    if (aersl5_1.nptbet == 0) {
	goto L4;
    }
    io___54.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___54);
    e_wsfe();
    i__1 = aersl5_1.nvar;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = aersl5_1.nptbet;
	for (j = 1; j <= i__2; ++j) {
	    wbet = aersl2_1.w[i__ - 2] + (aersl2_1.w[i__ - 1] - aersl2_1.w[
		    i__ - 2]) * (doublereal) j / (doublereal) (
		    aersl5_1.nptbet + 1);
	    xbet = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &wbet);
	    d__1 = xbet * 1.90986;
	    d__ = pow_dd(&d__1, &c_b181);
	    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &c__1, &wbet, zdist);
	    distd = zdist[0] * 6.90776 / aersl6_2.algxba;
	    dists = d__ * 3.14159 * d__ * distd;
	    distv = xbet * distd;
/* L3: */
	    io___59.ciunit = aersl5_1.iprnt;
	    s_wsfe(&io___59);
	    i__3 = j + aersl5_1.nptbet * (i__ - 2);
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&xbet, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&distd, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dists, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&distv, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }

/*        COMPUTE INTEGRAL MOMENTS OF THE DISTRIBUTION IN THE CONTINUOUS */
/*        REGIME */

L4:
    totnum = 0.;
    dmean = 0.;
    vmean = 0.;
    totsur = 0.;
    totvol = 0.;
    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, aersl3_2.z__,
	     zdist);
    i__3 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__3; ++i__) {
	zdist[i__ - 1] *= aersl3_2.factor[i__ - 1];
	totnum += zdist[i__ - 1];
	d__1 = aersl3_2.zvolm[i__ - 1] * 1.90986;
	d__ = pow_dd(&d__1, &c_b181);
	dmean += d__ * zdist[i__ - 1];
	totsur += d__ * 3.14159 * d__ * zdist[i__ - 1];
/* L10: */
	totvol += aersl3_2.zvolm[i__ - 1] * zdist[i__ - 1];
    }

    if (totnum < 1e-35) {
	goto L11;
    }
    dmean /= totnum;
    vmean = totvol / totnum;
L11:
    io___65.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&vmean, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dmean, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totnum, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totsur, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totvol, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*       NOTE THAT TOTAL CAN NOT BE CALLED TWICE FOR THE SAME IPTIM */

    if (*istop != 0 || aersl5_1.icondn == 1) {
	goto L8;
    }
    nspnr = aersl5_1.ncsorc + aersl5_1.ncremv;
    if (nspnr + aersl5_1.ndsorc + aersl5_1.ndremv > 0) {
	total_(totalq);
    }
    if (nspnr != 0) {
	io___68.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___68);
	i__3 = nspnr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&totalq[i__ - 1], (ftnlen)sizeof(doublereal)
		    );
	}
	e_wsfe();
    }
L8:
    if (aersl5_1.ndiscr == 0) {
	goto L303;
    }

/*            OUTPUT DISCRETE REGIME */

L15:
    if (*istop == -2 || aersl5_1.issflg == 1) {
	io___69.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___69);
	e_wsfe();
    }
    io___70.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___70);
    e_wsfe();
    totnum = 0.;
    dmean = 0.;
    totsur = 0.;
    totvol = 0.;
    i__3 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__3; ++i__) {
	d__1 = aersl2_1.x[i__ - 1] * 1.90986;
	d__ = pow_dd(&d__1, &c_b181);
	dists = d__ * 3.14159 * d__ * aersl2_1.a[i__ - 1];
	distv = aersl2_1.x[i__ - 1] * aersl2_1.a[i__ - 1];
	totnum += aersl2_1.a[i__ - 1];
	dmean += d__ * aersl2_1.a[i__ - 1];
	totsur += dists;
	totvol += distv;
/* L5: */
	io___71.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___71);
	i__2 = i__ - aersl5_1.nvar;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aersl2_1.x[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&d__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&aersl2_1.a[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, (char *)&dists, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&distv, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

    if (totnum < 1e-35) {
	goto L12;
    }
    vmean = totvol / totnum;
    dmean /= totnum;
L12:
    io___72.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___72);
    do_fio(&c__1, (char *)&vmean, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&dmean, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totnum, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totsur, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&totvol, (ftnlen)sizeof(doublereal));
    e_wsfe();
    ndspnr = aersl5_1.ndsorc + aersl5_1.ndremv;
    if (ndspnr != 0 && *istop == 0) {
	io___74.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___74);
	i__2 = ndspnr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&totalq[i__ + nspnr - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
    }
L303:
    return 0;
} /* output_ */

/* Subroutine */ int chktbl_(integer *istop)
{
    /* Format strings */
    static char fmt_6[] = "(\0021\002//40x,\002*** CURVE FIT NEGATIVE AT FOL\
LOWING POINTS\002//(55x,i4,e13.4))";
    static char fmt_23[] = "(\0020UNABLE TO ADD GRID POINTS SINCE SPECIFIED \
MAXIMUM HAS BEEN REACHED\002)";
    static char fmt_22[] = "(\0020*** NEGATIVE DISTRIBUTION IN USER SUPPLIED\
 MULTIPLET REGION, USER MUST RESTART\002)";
    static char fmt_11[] = "(15x,\002TIME STEP WILL BE RETRIED WITH THE ADDI\
TION OF THE FIRST\002,i3,\002 OF THE ABOVE POINTS IN PARTICLE SIZE\002)";
    static char fmt_4[] = "(\0021\002/55x,\002TIME =\002,e12.4/40x,\002*** D\
ISTRIBUTION NEGATIVE AT\rFOLLOWING GRID POINTS\002//(50x,i4,e13.4))";
    static char fmt_19[] = "(\0021\002/55x,\002TIME =\002,e12.4/40x,\002*** \
CLUSTER CONCENTRATION NEGATIVE AT FOLLOWING SIZES\002//(50x,i4,e13.4))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    double log(doublereal);

    /* Local variables */
    static integer i__, j, k, jj;
    static doublereal wbet[9], asave[60];
    extern /* Subroutine */ int distw_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal zdist[9];
    extern /* Subroutine */ int xintl_(integer *, doublereal *, doublereal *),
	     setup_(integer *, integer *, integer *);
    static integer nreplc, istart, ntroub;
    static doublereal wtroub[40], xtroub[60];
    static integer nvptrl;
    extern /* Subroutine */ int output_(integer *);

    /* Fortran I/O blocks */
    static cilist io___82 = { 0, 0, 0, fmt_6, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_23, 0 };
    static cilist io___84 = { 0, 0, 0, fmt_22, 0 };
    static cilist io___91 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_4, 0 };
    static cilist io___93 = { 0, 0, 0, fmt_19, 0 };


/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE CHECKS THE INTERPOLATION FOR NEGATIVE VALUES     * */
/* *    OF THE DISTRIBUTION AND WILL REGRESS TO THE PREVIOUS OUTPUT      * */
/* *    TIME AND ADD GRID POINTS AT THOSE LOCATIONS WHERE THE            * */
/* *    DISTRIBUTION IS NEGATIVE.  ISTOP WILL BE SET TO -1 TO INDICATE   * */
/* *    TO THE MAIN PROGRAM TO RESTART AND TRY AGAIN.  IF THE MAXIMUM    * */
/* *    NUMBER OF GRID POINTS HAS ALREADY BEEN REACHED, ISTOP WILL BE    * */
/* *    SET TO 1 AND THE CASE WILL BE SKIPPED IN THE MAIN PROGRAM.  IF   * */
/* *    ALL THE INTERPOLATED POINTS AT WHICH THE DISTRIBUTION IS         * */
/* *    NEGATIVE CAN NOT BE ADDED, AS MANY GRID POINTS AS THERE IS ROOM  * */
/* *    FOR WILL BE ADDED AND THE TIME STEP WILL BE RETIRED.  IF THE     * */
/* *    DISTRIBUTION IS NEGATIVE AT ONE OF THE GRID POINTS, THIS         * */
/* *    USUALLY INDICATES THAT THE TIME STEP WAS TOO LARGE AND SOME OF   * */
/* *    THE INTERPOLATED POINTS WERE NEGATIVE DURING THE TIME STEP WHICH * */
/* *    RESULTED IN A NEGATIVE DISTRIBUTION AT THE GRID POINTS.  IF THIS * */
/* *    HAPPENS AND NONE OF THE INTERPOLATED DISTRIBUTION POINTS ARE     * */
/* *    NEGATIVE OR THEY WERE NOT EVALUATED, ISTOP WILL BE SET TO 1 AND  * */
/* *    THE CASE WILL BE SKIPPED IN THE MAIN PROGRAM.                    * */
/* ----------------------------------------------------------------------- */
    *istop = 0;
    ntroub = 0;
    if (aersl5_1.nptbet == 0) {
	goto L1;
    }

/*        CHECK INTERPOLATED DISTRIBUTION FOR NEGATIVE VALUES */

    i__1 = aersl5_1.nvar;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = aersl5_1.nptbet;
	for (j = 1; j <= i__2; ++j) {
/* L3: */
	    wbet[j - 1] = aersl2_1.w[i__ - 2] + (aersl2_1.w[i__ - 1] - 
		    aersl2_1.w[i__ - 2]) * (doublereal) j / (doublereal) (
		    aersl5_1.nptbet + 1);
	}
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nptbet, wbet, 
		zdist);
	i__2 = aersl5_1.nptbet;
	for (j = 1; j <= i__2; ++j) {
	    if (zdist[j - 1] > 0.) {
		goto L27;
	    }
	    ++ntroub;
	    wtroub[ntroub - 1] = wbet[j - 1];
	    xtroub[ntroub - 1] = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &wbet[
		    j - 1]);
L27:
	    ;
	}
    }
    if (ntroub == 0) {
	goto L1;
    }
    *istop = 1;
    output_(istop);

/*        IF NO GRID POINTS CAN BE ADDED, SET ISTOP TO 1 AND RETURN */

    io___82.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___82);
    i__2 = ntroub;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xtroub[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    if (aersl5_1.nvar < aersl5_1.nmxvar) {
	goto L5;
    }
    io___83.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___83);
    e_wsfe();
    return 0;
L5:
    if (xtroub[0] > aersl2_1.x[aersl5_1.nmulp1 - 1]) {
	goto L16;
    }
    io___84.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___84);
    e_wsfe();
    return 0;
L16:
    nreplc = ntroub;

/*        ADD AS MANY GRID POINTS AS POSSIBLE.  IN CASE THERE IS TROUBLE */
/*        WITH THE FIRST DISTRIBUTION, PLACE THE A ARRAY IN THE ASAVE */
/*        ARRAY FOR THE FIRST DISTRIBUTION. */

    if (aersl5_1.nvar + ntroub > aersl5_1.nmxvar) {
	nreplc = aersl5_1.nmxvar - aersl5_1.nvar;
    }
    if (aersl5_1.iptim > 1) {
	goto L30;
    }
    i__2 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L24: */
	asave[i__ - 1] = aersl2_1.a[i__ - 1];
    }
L30:
    i__2 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L2: */
	aersl2_1.a[i__ - 1] = asave[i__ - 1];
    }

/*        ADD GRID POINTS ONE AT A TIME BY SHIFTING HIGHER GRID POINTS */
/*        UP. SINCE THE A ARRAY IS NEEDED FOR INTERPOLATION, PLACE THE */
/*        THE CORRESPONDING VALUES OF A FOR THE NEW GRID POINTS IN THE */
/*        ASAVE ARRAY.  SIMILARLY SHIFT THE X ARRAY BUT NOT THE W */
/*        ARRAY, SINCE THE W ARRAY IS NEEDED FOR INTERPOLATION */

    if (aersl5_1.isetup == 1 && aersl5_1.iptim > 1 || aersl5_1.isetup == 1 && 
	    aersl5_1.iptim == 1 && aersl6_2.tout[0] > 0.) {
	setup_(&aersl5_1.isplin, &c__0, &aersl5_1.nvar);
    }
    nvptrl = aersl5_1.nvar + nreplc;
    i__ = 1;
    i__2 = nreplc;
    for (k = 1; k <= i__2; ++k) {
	istart = i__ + 1;
	i__1 = nvptrl;
	for (i__ = istart; i__ <= i__1; ++i__) {
	    if (xtroub[k - 1] < aersl2_1.x[i__ - 1]) {
		goto L14;
	    }
/* L8: */
	}

L14:
	i__1 = aersl5_1.nvart;
	for (jj = i__; jj <= i__1; ++jj) {
	    j = aersl5_1.nvart + i__ - jj;
	    asave[j] = asave[j - 1];
/* L9: */
	    aersl2_1.x[j] = aersl2_1.x[j - 1];
	}
	aersl2_1.x[i__ - 1] = xtroub[k - 1];
	if (aersl5_1.iptim > 1 || aersl5_1.iptim == 1 && aersl6_2.tout[0] > 
		0.) {
	    goto L25;
	}
	xintl_(&c__1, &xtroub[k - 1], &asave[i__ - 1]);
	asave[i__ - 1] = asave[i__ - 1] * xtroub[k - 1] * aersl6_2.algxba;
	goto L26;
L25:
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &c__1, &wtroub[k - 1], 
		zdist);
	asave[i__ - 1] = zdist[0];
L26:
	++aersl5_1.nvart;
/* L13: */
    }

/*        NOW THAT ALL THE GRID POINTS HAVE BEEN ADDED, COMPUTE THE NEW */
/*        W ARRAY AND REPLACE THE A ARRAY.  ALSO INDICATE TO THE USER */
/*        WHAT CHANGES HAVE BEEN MADE AND RESET NVAR */

    aersl5_1.nvar = nvptrl;
    aersl5_1.nvarm1 = aersl5_1.nvar - 1;
    aersl5_1.nvarp1 = aersl5_1.nvar + 1;
    aersl5_1.nvarp2 = aersl5_1.nvar + 2;
    i__2 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	aersl2_1.w[i__ - 1] = log(aersl2_1.x[i__ - 1] / aersl6_2.xa) / 
		aersl6_2.algxba;
/* L10: */
	aersl2_1.a[i__ - 1] = asave[i__ - 1];
    }
    if (aersl5_1.ndiscr == 0) {
	goto L20;
    }
    i__2 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__2; ++i__) {
/* L21: */
	aersl2_1.a[i__ - 1] = asave[i__ - 1];
    }
L20:
    io___91.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___91);
    do_fio(&c__1, (char *)&nreplc, (ftnlen)sizeof(integer));
    e_wsfe();
    *istop = -1;
    if (aersl5_1.isetup == 1) {
	setup_(&aersl5_1.isplin, &c__1, &aersl5_1.nvar);
    }
    return 0;

/*        CHECK DISTRIBUTION AT GRID POINTS */

L1:
    i__2 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (aersl2_1.a[i__ - 1] > 0.) {
	    goto L12;
	}
	++ntroub;
	xtroub[ntroub - 1] = aersl2_1.x[i__ - 1];
L12:
	;
    }
    if (ntroub == 0) {
	goto L7;
    }
    *istop = 1;
    io___92.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___92);
    do_fio(&c__1, (char *)&aersl6_2.time, (ftnlen)sizeof(doublereal));
    i__2 = ntroub;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xtroub[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    output_(istop);

L7:
    if (aersl5_1.ndiscr == 0) {
	goto L17;
    }

/*        CHECK CLUSTER CONCENTRATIONS NDCLS ARRAY IS DESTROYED IF */
/*        ERRORS ARE FOUND */

    ntroub = 0;
    i__2 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__2; ++i__) {
	if (aersl2_1.a[i__ - 1] >= 0.) {
	    goto L18;
	}
	++ntroub;
	aersl5_1.ndcls[ntroub - 1] = i__ - aersl5_1.nvar;
	xtroub[ntroub - 1] = aersl2_1.a[i__ - 1];
L18:
	;
    }
    if (ntroub == 0) {
	goto L17;
    }
    *istop = 1;
    io___93.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___93);
    do_fio(&c__1, (char *)&aersl6_2.time, (ftnlen)sizeof(doublereal));
    i__2 = ntroub;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&aersl5_1.ndcls[i__ - 1], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, (char *)&xtroub[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

/*        NOW THAT ALL THE GRID POINTS AND INTERPOLATED POINTS HAVE BEEN */
/*        CHECKED, SAVE THE A ARRAY IN THE ASAVE ARRAY IN CASE THE */
/*        DISTRIBUTION AT THE NEXT OUTPUT TIME HAS ANY TROUBLE. */

L17:
    if (aersl5_1.nptbet == 0 || *istop == 1) {
	return 0;
    }
    i__2 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L15: */
	asave[i__ - 1] = aersl2_1.a[i__ - 1];
    }
    return 0;
} /* chktbl_ */

/* Subroutine */ int total_(doublereal *totalq)
{
    /* Initialized data */

    static doublereal ztime[4] = { -.86113631159405,-.33998104358486,
	    .33998104358486,.86113631159405 };
    static doublereal tfact[4] = { .34785484513745,.65214515486255,
	    .65214515486255,.34785484513745 };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum, sum2[6], avrgr[6], avrgx[3];
    extern /* Subroutine */ int distw_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal zdist[24];
    static integer ntotq, ncrem1;
    static doublereal balanc, remove[72]	/* was [3][24] */, source[72]	
	    /* was [3][24] */, tfinsh;
    extern /* Subroutine */ int csourc_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    static integer ntsorc;
    static doublereal dumtim;
    extern /* Subroutine */ int dsourc_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer ntremv;
    static doublereal totmas, startm, tstart;

/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE CALCULATES THE TOTAL MASS OF AEROSOL ADDED AND   * */
/* *    REMOVED BY EACH MECHANISM FROM TIME=TOUT(IPTIM-1) TO TIME=       * */
/* *    TOUT(IPTIM).  IF THIS ROUTINE IS CALLED AT THE INITIAL TIME (I.E.* */
/* *    IPTIM=1), THE TOTAL MASSES ADDED/REMOVED ARE INITIALIZED TO ZERO * */
/* *    AND THE REMOVAL RATES ARE CALCULATED AND STORED.                 * */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --totalq;

    /* Function Body */

    ntotq = aersl5_1.ncsorc + aersl5_1.ncremv + aersl5_1.ndsorc + 
	    aersl5_1.ndremv;
    if (ntotq == 0) {
	return 0;
    }

/*        INITIALIZE TOTAL MASS ADDED AND REMOVED */

    if (aersl5_1.iptim != 1) {
	goto L2;
    }
    i__1 = ntotq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L4: */
	totalq[i__] = 0.;
    }



/*        COMPUTE MASS ADDED TO THE CONTINUOUS REGIME FROM THE LAST */
/*        OUTPUT TIME TO THE CURRENT OUTPUT TIME BY USING A NQUADP */
/*        QUADRATURE FORMULA IN X AND A 4 POINT QUADRATURE FORMULA */
/*        IN TIME */

L2:
    ntsorc = aersl5_1.ncsorc + aersl5_1.ndsorc;
    if (aersl5_1.ncsorc + aersl5_1.ncremv == 0) {
	goto L10;
    }
    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, aersl3_1.z__,
	     zdist);
    i__1 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	zdist[i__ - 1] /= aersl3_1.zmass[i__ - 1] * aersl6_2.algxba;
    }
L10:
    if (aersl5_1.iptim == 1 || ntsorc == 0) {
	goto L17;
    }
    tstart = aersl6_2.tout[aersl5_1.iptim - 2];
    tfinsh = aersl6_2.tout[aersl5_1.iptim - 1];
    if (aersl5_1.ncsorc == 0) {
	goto L11;
    }
    i__1 = aersl5_1.ncsorc;
    for (j = 1; j <= i__1; ++j) {
/* L31: */
	sum2[j - 1] = 0.;
    }

    for (k = 1; k <= 4; ++k) {
	i__1 = aersl5_1.ncsorc;
	for (j = 1; j <= i__1; ++j) {
/* L33: */
	    avrgx[j - 1] = 0.;
	}
	d__1 = (ztime[k - 1] * (tfinsh - tstart) + tfinsh + tstart) * .5;
	csourc_(&aersl5_1.nquadp, &d__1, aersl3_1.zmass, zdist, source, 
		remove);
	i__1 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = aersl5_1.ncsorc;
	    for (j = 1; j <= i__2; ++j) {
/* L34: */
/* Computing 2nd power */
		d__1 = aersl3_1.zmass[i__ - 1];
		avrgx[j - 1] += source[j + i__ * 3 - 4] * aersl3_1.factor[i__ 
			- 1] * (d__1 * d__1);
	    }
	}
	i__2 = aersl5_1.ncsorc;
	for (j = 1; j <= i__2; ++j) {
/* L32: */
	    sum2[j - 1] += tfact[k - 1] * avrgx[j - 1];
	}
    }

    i__2 = aersl5_1.ncsorc;
    for (j = 1; j <= i__2; ++j) {
/* L35: */
	totalq[j] += (tfinsh - tstart) * .5 * sum2[j - 1] * aersl6_2.algxba;
    }
L11:
    if (aersl5_1.ndsorc == 0) {
	goto L17;
    }

/*        COMPUTE MASS ADDED IN THE DISCRETE REGIME FROM THE PAST */
/*        OUTPUT TIME TO THE CURRENT OUTPUT TIME */

    i__2 = aersl5_1.ndsorc;
    for (j = 1; j <= i__2; ++j) {
/* L24: */
	sum2[j - 1] = 0.;
    }

    for (k = 1; k <= 4; ++k) {
	dumtim = (ztime[k - 1] * (tfinsh - tstart) + tfinsh + tstart) * .5;
	i__2 = aersl5_1.ndsorc;
	for (j = 1; j <= i__2; ++j) {
/* L38: */
	    avrgx[j - 1] = 0.;
	}
	dsourc_(&aersl5_1.ndiscr, &dumtim, &aersl2_1.a[aersl5_1.nvarp1 - 1], 
		source, remove);
	i__2 = aersl5_1.ndiscr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = aersl5_1.ndsorc;
	    for (j = 1; j <= i__1; ++j) {
/* L36: */
		avrgx[j - 1] += source[j + i__ * 3 - 4] * aersl2_1.x[i__ + 
			aersl5_1.nvar - 1];
	    }
	}
	i__1 = aersl5_1.ndsorc;
	for (j = 1; j <= i__1; ++j) {
/* L37: */
	    sum2[j - 1] += tfact[k - 1] * avrgx[j - 1];
	}
    }

    i__1 = aersl5_1.ndsorc;
    for (j = 1; j <= i__1; ++j) {
/* L26: */
	totalq[j + aersl5_1.ncsorc + aersl5_1.ncremv] += sum2[j - 1] * .5 * (
		tfinsh - tstart);
    }

/*        COMPUTE THE AVERAGE MASS REMOVAL RATE OF EACH MECHANISM FOR */
/*        TIME STEP BY SUMMING THE REMOVAL RATES OF EACH MECHANISM AT */
/*        THE PREVIOUS OUTPUT TIME AND THE CURRENT OUTPUT TIME, AND */
/*        DIVIDING BY 2.  THIS AVERAGE REMOVAL RATE IS STORED IN */
/*        AVRGR(I) FOR THE I-TH MECHANISM.  TO DETERMINE THE MASS */
/*        REMOVED FOR THE TIME STEP, COMPUTE THE TOTAL MASS REMOVED FOR */
/*        TIME STEP AND STORE THAT VALUE IN BALANC.  THEN THE MASS */
/*        REMOVED BY THE I-TH MECHANISM IS TOTALQ(I+NCSORC)=AVRGR(I)* */
/*        BALANC DIVIDED BY THE SUM OF ALL AVRGR(J),J=1,NCREMV.  FOR THE */
/*        INITIAL CALL (I.E. IPTIM=1), ONLY STORE THE REMOVAL RATE FOR */
/*        EACH MECHANISM IN THE ARRAY AVRGR.  FOR SUBSEQUENT TIMES, */
/*        SUM2(I) IS THE REMOVAL RATE OF THE I-TH MECHANISM AT THE */
/*        CURRENT TIME.  AVRGR(I) IS USED TO TEMPORARILY STORE THE */
/*        AVERAGE REMOVAL RATE OVER THE TIME STEP, BUT BEFORE */
/*        RETURNING, AVRGR(I) STORES THE REMOVAL RATE AT THE CURRENT */
/*        OUTPUT TIME.  THE SAME METHOD IS USED IN THE DISCRETE REGIME. */

L17:
    ntremv = aersl5_1.ncremv + aersl5_1.ndremv;
    if (ntremv == 0) {
	return 0;
    }
    totmas = 0.;
    i__1 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L18: */
/* Computing 2nd power */
	d__1 = aersl3_1.zmass[i__ - 1];
	totmas += aersl3_1.factor[i__ - 1] * zdist[i__ - 1] * (d__1 * d__1);
    }
    totmas *= aersl6_2.algxba;

    if (aersl5_1.ndiscr == 0) {
	goto L19;
    }
    i__1 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__1; ++i__) {
/* L20: */
	totmas += aersl2_1.x[i__ - 1] * aersl2_1.a[i__ - 1];
    }

L19:
    if (aersl5_1.iptim == 1) {
	startm = totmas;
    }
    i__1 = ntremv;
    for (j = 1; j <= i__1; ++j) {
/* L21: */
	sum2[j - 1] = 0.;
    }
    if (aersl5_1.ncremv == 0) {
	goto L39;
    }

    csourc_(&aersl5_1.nquadp, &aersl6_2.time, aersl3_1.zmass, zdist, source, 
	    remove);
    i__1 = aersl5_1.nquadp;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = aersl3_1.zmass[j - 1];
	zdist[j - 1] = aersl3_1.factor[j - 1] * (d__1 * d__1);
	i__2 = aersl5_1.ncremv;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L22: */
	    sum2[i__ - 1] += remove[i__ + j * 3 - 4] * zdist[j - 1];
	}
    }

    if (aersl5_1.iptim > 1) {
	goto L5;
    }
    i__2 = aersl5_1.ncremv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L6: */
	avrgr[i__ - 1] = sum2[i__ - 1];
    }
    goto L16;
L5:
    i__2 = aersl5_1.ncremv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L9: */
	avrgr[i__ - 1] = (sum2[i__ - 1] + avrgr[i__ - 1]) * .5;
    }

/*        COMPUTE MASS REMOVED IN DISCRETE REGIME */

L16:
    if (aersl5_1.ndiscr == 0) {
	goto L12;
    }
L39:
    dsourc_(&aersl5_1.ndiscr, &aersl6_2.time, &aersl2_1.a[aersl5_1.nvarp1 - 1]
	    , source, remove);
    i__2 = aersl5_1.ndiscr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = aersl5_1.ndremv;
	for (j = 1; j <= i__1; ++j) {
/* L42: */
	    sum2[j + aersl5_1.ncremv - 1] += aersl2_1.x[i__ + aersl5_1.nvar - 
		    1] * remove[j + i__ * 3 - 4];
	}
    }

    if (aersl5_1.iptim > 1) {
	goto L7;
    }
    ncrem1 = aersl5_1.ncremv + 1;
    i__1 = ntremv;
    for (i__ = ncrem1; i__ <= i__1; ++i__) {
/* L8: */
	avrgr[i__ - 1] = sum2[i__ - 1];
    }
    return 0;
L7:
    i__1 = ntremv;
    for (i__ = ncrem1; i__ <= i__1; ++i__) {
/* L3: */
	avrgr[i__ - 1] = (sum2[i__ - 1] + avrgr[i__ - 1]) * .5;
    }

/*        FROM A MASS BALANCE DETERMINE THE AEROSOL MASS PER UNIT VOLUME */
/*        WHICH WAS REMOVED SINCE THE LAST STEP AND STORE THAT IN BALANC */

L12:
    if (aersl5_1.iptim == 1) {
	return 0;
    }
    balanc = totmas - startm;
    i__1 = ntotq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	balanc -= totalq[i__];
    }

/*        IF MASS WAS REMOVED (I.E. BALANC.LT.0.) DETERMINE THE AMOUNT */
/*        REMOVED BY EACH MECHANISM */

    if (balanc > 0.) {
	goto L13;
    }
    sum = 0.;
    i__1 = ntremv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L27: */
	sum += avrgr[i__ - 1];
    }
    i__1 = aersl5_1.ncremv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	totalq[i__ + aersl5_1.ncsorc] = avrgr[i__ - 1] * balanc / sum + 
		totalq[i__ + aersl5_1.ncsorc];
    }
    if (aersl5_1.ndremv == 0) {
	goto L13;
    }
    i__1 = aersl5_1.ndremv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L23: */
	totalq[i__ + ntsorc + aersl5_1.ncremv] = avrgr[i__ + aersl5_1.ncremv 
		- 1] * balanc / sum + totalq[i__ + ntsorc + aersl5_1.ncremv];
    }

/*        STORE REMOVAL RATE OF THE CURRENT OUTPUT TIME */

L13:
    i__1 = ntremv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L14: */
	avrgr[i__ - 1] = sum2[i__ - 1];
    }
    return 0;
} /* total_ */

/* Subroutine */ int coagcf_(integer *istop)
{
    /* Initialized data */

    static doublereal xe[2] = { .33998104358486,.86113631159405 };
    static doublereal we[2] = { .65214515486255,.34785484513745 };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), d_sign(doublereal *, doublereal *), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, mw;
    static doublereal dum, txa;
    extern /* Subroutine */ int beta_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), evap_(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer maxg;
    static doublereal bwmz1;
    static integer nbefr;
    static doublereal xbxaw, upper;
    static integer nquad2, jpnvar;
    static doublereal bottom;

/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE COMPUTES THE QUADRATURE POINTS AND COAGULATION   * */
/* *    COEFFICIENTS AT THE QUADRATURE POINTS FOR EVALUATING THE         * */
/* *    COAGULATION INTEGRALS.  THE EVAPORATION COEFFICIENTS ARE ALSO    * */
/* *    EVALUATED AND STORED                                             * */
/* ----------------------------------------------------------------------- */



/*        DETERMINE THE QUADRATURE POINTS AND FACTORS NEEDED TO COMPUTE */
/*        THE FIRST COAGUALTION INTEGRAL, I.E. THE INTEGRAL OF */

/*                        Y           Z1                        W(MW) */
/*         BETA(XA*(XB/XA) ,XA*(XB/XA)  )*M(Y,T)*M(Z1,T)*(XB/XA) */
/*         ---------------------------------------------------------- */
/*                          W(MW)                      Z1 */
/*                   (XB/XA)         -          (XB/XA) */


/*                                   Z1 */
/*              - BETA(X(MW),XA*XBOXA  )*M(Z1,T)*M(W(MW),T) */


/*        AS Z1 VARIES FROM ZERO TO  W(MW)-LN(2)/LN(XB/XA) */
/*        WHERE */
/*                      (       W(MW)                Z1) */
/*                 Y= LN((XB/XA)       -      (XB/XA)  ) */
/*                    ---------------------------------- */
/*                                   LN(XB/XA) */

/*        M(W,T)=DISTRIBUTION FUNCTION OF DIMENSIONLESS PARTICLE SIZE W */
/*          AND TIME T */
/*        A NQUADP POINT GAUSSIAN LEGENDRE QUADRATURE IS USED, TAKEN FROM */
/*        THE HANDBOOK OF MATHEMATICAL FUNCTIONS BY MILTON ABRAMOWITZ AND */
/*        IRENE STEGUN, DOVER PUBLICATIONS INC., N.Y. TABLE 25.4 PAGE */
/*        917. FOR A COMPLETE DESCRIPTION OF QUADRATURE SEE B.CARNAHAN, */
/*        H.A.LUTHER AND J.O.WILKES,'APPLIED NUMERICAL METHODS',JOHN */
/*        WILEY AND SONS, PAGES 100-105 (1969) */

    nquad2 = aersl5_1.nquadp / 2;
    txa = log(2.) / aersl6_2.algxba;
    i__1 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L14: */
	if (aersl2_1.w[aersl5_1.nvar - i__] > txa) {
	    aersl5_1.minw = aersl5_1.nvar - i__ + 1;
	}
    }

    i__1 = aersl5_1.nvar;
    for (mw = aersl5_1.minw; mw <= i__1; ++mw) {
	upper = (aersl2_1.w[mw - 1] - txa) * .5;
	xbxaw = aersl2_1.x[mw - 1] / aersl6_2.xa;
	i__2 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = nquad2 - i__ + 1, i__4 = i__ - nquad2;
	    j = max(i__3,i__4);
/* L3: */
	    d__1 = (doublereal) (i__ - nquad2) - .5;
	    aersl4_1.z1[i__ + mw * 24 - 25] = upper * (d_sign(&c_b174, &d__1) 
		    * aersl8_1.wquad[j - 1] + 1.);
	}
	i__2 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L12: */
	    aersl4_1.z2[i__ + mw * 24 - 25] = aersl6_2.xa * pow_dd(&
		    aersl6_2.xboxa, &aersl4_1.z1[i__ + mw * 24 - 25]);
	}
	beta_(&aersl5_1.nquadp, &aersl6_2.time, &aersl2_1.x[mw - 1], &
		aersl4_1.z2[mw * 24 - 24], &aersl4_1.fact1b[mw * 24 - 24]);
	i__2 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = nquad2 - i__ + 1, i__4 = i__ - nquad2;
	    j = max(i__3,i__4);
	    aersl4_1.fact1b[i__ + mw * 24 - 25] = aersl4_1.fact1b[i__ + mw * 
		    24 - 25] * aersl8_1.wtquad[j - 1] * upper;
	    bwmz1 = xbxaw - pow_dd(&aersl6_2.xboxa, &aersl4_1.z1[i__ + mw * 
		    24 - 25]);
	    aersl4_1.wmz1[aersl5_1.nquadp - i__ + 1 + mw * 24 - 25] = log(
		    bwmz1) / aersl6_2.algxba;
	    d__1 = aersl6_2.xa * bwmz1;
	    beta_(&c__1, &aersl6_2.time, &aersl4_1.z2[i__ + mw * 24 - 25], &
		    d__1, &aersl4_1.fact1a[i__ + mw * 24 - 25]);
/* L5: */
	    aersl4_1.fact1a[i__ + mw * 24 - 25] = aersl4_1.fact1a[i__ + mw * 
		    24 - 25] * aersl8_1.wtquad[j - 1] * upper * xbxaw / bwmz1;
	}
    }

/*        DETERMINE THE QUADRATURE POINTS AND FACTORS NEEDED TO COMPUTE */
/*        SECOND COAGUATION INTEGRAL, I.E. THE INTEGRAL OF */

/*         BETA(X(MW),XA*XBOXA**Z2,TIME)*M(Z2,TIME)*M(W(MW),TIME) */

/*        AS Z2 VARIES FROM W(MW)-LN(2)/LN(XB/XA) TO ONE.  IF THIS LOWER */
/*        LIMIT IS NEGATIVE, THEN SUBSTITUTE ZERO FOR THE LOWER LIMIT. */

    i__2 = aersl5_1.nvar;
    for (mw = 1; mw <= i__2; ++mw) {
/* Computing MAX */
	d__1 = 0., d__2 = aersl2_1.w[mw - 1] - txa;
	bottom = max(d__1,d__2);
	i__1 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    i__3 = nquad2 - i__ + 1, i__4 = i__ - nquad2;
	    j = max(i__3,i__4);
	    d__1 = (doublereal) (i__ - nquad2) - .5;
	    aersl4_1.z2[i__ + mw * 24 - 25] = (aersl8_1.wquad[j - 1] * d_sign(
		    &c_b174, &d__1) * (1. - bottom) + (float)1. + bottom) * 
		    .5;
/* L15: */
	    aersl4_1.facte1[i__ - 1] = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &
		    aersl4_1.z2[i__ + mw * 24 - 25]);
	}
	beta_(&aersl5_1.nquadp, &aersl6_2.time, &aersl2_1.x[mw - 1], 
		aersl4_1.facte1, &aersl4_1.fact2[mw * 24 - 24]);
	i__1 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    i__3 = nquad2 - i__ + 1, i__4 = i__ - nquad2;
	    j = max(i__3,i__4);
/* L8: */
	    aersl4_1.fact2[i__ + mw * 24 - 25] = aersl4_1.fact2[i__ + mw * 24 
		    - 25] * aersl8_1.wtquad[j - 1] * (1. - bottom) * .5;
	}
    }

/*                  **************************************** */
/*        THE REST OF THE ROUTINE IS USED ONLY IF THE DISCRETE - */
/*        CONTINUOUS G-D-E IS REQUESTED */
/*                  **************************************** */

    if (aersl5_1.ndiscr == 0) {
	return 0;
    }

/*        IN CASE GRID POINTS HAVE BEEN ADDED (I.E. ISTOP=-1), DO NOT */
/*        WASTE COMPUTER TIME RECOMPUTING CLUSTER-CLUSTER COAGULATION */
/*        AND EVAPORATION COEFFICIENTS */

    if (*istop == -1) {
	goto L2;
    }

/*        COMPUTE THE CLUSTER - CLUSTER COAGULATION COEFFICIENTS, NOTE */
/*        THAT THE COEFFICIENT IS A SYMMETRIC FUNCTION OF ITS ARGUMENTS */
/*        AND THUS A SIGNIFICANT AMOUNT OF STORAGE CAN BE SAVED */

    i__1 = aersl5_1.ndiscr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nbefr = (i__ - 1) * aersl5_1.ndiscr - (i__ - 1) * (i__ - 2) / 2;
/* L6: */
	i__2 = aersl5_1.ndiscr - i__ + 1;
	beta_(&i__2, &aersl6_2.time, &aersl2_1.x[i__ + aersl5_1.nvar - 1], &
		aersl2_1.x[i__ + aersl5_1.nvar - 1], &aersl4_1.betdsr[nbefr]);
    }

/*        COMPUTE THE FOUR QUADRATURE POINTS AND FACTORS NEEDED FOR */
/*        DETERMINING THE FORMATION RATE OF THE LARGEST CLUSTER DUE TO */
/*        EVAPORATION IN THE DISCRETE REGIME */

    upper = log((aersl6_2.xa + aersl2_1.x[aersl5_1.nvarp1 - 1]) / aersl6_2.xa)
	     / aersl6_2.algxba;
    for (i__ = 1; i__ <= 4; ++i__) {
/* Computing MAX */
	i__2 = 3 - i__, i__1 = i__ - 2;
	j = max(i__2,i__1);
	d__1 = (doublereal) i__ - 2.5;
	aersl4_1.ze[i__ - 1] = upper * .5 * (d_sign(&c_b174, &d__1) * xe[j - 
		1] + 1.);
	d__1 = aersl6_2.xa * pow_dd(&aersl6_2.xboxa, &aersl4_1.ze[i__ - 1]);
	evap_(&c__1, &aersl6_2.time, &d__1, &aersl4_1.facte3[i__ - 1]);
/* L11: */
	aersl4_1.facte3[i__ - 1] = aersl4_1.facte3[i__ - 1] * upper * .5 * we[
		j - 1];
    }

/*        COMPUTE THE QUADRATURE FACTORS FOR THE NET RATE OF MONOMER */
/*        EVAPORATION FROM THE CONTINUOUS REGIME */

    evap_(&aersl5_1.nquadp, &aersl6_2.time, aersl3_1.zmass, aersl4_1.facte1);
    i__2 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L7: */
	aersl4_1.facte1[i__ - 1] *= aersl3_1.factor[i__ - 1];
    }

/*       COMPUTE THE QUADRATURE FACTORS NEEDED FOR THE SCAVENGING OF */
/*       CLUSTERS BY PARTICLES IN THE CONTINUOUS REGIME */

L2:
    i__2 = aersl5_1.nvart;
    for (mw = aersl5_1.nvarp1; mw <= i__2; ++mw) {
	beta_(&aersl5_1.nquadp, &aersl6_2.time, &aersl2_1.x[mw - 1], 
		aersl3_1.zmass, &aersl4_1.fact2[mw * 24 - 24]);
	i__1 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L4: */
	    aersl4_1.fact2[i__ + mw * 24 - 25] *= aersl3_1.factor[i__ - 1];
	}
    }

/*        COMPUTE THE SIZE AND EVAPORATION COEFFICIENTS OF PARTICLES ONE */
/*        MONOMER LARGER THAN THE PARTICLES AT THE GRID POINTS */

    aersl5_1.maxevp = 0;
    i__1 = aersl5_1.nvarm1;
    for (mw = 1; mw <= i__1; ++mw) {
	aersl4_1.evpstr[mw - 1] = aersl2_1.x[mw - 1] + aersl2_1.x[
		aersl5_1.nvarp1 - 1];
	if (aersl4_1.evpstr[mw - 1] > aersl6_2.xb) {
	    goto L9;
	}
	aersl4_1.wpls[mw - 1] = log(aersl4_1.evpstr[mw - 1] / aersl6_2.xa) / 
		aersl6_2.algxba;
	aersl5_1.maxevp = mw;
/* L29: */
    }
L9:
    evap_(&aersl5_1.maxevp, &aersl6_2.time, aersl4_1.evpstr, aersl4_1.evppls);
    i__1 = aersl5_1.maxevp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L16: */
	aersl4_1.evppls[i__ - 1] = aersl2_1.x[i__ - 1] * aersl4_1.evppls[i__ 
		- 1] / aersl4_1.evpstr[i__ - 1];
    }

/*        STORE THE EVAPORATION COEFFICIENTS AT THE GRID POINTS AND FOR */
/*        FOR CLUSTERS. (NOTE THAT EVPSTR(NVAR+1) SHOULD NEVER BE USED, */
/*        HOWEVER FOR PROGRAM CLARITY EVPSTR(NVAR+I) IS THE EVAPORATION */
/*        COEFFICIENT OF THE I-TH CLUSTER). */

    evap_(&aersl5_1.nvar, &aersl6_2.time, aersl2_1.x, aersl4_1.evpstr);
    i__1 = aersl5_1.ndiscr - 1;
    evap_(&i__1, &aersl6_2.time, &aersl2_1.x[aersl5_1.nvarp2 - 1], &
	    aersl4_1.evpstr[aersl5_1.nvarp2 - 1]);

/*        COMPUTE THE W SIZES (I.E. WD) AND COAGULATION COEFFICIENTS (I.E */
/*        BESR) FOR THE NET RATE OF GROWTH DUE TO A DISCRETE CLUSTER */
/*        COAGULATING WITH A PARTICLE IN THE CONTINUOUS REGIME */

    aersl4_1.maxgdc[0] = 0;
    i__1 = aersl5_1.nvar;
    for (mw = 2; mw <= i__1; ++mw) {
	maxg = aersl5_1.ndiscr;
	if (mw <= aersl5_1.nmulp1) {
	    maxg = aersl5_1.ndcls[mw - 1] - aersl5_1.ndiscr - 1;
	}
	aersl4_1.maxgdc[mw - 1] = maxg;
	i__2 = maxg;
	for (j = 1; j <= i__2; ++j) {
	    jpnvar = j + aersl5_1.nvar;
	    dum = aersl2_1.x[mw - 1] - aersl2_1.x[jpnvar - 1];
	    beta_(&c__1, &aersl6_2.time, &dum, &aersl2_1.x[jpnvar - 1], &
		    aersl4_1.besr[j + mw * 20 - 21]);
	    aersl4_1.besr[j + mw * 20 - 21] = aersl4_1.besr[j + mw * 20 - 21] 
		    * aersl2_1.x[mw - 1] / dum;
/* L27: */
	    aersl4_1.wd[j + mw * 20 - 21] = log(dum / aersl6_2.xa) / 
		    aersl6_2.algxba;
	}
    }

/*       COMPUTE THE COAGULATION COEFFICIENTS BETWEEN CLUSTERS AND */
/*       PARTICLES IN THE CONTINUOUS REGIME */

    i__2 = aersl5_1.nvar;
    for (j = 1; j <= i__2; ++j) {
/* L30: */
	beta_(&aersl5_1.ndiscr, &aersl6_2.time, &aersl2_1.x[j - 1], &
		aersl2_1.x[aersl5_1.nvarp1 - 1], &aersl4_1.betadc[j * 20 - 20]
		);
    }
    return 0;
} /* coagcf_ */

/* ----------------------------------------------------------------------- */
/* *    THIS DATA BLOCK CONTAIN THE QUADRATURE POINTS AND WEIGHT         * */
/* *    FACTORS FOR THE GUASSIAN-LEGENDRE QUADRATURES.  NQUADP IS        * */
/* *    THE NUMBER OF QUADRATURE POINTS AND IT MUST BE AN EVEN NUMBER    * */
/* ----------------------------------------------------------------------- */

/*     DATA NQUADP/12/ */
/*     DATA WQUAD/.125233408511469  ,.367831498998180  , */
/*    $.587317954286617  ,.769902674194305  , */
/*    $.904117256370475  ,.981560634246719  / */
/*     DATA WTQUAD/.249147045813403  ,.233492536538355  , */
/*    $.203167426723066  ,.160078328543346  , */
/*    $.106939325995318  ,.047175336386512  / */

/*     DATA NQUADP/40/ */
/*     DATA WQUAD/.038772417506050821933  ,.116084070675255208483  , */
/*    $.192697580701371099716  ,.268152185007253681141  , */
/*    $.341994090825758473007  ,.413779204371605001525  , */
/*    $.483075801686178712909  ,.549467125095128202076  , */
/*    $.612553889667980237953  ,.671956684614179548379  , */
/*    $.727318255189927103281  ,.778305651426519387695  , */
/*    $.824612230833311663196  ,.865959503212259503821  , */
/*    $.902098806968874296728  ,.932812808278676533361  , */
/*    $.957916819213791655805  ,.977259949983774262663  , */
/*    $.990726238699457006453  ,.998237709710559200350  / */
/*     DATA WTQUAD/.077505947978424811264  ,.077039818164247965588  , */
/*    $.076110361900626242372  ,.074723169057968264200  , */
/*    $.072886582395804059061  ,.070611647391286779695  , */
/*    $.067912045815233903826  ,.064804013456601038075  , */
/*    $.061306242492928939167  ,.057439769099391551367  , */
/*    $.053227846983936824355  ,.048695807635072232061  , */
/*    $.043870908185673271992  ,.038782167974472017640  , */
/*    $.033460195282547847393  ,.027937006980023401098  , */
/*    $.022245849194166957262  ,.016421058381907888713  , */
/*    $.010498284531152813615  ,.004521277098533191258  / */


/* Subroutine */ int diffun_(doublereal *time, doublereal *adif, doublereal *
	dadt)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, ii, mw;
    static doublereal sum, sum2, term;
    static integer msize;
    extern /* Subroutine */ int distw_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal zdist[40];
    extern /* Subroutine */ int setup_(integer *, integer *, integer *);
    static integer msmls, maxsz, nvarp3;
    static doublereal zdist2[40];
    extern /* Subroutine */ int coagcf_(integer *);
    static doublereal remove[120]	/* was [3][40] */, source[120]	/* 
	    was [3][40] */;
    extern /* Subroutine */ int steady_(doublereal *), derivt_(integer *, 
	    integer *, doublereal *);
    extern doublereal boundr_(doublereal *);
    extern /* Subroutine */ int growth_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *), csourc_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dsourc_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static integer mwmnvr;

/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE COMPUTES THE DERIVATIVE OF ADIF(I) WITH RESPECT  * */
/* *    TO TIME AND STORES THE RESULT IN DADT(I).  FOR THE DISCRETE -    * */
/* *    CONTINUOUS G-D-E, DADT(I) IS THE DERIVATIVE OF THE CLUSTER       * */
/* *    CONCENTRATIONS WITH RESPECT TO TIME FOR I=NVARP1 TO NVART, WHERE * */
/* *    DADT(NVARP1) IS FOR THE MONOMER AND DADT(NVARP1+1) IS FOR THE    * */
/* *    DIMER, ETC.)                                                     * */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --dadt;
    --adif;

    /* Function Body */
    i__1 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L3: */
	aersl2_1.a[i__ - 1] = adif[i__];
    }
    aersl6_1.dumtim = *time;
    if (aersl5_1.isetup == 1) {
	setup_(&aersl5_1.isplin, &c__0, &aersl5_1.nvar);
    }
    if (aersl5_1.itimdp == 1) {
	coagcf_(&c__0);
    }
    if (aersl5_1.issflg == 1) {
	steady_(&aersl6_1.dumtim);
    }

    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
/* L13: */
	dadt[mw] = 0.;
    }

/*        COMPUTE THE FIRST COAGULATION INTEGRAL BY QUADRATURE */

    i__1 = aersl5_1.nvar;
    for (mw = aersl5_1.minw; mw <= i__1; ++mw) {
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, &
		aersl4_1.z1[mw * 24 - 24], zdist);
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, &
		aersl4_1.wmz1[mw * 24 - 24], zdist2);
	sum = 0.;
	i__2 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L4: */
	    sum += zdist[i__ - 1] * (aersl4_1.fact1a[i__ + mw * 24 - 25] * 
		    zdist2[aersl5_1.nquadp - i__] - aersl4_1.fact1b[i__ + mw *
		     24 - 25] * aersl2_1.a[mw - 1]);
	}
/* L5: */
	dadt[mw] = sum;
    }

/*        COMPUTE THE SECOND COAGULATION INTEGRAL BY QUADRATURE */

    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, &
		aersl4_1.z2[mw * 24 - 24], zdist);
	sum = 0.;
	i__2 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	    sum += aersl4_1.fact2[i__ + mw * 24 - 25] * zdist[i__ - 1];
	}
/* L35: */
	dadt[mw] -= sum * aersl2_1.a[mw - 1];
    }

/*        IF CONDENSATION IS INCLUDED COMPUTE THE DERIVATIVE WITH RESPECT */
/*        TO W OF THE RATE OF CHANGE OF W WITH RESPECT TO TIME, TIMES THE */
/*        'M' DISTRIBUTION FUNCTION.  IN TERMS OF THE PARTICLE SIZE */
/*        X, PARTICLE SIZE GROWTH RATE 'RATE',THE DERIVATIVE OF 'RATE' */
/*        WITH RESPECT TO PARTICLE SIZE 'PARTAL', THE DISTRIBUTION */
/*        FUNCTION M(W,TIME), AND THE DERIVATIVE OF M(W,TIME) WITH */
/*        RESPECT TO W, DM/DW, THE EQUIVALENT EXPRESSION IS: */

/*                RATE       DM               (         RATE) */
/*            ------------ * --  +  M(W,TIME)*(PARTAL - ----) */
/*            X  *  ALGXBA   DW               (          X  ) */

    if (aersl5_1.icondn == 0) {
	goto L30;
    }
    growth_(&aersl5_1.nvar, &aersl6_1.dumtim, aersl2_1.x, zdist, zdist2);
    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
/* L33: */
	dadt[mw] -= aersl2_1.a[mw - 1] * (zdist2[mw - 1] - zdist[mw - 1] / 
		aersl2_1.x[mw - 1]);
    }
    derivt_(&aersl5_1.isplin, &aersl5_1.nvar, zdist2);
    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
/* L34: */
	dadt[mw] -= zdist[mw - 1] * zdist2[mw - 1] / (aersl2_1.x[mw - 1] * 
		aersl6_1.algxba);
    }
    dadt[1] += boundr_(&aersl6_1.dumtim) * aersl6_1.algxba * aersl2_1.x[0];

/*       COMPUTE THE EFFECT OF ANY SOURCE OR REMOVAL MECHANISMS */

L30:
    if (aersl5_1.ncsorc + aersl5_1.ncremv == 0) {
	goto L16;
    }
    i__1 = aersl5_1.nvar;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L31: */
	zdist[i__ - 1] = aersl2_1.a[i__ - 1] / (aersl6_1.algxba * aersl2_1.x[
		i__ - 1]);
    }
    csourc_(&aersl5_1.nvar, &aersl6_1.dumtim, aersl2_1.x, zdist, source, 
	    remove);
    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
	sum = 0.;
	if (aersl5_1.ncremv == 0) {
	    goto L17;
	}
	i__2 = aersl5_1.ncremv;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L18: */
	    sum += remove[i__ + mw * 3 - 4];
	}
L17:
	if (aersl5_1.ncsorc == 0) {
	    goto L19;
	}
	i__2 = aersl5_1.ncsorc;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    sum += source[i__ + mw * 3 - 4];
	}
L19:
	dadt[mw] += aersl6_1.algxba * aersl2_1.x[mw - 1] * sum;
    }
L16:
    if (aersl5_1.ndiscr == 0) {
	return 0;
    }

/*        COMPUTE THE NET RATE OF GROWTH IN THE CONTINUOUS REGIME DUE TO */
/*        DISCRETE CLUSTERS COAGULATING WITH CONTINUOUS PARTICLES.  TO */
/*        REDUCE ROUND OFF ERROR START WITH THE LARGEST DISCRETE PARTICLE */
/*        SINCE IT PROBABLY HAS THE SMALLEST EFFECT */

    i__1 = aersl5_1.nvar;
    for (mw = 1; mw <= i__1; ++mw) {
	if (mw > 1) {
	    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl4_1.maxgdc[mw - 
		    1], &aersl4_1.wd[mw * 20 - 20], zdist);
	}
	i__2 = aersl5_1.ndiscr;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = aersl5_1.ndiscr - ii + 1;
	    term = -aersl2_1.a[mw - 1] * aersl4_1.betadc[i__ + mw * 20 - 21];
	    if (i__ <= aersl4_1.maxgdc[mw - 1]) {
		term += aersl4_1.besr[i__ + mw * 20 - 21] * zdist[i__ - 1];
	    }
/* L8: */
	    dadt[mw] += aersl2_1.a[i__ + aersl5_1.nvar - 1] * term;
	}
    }

/*       COMPUTE PARTICLE FORMATION RATE DUE TO TWO DISCRETE CLUSTERS */
/*       COAGULATING TO FORM A PARTICLE IN THE DISCRETE REGIME */

    i__2 = aersl5_1.nmulp1;
    for (mw = 1; mw <= i__2; ++mw) {
	sum = 0.;
	msmls = aersl5_1.ndcls[mw - 1] - aersl5_1.ndiscr;
	msize = aersl5_1.ndcls[mw - 1] / 2;
	i__1 = msize;
	for (i__ = msmls; i__ <= i__1; ++i__) {
	    j = aersl5_1.ndcls[mw - 1] - i__;
	    term = aersl4_1.betdsr[j + (i__ - 1) * aersl5_1.ndiscr - i__ * (
		    i__ - 1) / 2 - 1] * aersl2_1.a[i__ + aersl5_1.nvar - 1] * 
		    aersl2_1.a[j + aersl5_1.nvar - 1] * (1. - (doublereal) (
		    i__ / j) * .5);
/* L14: */
	    sum += term;
	}
/* L15: */
	dadt[mw] += sum * aersl2_1.x[i__ - 1] * aersl6_1.algxba / aersl2_1.x[
		aersl5_1.nvarp1 - 1];
    }

/*        COMPUTE THE NET RATE OF LOSS DUE TO EVAPORATION OF MONOMER */

    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.maxevp, 
	    aersl4_1.wpls, zdist);
    i__2 = aersl5_1.maxevp;
    for (mw = 1; mw <= i__2; ++mw) {
/* L32: */
	dadt[mw] += aersl4_1.evppls[mw - 1] * zdist[mw - 1] - aersl4_1.evpstr[
		mw - 1] * aersl2_1.a[mw - 1];
    }

/*                      ******************* */
/*       THE REMAINDER OF THE SUBROUTINE DETERMINES THE RATE OF CHANGE OF */
/*       THE DISCRETE CLUSTER CONCENTRATIONS WITH RESPECT TO TIME */
/*                      ******************* */

    if (aersl5_1.issflg == 1) {
	return 0;
    }

/*        COMPUTE THE NET RATE OF MONOMER FORMATION DUE TO EVAPORATION */
/*        FROM CLUSTERS IN THE DISCRETE REGIME */

    sum = 0.;
    nvarp3 = aersl5_1.nvar + 3;
    if (nvarp3 > aersl5_1.nvart) {
	goto L28;
    }
    i__2 = aersl5_1.nvart;
    for (i__ = nvarp3; i__ <= i__2; ++i__) {
/* L2: */
	sum += (aersl4_1.evpstr[i__ - 1] - aersl4_1.betdsr[i__ - 
		aersl5_1.nvar - 1] * aersl2_1.a[aersl5_1.nvarp1 - 1]) * 
		aersl2_1.a[i__ - 1];
    }

/*        COMPUTE THE NET RATE OF MONOMER FORMATION DUE TO EVAPORATION */
/*        FROM PARTICLES IN THE CONTINUOUS REGIME */

L28:
    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, aersl3_1.z__,
	     zdist);
    i__2 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L6: */
	sum += (aersl4_1.facte1[i__ - 1] - aersl4_1.fact2[i__ + 
		aersl5_1.nvarp1 * 24 - 25] * aersl2_1.a[aersl5_1.nvarp1 - 1]) 
		* zdist[i__ - 1];
    }

/*        COMPUTE THE EFFECT OF ANY SOURCES AND OR REMOVAL MECHANISMS */

    if (aersl5_1.ndsorc + aersl5_1.ndremv == 0) {
	goto L21;
    }
    dsourc_(&aersl5_1.ndiscr, &aersl6_1.dumtim, &aersl2_1.a[aersl5_1.nvarp1 - 
	    1], source, remove);
    if (aersl5_1.ndremv == 0) {
	goto L22;
    }
    i__2 = aersl5_1.ndremv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L23: */
	sum += remove[i__ - 1];
    }
L22:
    if (aersl5_1.ndsorc == 0) {
	goto L21;
    }
    i__2 = aersl5_1.ndsorc;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L24: */
	sum += source[i__ - 1];
    }
L21:
/* Computing 2nd power */
    d__1 = aersl2_1.a[aersl5_1.nvarp1 - 1];
    dadt[aersl5_1.nvarp1] = sum + (aersl4_1.evpstr[aersl5_1.nvarp2 - 1] + 
	    aersl4_1.evpstr[aersl5_1.nvarp2 - 1] - aersl4_1.betdsr[1] * 
	    aersl2_1.a[aersl5_1.nvarp1 - 1]) * aersl2_1.a[aersl5_1.nvarp2 - 1]
	     - aersl4_1.betdsr[0] * (d__1 * d__1);

/*        COMPUTE THE NET RATE OF CHANGE OF THE DISCRETE CLUSTERS DUE TO */
/*        COAGUALTION WITHIN THE DISCRETE REGIME.  TO REDUCE ROUND OFF */
/*        ERROR, START WITH THE LARGEST CLUSTER */

    i__2 = aersl5_1.nvart;
    for (mw = aersl5_1.nvarp2; mw <= i__2; ++mw) {
	mwmnvr = mw - aersl5_1.nvar;
	sum = 0.;
	maxsz = mwmnvr / 2;
	i__1 = aersl5_1.ndiscr;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__ = aersl5_1.ndiscr - ii + 1;
	    term = 0.;
	    if (i__ <= maxsz) {
		i__3 = mwmnvr - i__;
		term = aersl4_1.betdsr[i__3 + (i__ - 1) * aersl5_1.ndiscr - 
			i__ * (i__ - 1) / 2 - 1] * aersl2_1.a[mw - i__ - 1] * 
			(1. - (doublereal) (i__ / (mwmnvr - i__)) * .5);
	    }
/* L9: */
	    i__3 = min(i__,mwmnvr);
	    i__4 = max(i__,mwmnvr);
	    sum += aersl2_1.a[i__ + aersl5_1.nvar - 1] * (term - 
		    aersl4_1.betdsr[i__4 + (i__3 - 1) * aersl5_1.ndiscr - 
		    i__3 * (i__3 - 1) / 2 - 1] * aersl2_1.a[mw - 1]);
	}

/*        COMPUTE THE LOSS OF CLUSTERS DUE TO COAGULATION WITH PARTICLES */
/*        IN THE CONTINUOUS REGIME */

	sum2 = 0.;
	i__3 = aersl5_1.nquadp;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L11: */
	    sum2 += aersl4_1.fact2[i__ + mw * 24 - 25] * zdist[i__ - 1];
	}

/*        COMPUTE THE LOSS OF CLUSTERS DUE TO EVAPORATION.  FOR THE */
/*        LARGEST CLUSTER INTEGRATE THE CONTINUOUS DISTRIBUTION FROM THE */
/*        SMALLEST CONTINUOUS PARTICLE TO THE SMALLEST CONTINUOUS */
/*        PARTICLE PLUS A MONOMER.  SINCE A VERY SMALL REGION IS BEING */
/*        INTEGRATED ONLY A FOUR POINT QUADRATURE IS REQUIRED */

	sum -= aersl2_1.a[mw - 1] * (aersl4_1.evpstr[mw - 1] + sum2);
	if (mw < aersl5_1.nvart) {
	    goto L7;
	}
	distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &c__4, aersl4_1.ze, zdist);
	for (i__ = 1; i__ <= 4; ++i__) {
/* L12: */
	    sum += aersl4_1.facte3[i__ - 1] * zdist[i__ - 1];
	}
	goto L29;
L7:
	sum += aersl2_1.a[mw] * aersl4_1.evpstr[mw];
L29:
	if (aersl5_1.ndsorc + aersl5_1.ndremv == 0) {
	    goto L1;
	}

/*        COMPUTE THE EFFECT OF ANY SOURCE AND/OR REMOVAL MECHANISMS */

	if (aersl5_1.ndremv == 0) {
	    goto L25;
	}
	i__3 = aersl5_1.ndremv;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L26: */
	    sum += remove[i__ + mwmnvr * 3 - 4];
	}
L25:
	if (aersl5_1.ndsorc == 0) {
	    goto L1;
	}
	i__3 = aersl5_1.ndsorc;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L27: */
	    sum += source[i__ + mwmnvr * 3 - 4];
	}
L1:
	dadt[mw] = sum;
    }
    return 0;
} /* diffun_ */

/* Subroutine */ int steady_(doublereal *time)
{
    /* Format strings */
    static char fmt_7[] = "(\0020*** STEADY STATE CONVERGENCE TROUBLE, CONTA\
CT FRED GELBARD\002/10x,\002ITERATION\002,i5,3x,\002TIME\002,e12.4/(8e15.5))";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal b, c__;
    static integer i__, j, im1, ip1;
    static doublereal sum, atry[40], evapr, eflux;
    static integer jrevs;
    extern /* Subroutine */ int distw_(integer *, integer *, integer *, 
	    doublereal *, doublereal *);
    static doublereal zdist[24];
    static integer maxsz;
    static doublereal scavng[40], discrm, remove[60]	/* was [3][20] */, 
	    source[60]	/* was [3][20] */;
    extern /* Subroutine */ int dsourc_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer jcount, imnvar, ierror;
    static doublereal strsrc[40];

    /* Fortran I/O blocks */
    static cilist io___168 = { 0, 0, 0, fmt_7, 0 };


/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS IN THE  * */
/* *    DISCRETE REGIME USING SUCESSIVE SUBSTITUTION AS THE ITERATION    * */
/* *    TECHNIQUE                                                        * */
/* ----------------------------------------------------------------------- */


/*        COMPUTE THE VALUE OF THE M DISTRIBUTION AT THE I-TH QUADRATURE */
/*        POINT AND STORE THE RESULT IN ZDIST(I). */

    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &aersl5_1.nquadp, aersl3_1.z__,
	     zdist);

/*        COMPUTE THE SCAVENGING RATE OF THE I-TH CLUSTER */
/*        *** IF THERE ARE REMOVAL MECHANISMS IN THE DISCRETE REGIME */
/*            THEY MUST BE FIRST ORDER *** */

    if (aersl5_1.ndsorc + aersl5_1.ndremv > 0) {
	dsourc_(&aersl5_1.ndiscr, time, &aersl2_1.a[aersl5_1.nvarp1 - 1], 
		source, remove);
    }
    i__1 = aersl5_1.ndiscr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	scavng[i__ - 1] = 0.;
	strsrc[i__ - 1] = 0.;
	if (aersl5_1.ndsorc == 0) {
	    goto L9;
	}
	sum = 0.;
	i__2 = aersl5_1.ndsorc;
	for (j = 1; j <= i__2; ++j) {
/* L10: */
	    sum += source[j + i__ * 3 - 4];
	}
	strsrc[i__ - 1] = sum;
	sum = 0.;
L9:
	if (aersl5_1.ndremv == 0) {
	    goto L8;
	}
	i__2 = aersl5_1.ndremv;
	for (j = 1; j <= i__2; ++j) {
/* L11: */
	    sum += remove[j + i__ * 3 - 4];
	}
	sum = -sum / aersl2_1.a[i__ + aersl5_1.nvar - 1];
L8:
	i__2 = aersl5_1.nquadp;
	for (j = 1; j <= i__2; ++j) {
/* L1: */
	    sum += aersl4_1.fact2[j + (i__ + aersl5_1.nvar) * 24 - 25] * 
		    zdist[j - 1];
	}
/* L21: */
	scavng[i__ - 1] = sum;
    }

/*        COMPUTE THE MONOMER EVAPORATION RATE FROM THE CONTINUOUS REGIME */

    evapr = 0.;
    i__1 = aersl5_1.nquadp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L22: */
	evapr += aersl4_1.facte1[i__ - 1] * zdist[i__ - 1];
    }
    distw_(&aersl5_1.isplin, &aersl5_1.nvarm1, &c__4, aersl4_1.ze, zdist);
    eflux = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L23: */
	eflux += aersl4_1.facte3[i__ - 1] * zdist[i__ - 1];
    }
    for (jcount = 1; jcount <= 50; ++jcount) {

/*       DETERMINE THE MONOMER CONCENTRATION */

	b = scavng[0];
	c__ = evapr + aersl4_1.evpstr[aersl5_1.nvarp2 - 1] * aersl2_1.a[
		aersl5_1.nvarp2 - 1] + strsrc[0];
	i__1 = aersl5_1.nvart;
	for (i__ = aersl5_1.nvarp2; i__ <= i__1; ++i__) {
	    b += aersl4_1.betdsr[i__ - aersl5_1.nvar - 1] * aersl2_1.a[i__ - 
		    1];
/* L6: */
	    c__ += aersl4_1.evpstr[i__ - 1] * aersl2_1.a[i__ - 1];
	}
	discrm = b * b + aersl4_1.betdsr[0] * 4. * c__;
	if (discrm < 0.) {
	    goto L12;
	}
	atry[0] = (-b + sqrt(discrm)) / (aersl4_1.betdsr[0] * 2.);

/*       DETERMINE THE CLUSTER CONCENTRATIONS FROM THE DIMER TO THE */
/*       LARGEST CLUSTER IN THE DISCRETE REGIME. */

	i__1 = aersl5_1.nvart;
	for (i__ = aersl5_1.nvarp2; i__ <= i__1; ++i__) {
	    imnvar = i__ - aersl5_1.nvar;
	    b = scavng[imnvar - 1] + aersl4_1.evpstr[i__ - 1];
	    im1 = i__ - 1;
	    i__2 = im1;
	    for (j = aersl5_1.nvarp1; j <= i__2; ++j) {
/* L14: */
		i__3 = j - aersl5_1.nvar;
		b += aersl4_1.betdsr[imnvar + (i__3 - 1) * aersl5_1.ndiscr - 
			i__3 * (i__3 - 1) / 2 - 1] * aersl2_1.a[j - 1];
	    }
	    if (i__ == aersl5_1.nvart) {
		goto L17;
	    }
	    ip1 = i__ + 1;
	    i__3 = aersl5_1.nvart;
	    for (j = ip1; j <= i__3; ++j) {
/* L16: */
		i__2 = j - aersl5_1.nvar;
		b += aersl4_1.betdsr[i__2 + (imnvar - 1) * aersl5_1.ndiscr - 
			imnvar * (imnvar - 1) / 2 - 1] * aersl2_1.a[j - 1];
	    }
	    c__ = aersl4_1.evpstr[ip1 - 1] * aersl2_1.a[ip1 - 1] + strsrc[
		    imnvar - 1];
	    goto L18;
L17:
	    c__ = eflux + strsrc[imnvar - 1];
L18:
	    maxsz = imnvar / 2;
	    i__2 = maxsz;
	    for (jrevs = 1; jrevs <= i__2; ++jrevs) {
		j = maxsz - jrevs + 1;
/* L19: */
		i__3 = imnvar - j;
		c__ += aersl4_1.betdsr[i__3 + (j - 1) * aersl5_1.ndiscr - j * 
			(j - 1) / 2 - 1] * aersl2_1.a[j + aersl5_1.nvar - 1] *
			 aersl2_1.a[i__ - j - 1] * (1. - (doublereal) (j / (
			imnvar - j)) * .5);
	    }
	    discrm = b * b + aersl4_1.betdsr[imnvar + (imnvar - 1) * 
		    aersl5_1.ndiscr - imnvar * (imnvar - 1) / 2 - 1] * 4. * 
		    c__;
	    if (discrm < 0.) {
		goto L12;
	    }
/* L15: */
	    atry[imnvar - 1] = (-b + sqrt(discrm)) / (aersl4_1.betdsr[imnvar 
		    + (imnvar - 1) * aersl5_1.ndiscr - imnvar * (imnvar - 1) /
		     2 - 1] * 2.);
	}

/*       CHECK FOR CONVERGENCE OF THE STEADY STATE PROFILE */

	ierror = 0;
	i__1 = aersl5_1.ndiscr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if ((d__1 = aersl2_1.a[i__ + aersl5_1.nvar - 1] - atry[i__ - 1], 
		    abs(d__1)) * aersl6_3.relerr > atry[i__ - 1] + 
		    aersl6_3.abserr) {
		goto L3;
	    }
/* L2: */
	}
	goto L4;
L3:
	ierror = 1;
L4:
	i__1 = aersl5_1.ndiscr;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	    aersl2_1.a[i__ + aersl5_1.nvar - 1] = atry[i__ - 1];
	}
	if (ierror == 0) {
	    return 0;
	}
/* L25: */
    }
L12:
    io___168.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___168);
    do_fio(&c__1, (char *)&jcount, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*time), (ftnlen)sizeof(doublereal));
    i__1 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&aersl2_1.a[i__ - 1], (ftnlen)sizeof(doublereal)
		);
    }
    i__3 = aersl5_1.ndiscr;
    for (i__ = 1; i__ <= i__3; ++i__) {
	do_fio(&c__1, (char *)&atry[i__ - 1], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* steady_ */

/* Subroutine */ int check_(integer *istop)
{
    /* Initialized data */

    static integer maxvar = 40;
    static integer maxdis = 20;
    static integer maxtim = 16;
    static integer maxsrc = 3;
    static integer maxrmv = 3;
    static integer maxdsr = 3;
    static integer maxdrm = 3;

    /* Format strings */
    static char fmt_30[] = "(\0020---SMALLEST SIZE MUST BE LESS THAN LARGEST\
 SIZE\002)";
    static char fmt_31[] = "(\0020---SMALLEST SIZE MUST BE POSITIVE\002)";
    static char fmt_26[] = "(\0020---MAXIMUM NUMBER OF GRID POINTS IS GREATE\
R THAN\002,i4)";
    static char fmt_33[] = "(\0020---NUMBER OF DISCRETE SIZES IS GREATER T\
HAN\002,i4)";
    static char fmt_45[] = "(\0020---OUTPUT TIMES MUST BE NON-NEGATIVE AND I\
N ASCENDING ORDER\002)";
    static char fmt_20[] = "(\0020---NUMBER OF OUTPUT TIMES IS GREATER THA\
N\002,i4)";
    static char fmt_36[] = "(\0020---NUMBER OF CONTINUOUS SOURCES IS GREATER\
 THAN\002,i4)";
    static char fmt_40[] = "(\0020---NUMBER OF CONTINUOUS REMOVALS IS GREATE\
R THAN\002,i4)";
    static char fmt_44[] = "(\0020---NUMBER OF DISCRETE SOURCES IS GREATER T\
HAN\002,i4)";
    static char fmt_49[] = "(\0020---NUMBER OF DISCRETE REMOVALS IS GREATER \
THAN \002,i4)";
    static char fmt_42[] = "(\0020---IF THERE ARE DISCRETE SIZES AT LEAST 1 \
MULTIPLET MUST BE GIVEN\002)";
    static char fmt_11[] = "(\0020---IF CONDENSATION IS INCLUDED THERE CAN B\
E NO DISCRETE SIZES\002)";
    static char fmt_35[] = "(\0020---CONDENSATION HAS NOT BEEN PROPERLY SPEC\
IFIED\002)";
    static char fmt_51[] = "(\0020---WHETHER TO CALL SETUP HAS NOT BEEN PROP\
ERLY SPECIFIED\002)";
    static char fmt_62[] = "(\0020---IF COAGULATION COEFFICIENT IS TIME DEPE\
NDENT HAS NOT BEEN PROPERLY SPECIFIED\002)";
    static char fmt_54[] = "(\0020---GROWTH RATE IS NEGATIVE\002)";
    static char fmt_13[] = "(\0020---NUMERICAL DERIVATIVE OF GROWTH RATE FOR\
 SMALLEST PARTICLE AT INITIAL TIME DOES NOT MATCH GIVEN DERIVATIVE\002/9x\
,\002NUMERICAL AND GIVEN DERIVATIVES RESPECTIVELY ARE\002,2e15.5)";
    static char fmt_58[] = "(\0020---INITIAL DISTRIBUTION NEGATIVE FOR A PAR\
TICLE SIZE=\002,e12.4)";
    static char fmt_7[] = "(\0020---BETA ROUTINE IS NOT SYMMETRIC\002)";
    static char fmt_4[] = "(\0020---BETA ROUTINE IS NOT POSITIVE\002)";
    static char fmt_61[] = "(\0020---FOR PROBLEMS WITH NO COAGULATION IT IS \
BETTER TO USE THE METHOD OF CHARACTERISTICS\002/\002   REF. J. COLLOID AND I\
NTERFACE SCIENCE, V.68, P.173 (1979)\002)";
    static char fmt_56[] = "(\0020---EVAP ROUTINE IS NEGATIVE\002)";
    static char fmt_60[] = "(\0020---IF THE DISCRETE CONTINUOUS EQUATIONS AR\
E USED, THERE MUST BE AT LEAST 2 DISCRETE SIZES\002)";
    static char fmt_47[] = "(\0020---MULTIPLET NUMBER\002,i4,\002 HAS A CONC\
ENTRATION OF\002,e13.4)";
    static char fmt_29[] = "(\0020---GRID POINT NUMBER\002,i3,\002 HAS A MUL\
TIPLET NUMBER LESS THAN OR EQUAL TO THE NUMBER OF DISCRETE SIZES + 1\002)";
    static char fmt_8[] = "(\0020---MULTIPLET NUMBER\002,i5,\002 IS GREATER \
THAN TWICE THE NUMBER OF DISCRETE SIZES\002)";
    static char fmt_25[] = "(\0020---MULTIPLET NUMBERS\002,i4,\002 AND\002,i\
4,\002 ARE OUT OF ORDER\002)";
    static char fmt_16[] = "(\0020---GRID POINTS MUST BE BETWEEN THE SMALLES\
T AND LARGEST PARTICLE SIZE\002)";
    static char fmt_39[] = "(\0020---NUMBER OF GRID POINTS GIVEN MUST BE LES\
S THAN OR EQUAL TO THE NUMBER OF TOTAL GRID POINTS\002/9x,\002MINUS THE NUMB\
ER OF MULTIPLETS GIVEN MINUS 2\002)";
    static char fmt_14[] = "(\0020---GRID POINTS MUST BE IN ORDER, CHECK GRI\
D POINT NO.\002,i4)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(), do_fio(integer *, char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal p1, r1, r2, r3;
    extern /* Subroutine */ int beta_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), evap_(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal xmean;
    extern /* Subroutine */ int growth_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___176 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___177 = { 0, 0, 0, fmt_31, 0 };
    static cilist io___178 = { 0, 0, 0, fmt_26, 0 };
    static cilist io___179 = { 0, 0, 0, fmt_33, 0 };
    static cilist io___181 = { 0, 0, 0, fmt_45, 0 };
    static cilist io___182 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___183 = { 0, 0, 0, fmt_36, 0 };
    static cilist io___184 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___185 = { 0, 0, 0, fmt_44, 0 };
    static cilist io___186 = { 0, 0, 0, fmt_49, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_42, 0 };
    static cilist io___188 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_35, 0 };
    static cilist io___190 = { 0, 0, 0, fmt_51, 0 };
    static cilist io___191 = { 0, 0, 0, fmt_62, 0 };
    static cilist io___197 = { 0, 0, 0, fmt_54, 0 };
    static cilist io___198 = { 0, 0, 0, fmt_13, 0 };
    static cilist io___199 = { 0, 0, 0, fmt_58, 0 };
    static cilist io___200 = { 0, 0, 0, fmt_7, 0 };
    static cilist io___201 = { 0, 0, 0, fmt_4, 0 };
    static cilist io___202 = { 0, 0, 0, fmt_61, 0 };
    static cilist io___203 = { 0, 0, 0, fmt_56, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___205 = { 0, 0, 0, fmt_47, 0 };
    static cilist io___206 = { 0, 0, 0, fmt_29, 0 };
    static cilist io___207 = { 0, 0, 0, fmt_8, 0 };
    static cilist io___208 = { 0, 0, 0, fmt_25, 0 };
    static cilist io___209 = { 0, 0, 0, fmt_16, 0 };
    static cilist io___210 = { 0, 0, 0, fmt_39, 0 };
    static cilist io___211 = { 0, 0, 0, fmt_14, 0 };


/* ----------------------------------------------------------------------- */
/* *    THIS SUBROUTINE CHECKS THE USER INPUT FOR OBVIOUS ERRORS AND     * */
/* *    PRINTS AN ERROR MESSAGE IF ANY ERRORS ARE ENCOUNTERED.  THE      * */
/* *    ROUTINE WILL RUN THROUGH ALL ERROR CHECKS REGARDLESS OF THE      * */
/* *    NUMBER OF ERRORS.  IF ANY ERRORS ARE FOUND ISTOP WILL BE SET     * */
/* *    TO 1 AND THE CASE WILL BE SKIPPED IN THE MAIN PROGRAM            * */
/* ----------------------------------------------------------------------- */
    *istop = 0;
/*                 ********************** */
/*                    CHECK INPUT DATA */
/*                 ********************** */
    if (aersl6_2.xa < aersl6_2.xb) 
    {
	goto L1;
    }
    io___176.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___176);
    e_wsfe();
    *istop = 1;
L1:
    if (aersl6_2.xa > 0.) 
    {
	goto L2;
    }
    io___177.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___177);
    e_wsfe();
    *istop = 1;
L2:
    if (aersl5_1.nmxvar <= maxvar) 
    {
	goto L21;
    }
    *istop = 1;
    io___178.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___178);
    do_fio(&c__1, (char *)&maxvar, (ftnlen)sizeof(integer));
    e_wsfe();
L21:
    if (aersl5_1.ndiscr <= maxdis) 
    {
	goto L28;
    }
    *istop = 1;
    io___179.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___179);
    do_fio(&c__1, (char *)&maxdis, (ftnlen)sizeof(integer));
    e_wsfe();
L28:
    if (aersl5_1.nprnt == 1) 
    {
	goto L9;
    }
    i__1 = aersl5_1.nprnt;
    for (i__ = 2; i__ <= i__1; ++i__) 
    {
		
	if (aersl6_2.tout[i__ - 2] > aersl6_2.tout[i__ - 1] || aersl6_2.tout[i__ - 2] < 0.) 
        {
	    goto L3;
	}
/* L43: */
    }
    goto L9;
L3:
    *istop = 1;
    io___181.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___181);
    e_wsfe();
L9:
    if (aersl5_1.nprnt <= maxtim) 
    {
	goto L18;
    }
    *istop = 1;
    io___182.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___182);
    do_fio(&c__1, (char *)&maxtim, (ftnlen)sizeof(integer));
    e_wsfe();
L18:
    if (aersl5_1.ncsorc <= maxsrc) {
	goto L34;
    }
    *istop = 1;
    io___183.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___183);
    do_fio(&c__1, (char *)&maxsrc, (ftnlen)sizeof(integer));
    e_wsfe();
L34:
    if (aersl5_1.ncremv <= maxrmv) {
	goto L37;
    }
    *istop = 1;
    io___184.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___184);
    do_fio(&c__1, (char *)&maxrmv, (ftnlen)sizeof(integer));
    e_wsfe();
L37:
    if (aersl5_1.ndsorc <= maxdsr) {
	goto L41;
    }
    *istop = 1;
    io___185.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___185);
    do_fio(&c__1, (char *)&maxdsr, (ftnlen)sizeof(integer));
    e_wsfe();
L41:
    if (aersl5_1.ndremv <= maxdrm) {
	goto L48;
    }
    *istop = 1;
    io___186.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___186);
    do_fio(&c__1, (char *)&maxdrm, (ftnlen)sizeof(integer));
    e_wsfe();
L48:
    if (aersl5_1.ndiscr == 0 && aersl5_1.nmulp1 == 1 || aersl5_1.ndiscr > 0 &&
	     aersl5_1.nmulp1 > 1) {
	goto L19;
    }
    io___187.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___187);
    e_wsfe();
    *istop = 1;
L19:
    if (aersl5_1.icondn == 1 && aersl5_1.ndiscr == 0 || aersl5_1.icondn == 0) 
	    {
	goto L6;
    }
    *istop = 1;
    io___188.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___188);
    e_wsfe();
L6:
    if (aersl5_1.icondn == 0 || aersl5_1.icondn == 1) {
	goto L32;
    }
    *istop = 1;
    io___189.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___189);
    e_wsfe();
L32:
    if (aersl5_1.isetup == 0 || aersl5_1.isetup == 1) {
	goto L52;
    }
    *istop = 1;
    io___190.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___190);
    e_wsfe();
L52:
    if (aersl5_1.itimdp == 0 || aersl5_1.itimdp == 1) {
	goto L24;
    }
    *istop = 1;
    io___191.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___191);
    e_wsfe();

/*          ************************************* */
/*             CHECK USER SUPPLIED SUBROUTINES */
/*          ************************************* */

L24:
    xmean = sqrt(aersl6_2.xa * aersl6_2.xb);
    if (aersl5_1.icondn == 0) {
	goto L10;
    }
    growth_(&c__1, &aersl6_2.tout[aersl5_1.nprnt - 1], &xmean, &r3, &p1);
    d__1 = aersl6_2.xa * 1.00001;
    growth_(&c__1, aersl6_2.tout, &d__1, &r2, &p1);
    growth_(&c__1, aersl6_2.tout, &aersl6_2.xa, &r1, &p1);
    if (r1 >= (float)0. && r2 >= (float)0. && r3 >= 0.) {
	goto L53;
    }
    *istop = 1;
    io___197.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___197);
    e_wsfe();
L53:
    r3 = (r2 - r1) * 1e5 / aersl6_2.xa;
    if ((d__1 = r3 - p1, abs(d__1)) * (float)100. < abs(p1)) {
	goto L10;
    }
    *istop = 1;
    io___198.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___198);
    do_fio(&c__1, (char *)&r3, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&p1, (ftnlen)sizeof(doublereal));
    e_wsfe();
L10:
    if (aersl6_2.tout[0] > 0.) {
	goto L57;
    }
    i__1 = aersl5_1.nvart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (aersl2_1.a[i__ - 1] < 0.) {
	    goto L64;
	}
/* L63: */
    }
    goto L57;
L64:
    *istop = 1;
    io___199.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___199);
    do_fio(&c__1, (char *)&aersl2_1.x[i__ - 1], (ftnlen)sizeof(doublereal));
    e_wsfe();
L57:
    beta_(&c__1, aersl6_2.tout, &aersl6_2.xa, &aersl6_2.xb, &r1);
    beta_(&c__1, aersl6_2.tout, &aersl6_2.xb, &aersl6_2.xa, &r2);
    if ((d__1 = r1 - r2, abs(d__1)) * (float)1e3 <= abs(r1)) {
	goto L50;
    }
    *istop = 1;
    io___200.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___200);
    e_wsfe();
L50:
    if (r1 > (float)0. && r2 > 0.) {
	goto L17;
    }
    *istop = 1;
    io___201.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___201);
    e_wsfe();
L17:
    if (r1 == 0. || r2 == 0.) {
	io___202.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___202);
	e_wsfe();
    }

/*      ************************************************************** */
/*      CHECK QUANTITIES ASSOCIATED WITH THE DISCRETE-CONTINUOUS G-D-E */
/*      ************************************************************** */

    if (aersl5_1.ndiscr == 0) {
	goto L5;
    }
    evap_(&c__1, aersl6_2.tout, &xmean, &r1);
    evap_(&c__1, aersl6_2.tout, &aersl6_2.xa, &r2);
    evap_(&c__1, aersl6_2.tout, &aersl6_2.xb, &r3);
    if (r1 >= (float)0. && r2 >= (float)0. && r3 >= 0.) {
	goto L55;
    }
    io___203.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___203);
    e_wsfe();
L55:
    if (aersl5_1.ndiscr >= 2) {
	goto L59;
    }
    *istop = 1;
    io___204.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___204);
    e_wsfe();
L59:
    i__1 = aersl5_1.nvart;
    for (i__ = aersl5_1.nvarp1; i__ <= i__1; ++i__) {
	if (aersl2_1.a[i__ - 1] >= 0.) {
	    goto L46;
	}
	*istop = 1;
	io___205.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___205);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aersl2_1.a[i__ - 1], (ftnlen)sizeof(doublereal)
		);
	e_wsfe();
L46:
	;
    }
    i__1 = aersl5_1.nmulp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (aersl5_1.ndcls[i__ - 1] > aersl5_1.ndiscr + 1) {
	    goto L27;
	}
	*istop = 1;
	io___206.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___206);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L27:
	if (aersl5_1.ndcls[i__ - 1] <= aersl5_1.ndiscr << 1) {
	    goto L22;
	}
	*istop = 1;
	io___207.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___207);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L22:
	if (aersl5_1.ndcls[i__ - 2] < aersl5_1.ndcls[i__ - 1]) {
	    goto L23;
	}
	*istop = 1;
	io___208.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___208);
	i__2 = i__ - 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L23:
	;
    }

/*               ***************** */
/*               CHECK GRID POINTS */
/*               ***************** */

L5:
    if (aersl5_1.igrid == 0 && aersl5_1.ndiscr == 0) {
	return 0;
    }
    if (aersl2_1.w[1] > (float)0. && aersl2_1.w[aersl5_1.nvarm1 - 1] < 1.) {
	goto L15;
    }
    io___209.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___209);
    e_wsfe();
    *istop = 1;
L15:
    if (abs(aersl5_1.igrid) <= aersl5_1.nvarm1 - aersl5_1.nmulp1) {
	goto L38;
    }
    io___210.ciunit = aersl5_1.iprnt;
    s_wsfe(&io___210);
    e_wsfe();
    *istop = 1;
L38:
    i__1 = aersl5_1.nvarm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (aersl2_1.w[i__ - 1] > aersl2_1.w[i__ - 2]) {
	    goto L12;
	}
	*istop = 1;
	io___211.ciunit = aersl5_1.iprnt;
	s_wsfe(&io___211);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	e_wsfe();
L12:
	;
    }
    return 0;
} /* check_ */

/* Subroutine */ int distw_(integer *isplin, integer *nvarm1, integer *npts, doublereal *wx, doublereal *zdist)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, jstart;

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

    /* Parameter adjustments */
    --zdist;
    --wx;

    /* Function Body */
    jstart = 1;
    i__1 = *npts;
    for (i__ = 1; i__ <= i__1; ++i__) 
    {
	i__2 = *nvarm1;
	for (j = jstart; j <= i__2; ++j) 
	{
	    if (wx[i__] <= aersl2_1.w[j]) 
	    {
		goto L20;
	    }
/* L25: */
	}
	j = *nvarm1;
L20:
	if (*isplin == 1) 
	{
	    goto L1;
	}
	if (*isplin == 2) 
	{
	    goto L2;
	}

/*            CUBIC SPLINES */

	if (j == 1) 
	{
	    goto L30;
	}
	if (j == *nvarm1) 
	{
	    goto L40;
	}
/* Computing 3rd power */
	d__1 = aersl2_1.w[j] - wx[i__];
/* Computing 3rd power */
	d__2 = wx[i__] - aersl2_1.w[j - 1];
/* Computing 2nd power */
	d__3 = aersl2_1.w[j] - aersl2_1.w[j - 1];
/* Computing 2nd power */
	d__4 = aersl2_1.w[j] - aersl2_1.w[j - 1];
	zdist[i__] = (aersl7_1.p[j - 2] * (d__1 * (d__1 * d__1)) + aersl7_1.p[j - 1] * (d__2 * (d__2 * d__2)) + (aersl2_1.a[j] * 6. - d__3 *
		 d__3 * aersl7_1.p[j - 1]) * (wx[i__] - aersl2_1.w[j - 1]) + (aersl2_1.a[j - 1] * 6. - d__4 * d__4 * aersl7_1.p[j - 2]) * (
		aersl2_1.w[j] - wx[i__])) / ((aersl2_1.w[j] - aersl2_1.w[j - 1]) * 6.);
	goto L3;
L30:
/* Computing 3rd power */
	d__1 = wx[i__] - aersl2_1.w[j - 1];
/* Computing 2nd power */
	d__2 = aersl2_1.w[j] - aersl2_1.w[j - 1];
	zdist[i__] = (aersl7_1.p[j - 1] * (d__1 * (d__1 * d__1)) + (
		aersl2_1.a[j] * 6. - d__2 * d__2 * aersl7_1.p[j - 1]) * (wx[
		i__] - aersl2_1.w[j - 1]) + aersl2_1.a[j - 1] * 6. * (
		aersl2_1.w[j] - wx[i__])) / ((aersl2_1.w[j] - aersl2_1.w[j - 
		1]) * 6.);
	goto L3;
L40:
/* Computing 3rd power */
	d__1 = aersl2_1.w[j] - wx[i__];
/* Computing 2nd power */
	d__2 = aersl2_1.w[j] - aersl2_1.w[j - 1];
	zdist[i__] = (aersl7_1.p[j - 2] * (d__1 * (d__1 * d__1)) + aersl2_1.a[
		j] * 6. * (wx[i__] - aersl2_1.w[j - 1]) + (aersl2_1.a[j - 1] *
		 6. - d__2 * d__2 * aersl7_1.p[j - 2]) * (aersl2_1.w[j] - wx[
		i__])) / ((aersl2_1.w[j] - aersl2_1.w[j - 1]) * 6.);
	goto L3;

/*            LINEAR SPLINES */

L1:
	zdist[i__] = aersl2_1.a[j] + (wx[i__] - aersl2_1.w[j]) / (aersl2_1.w[
		j - 1] - aersl2_1.w[j]) * (aersl2_1.a[j - 1] - aersl2_1.a[j]);
	goto L3;

/*            LOGARITHMIC SPLINES */

L2:
	zdist[i__] = log(aersl2_1.a[j]);
	zdist[i__] += (wx[i__] - aersl2_1.w[j]) / (aersl2_1.w[j - 1] - 
		aersl2_1.w[j]) * (log(aersl2_1.a[j - 1]) - zdist[i__]);
	zdist[i__] = exp(zdist[i__]);
L3:
	jstart = j;
    }
    return 0;
} /* distw_ */

/* Subroutine */ int derivt_(integer *isplin, integer *nvar, doublereal *
	deriv)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, nvarm1, nvarm2;

/* ----------------------------------------------------------------------- */
/* *    THIS FUNCTION CALCULATES THE DERIVATIVE OF THE 'M' DISTRIBUTION  * */
/* *    FUNCTION W.R.T. W, FOR THE NVAR GRID POINTS                      * */
/* *      REF. SEE REFERENCES TO ROUTINE DISTW                           * */
/* *                                                                     * */
/* *              VARIABLES                                              * */
/* *    A(I)='M' DISTRIBUTION AT W(I)                                    * */
/* *    DERIV=OUTPUT ARRAY OF DERIVATIVES AT GRID POINTS                 * */
/* *    NVAR=NUMBER OF GRID POINTS  (INPUT)                              * */
/* *    W(I)=LOGARITHMICALLY TRANSFORMED PARTICLE SIZE AT THE I-TH GRID  * */
/* *      POINT                                                          * */
/* *    X(I)=PARTICLE SIZE AT THE I-TH GRID POINT                        * */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --deriv;

    /* Function Body */
    nvarm1 = *nvar - 1;
    nvarm2 = *nvar - 2;
    deriv[1] = (aersl2_1.a[1] - aersl2_1.a[0]) / (aersl2_1.w[c__1] - 
	    aersl2_1.w[c__1 - 1]) - (aersl2_1.w[c__1] - aersl2_1.w[c__1 - 1]) 
	    * aersl7_1.p[0] / 6.;
    i__1 = nvarm2;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L1: */
	deriv[i__] = (aersl2_1.a[i__] - aersl2_1.a[i__ - 1]) / (aersl2_1.w[
		i__] - aersl2_1.w[i__ - 1]) - (aersl2_1.w[i__] - aersl2_1.w[
		i__ - 1]) * (aersl7_1.p[i__ - 2] * 2. + aersl7_1.p[i__ - 1]) /
		 6.;
    }
    deriv[nvarm1] = (aersl2_1.a[*nvar - 1] - aersl2_1.a[nvarm1 - 1]) / (
	    aersl2_1.w[nvarm1] - aersl2_1.w[nvarm1 - 1]) - (aersl2_1.w[nvarm1]
	     - aersl2_1.w[nvarm1 - 1]) * aersl7_1.p[nvarm2 - 1] / 3.;
    deriv[*nvar] = deriv[nvarm1] + (aersl2_1.w[nvarm1] - aersl2_1.w[nvarm1 - 
	    1]) * .5 * aersl7_1.p[nvarm2 - 1];
    return 0;
} /* derivt_ */

/* Subroutine */ int setup_(integer *isplin, integer *newset, integer *nvar)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal g[38];
    static integer i__, j;
    static doublereal sub[38], diag[38], super[38];
    static integer nvarm2, nvarm3;

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


/*        COMPUTE TRI-DIAGONAL MATRIX */

    if (*newset == 0) {
	goto L6;
    }
    nvarm2 = *nvar - 2;
    nvarm3 = nvarm2 - 1;
    i__1 = nvarm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = aersl2_1.w[i__] - aersl2_1.w[i__ - 1];
	sub[i__ - 1] = d__1 * d__1 / 6.;
	i__2 = i__ + 1;
	super[i__ - 1] = (aersl2_1.w[i__] - aersl2_1.w[i__ - 1]) * (
		aersl2_1.w[i__2] - aersl2_1.w[i__2 - 1]) / 6.;
/* L1: */
	diag[i__ - 1] = (sub[i__ - 1] + super[i__ - 1]) * 2.;
    }
    i__1 = nvarm2;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L2: */
	diag[i__ - 1] -= sub[i__ - 1] * super[i__ - 2] / diag[i__ - 2];
    }

/*        SOLVE THE LINEAR SYSTEM WITH THE TRI-DIAGONAL MATRIX */

L6:
    i__1 = nvarm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	i__2 = i__ + 1;
	aersl7_1.p[i__ - 1] = (aersl2_1.w[i__] - aersl2_1.w[i__ - 1]) * (
		aersl2_1.a[i__ + 1] - aersl2_1.a[i__]) / (aersl2_1.w[i__2] - 
		aersl2_1.w[i__2 - 1]) - aersl2_1.a[i__] + aersl2_1.a[i__ - 1];
    }
    g[0] = aersl7_1.p[0] / diag[0];
    i__2 = nvarm2;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* L3: */
	g[i__ - 1] = (aersl7_1.p[i__ - 1] - sub[i__ - 1] * g[i__ - 2]) / diag[
		i__ - 1];
    }
    aersl7_1.p[nvarm2 - 1] = g[nvarm2 - 1];
    i__2 = nvarm3;
    for (j = 1; j <= i__2; ++j) {
	i__ = nvarm2 - j;
/* L4: */
	aersl7_1.p[i__ - 1] = g[i__ - 1] - super[i__ - 1] * aersl7_1.p[i__] / 
		diag[i__ - 1];
    }
    return 0;
} /* setup_ */

/* Subroutine */ int csourc_(integer *npts, doublereal *time, doublereal *v, 
	doublereal *distv, doublereal *source, doublereal *remove)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;

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
    /* Parameter adjustments */
    remove -= 4;
    source -= 4;
    --distv;
    --v;

    /* Function Body */
    i__1 = *npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = log(v[i__] / c_b565) / c_b566;
	source[i__ * 3 + 1] = .3989 * c_b564 * exp(-.5 * (d__1 * d__1)) / (v[
		i__] * c_b566);
/* Computing 2nd power */
	d__1 = log(v[i__] / c_b568) / c_b569;
	source[i__ * 3 + 2] = .3989 * c_b564 * exp(-.5 * (d__1 * d__1)) / (v[
		i__] * c_b569);
	remove[i__ * 3 + 1] = v[i__] * -5e6 * exp(*time * -.005) * distv[i__];
/* L1: */
	remove[i__ * 3 + 2] = v[i__] * -6e4 * exp(*time * -.01) * distv[i__];
    }
    return 0;
} /* csourc_ */

/* Subroutine */ int dsourc_(integer *npts, doublereal *time, doublereal *
	clscon, doublereal *source, doublereal *remove)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;

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
    /* Parameter adjustments */
    remove -= 4;
    source -= 4;
    --clscon;

    /* Function Body */
    i__1 = *npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	remove[i__ * 3 + 1] = clscon[i__] * -.01;
	remove[i__ * 3 + 2] = exp(-((doublereal) i__)) * -.1 * clscon[i__] * 
		exp(*time * -.01);
	source[i__ * 3 + 1] = 0.;
	source[i__ * 3 + 2] = exp(-((doublereal) i__)) * 1e3;
	if ((doublereal) i__ > 1.) {
	    goto L1;
	}
	source[i__ * 3 + 1] = 1e7;
L1:
	;
    }
    return 0;
} /* dsourc_ */

/* Subroutine */ int evap_(integer *npts, doublereal *time, doublereal *v, 
	doublereal *coef)
{
    /* Initialized data */

    static doublereal surten = 70.;
    static doublereal vmolec = 6.54e-23;
    static doublereal boltz = 1.38e-16;
    static doublereal tempk = 298.;
    static doublereal eqnum = 1e6;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), exp(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__;
    extern /* Subroutine */ int beta_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

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
    /* Parameter adjustments */
    --coef;
    --v;

    /* Function Body */
    i__1 = *npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	coef[i__] = 0.;
	if ((doublereal) aersl1_1.iparm[0] == 1.) {
	    goto L1;
	}
	d__1 = v[i__] * 1.9099;
	d__ = pow_dd(&d__1, &c_b574);
	d__1 = v[i__] - vmolec;
	beta_(&c__1, time, &vmolec, &d__1, &coef[i__]);
	coef[i__] = coef[i__] * eqnum * exp(surten * 4. * vmolec / (d__ * 
		boltz * tempk));
L1:
	;
    }
    return 0;
} /* evap_ */

/* Subroutine */ int errchk_(integer *nchars, integer *narray)
{
    /* Format strings */
    static char fmt_5[] = "(/\002 ERRCHK WAS CALLED INCORRECTLY.\002)";
    static char fmt_10[] = "(\002 ERRCHK MESSAGE LIMIT REACHED.\002)";
    static char fmt_15[] = "(/\002 *** NOTE ***\002/\002 TO MAKE THE ERROR M\
ESSAGE PRINTED ABOVE BE NONFATAL,\002/\002 OR TO SUPPRESS THE MESSAGE COMPLE\
TELY,\002/\002 INSERT AN APPROPRIATE CALL TO ERXSET\002,\002 AT THE START OF\
 YOUR PROGRAM.\002/\002 FOR EXAMPLE, TO PRINT UP TO 10 NONFATAL WARNING MESS\
AGES, USE\002/\002          CALL ERXSET(10,0)\002)";
    static char fmt_20[] = "(/\002 PROGRAM ABORT DUE TO ERROR.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nf, nt;
    extern /* Subroutine */ int errget_(integer *, integer *), erxset_(
	    integer *, integer *), errprt_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___237 = { 0, 6, 0, fmt_5, 0 };
    static cilist io___238 = { 0, 6, 0, fmt_10, 0 };
    static cilist io___239 = { 0, 6, 0, fmt_15, 0 };
    static cilist io___240 = { 0, 6, 0, fmt_20, 0 };



/*     SANDIA MATHEMATICAL PROGRAM LIBRARY */
/*     APPLIED MATHEMATICS DIVISION 2642 */
/*     SANDIA LABORATORIES */
/*     ALBUQUERQUE, NEW MEXICO 87115 */

/*     SIMPLIFIED VERSION FOR STAND-ALONE USE.     APRIL 1977 */

/*     ABSTRACT */
/*         THE ROUTINES ERRCHK, ERXSET, AND ERRGET TOGETHER PROVIDE */
/*         A UNIFORM METHOD WITH SEVERAL OPTIONS FOR THE PROCESSING */
/*         OF DIAGNOSTICS AND WARNING MESSAGES WHICH ORIGINATE */
/*         IN THE MATHEMATICAL PROGRAM LIBRARY ROUTINES. */
/*         ERRCHK IS THE CENTRAL ROUTINE, WHICH ACTUALLY PROCESSES */
/*         MESSAGES. */

/*     DESCRIPTION OF ARGUMENTS */
/*         NCHARS - NUMBER OF CHARACTERS IN HOLLERITH MESSAGE. */
/*                  IF NCHARS IS NEGATED, ERRCHK WILL UNCONDITIONALLY */
/*                  PRINT THE MESSAGE AND STOP EXECUTION.  OTHERWISE, */
/*                  THE BEHAVIOR OF ERRCHK MAY BE CONTROLLED BY */
/*                  AN APPROPRIATE CALL TO ERXSET. */
/*         NARRAY - NAME OF ARRAY OR VARIABLE CONTAINING THE MESSAGE, */
/*                  OR ELSE A LITERAL HOLLERITH CONSTANT CONTAINING */
/*                  THE MESSAGE.  BY CONVENTION, ALL MESSAGES SHOULD */
/*                  BEGIN WITH *IN SUBNAM, ...*, WHERE SUBNAM IS THE */
/*                  NAME OF THE ROUTINE CALLING ERRCHK. */

/*     EXAMPLES */
/*         1. TO ALLOW CONTROL BY CALLING ERXSET, USE */
/*            CALL ERRCHK(30,30HIN QUAD, INVALID VALUE OF ERR.) */
/*         2. TO UNCONDITIONALLY PRINT A MESSAGE AND STOP EXECUTION, USE */
/*            CALL ERRCHK(-30,30HIN QUAD, INVALID VALUE OF ERR.) */


    /* Parameter adjustments */
    --narray;

    /* Function Body */
    errget_(&nf, &nt);
/*     IF ERRCHK WAS CALLED WITH NEGATIVE CHARACTER COUNT, SET FATAL FLAG */
    if (*nchars < 0) {
	nf = -1;
    }
/*     IF MESSAGES ARE TO BE SUPPRESSED, RETURN */
    if (nf == 0) {
	return 0;
    }
/*     IF CHARACTER COUNT IS INVALID, STOP */
    if (*nchars == 0) {
	s_wsfe(&io___237);
	e_wsfe();
    }
    if (*nchars == 0) {
	s_stop("", (ftnlen)0);
    }
/*     PRINT MESSAGE */
    i__1 = abs(*nchars);
    errprt_(&i__1, &narray[1]);
/*     IF LAST MESSAGE, SAY SO */
    if (nf == 1) {
	s_wsfe(&io___238);
	e_wsfe();
    }
/*     PRINT TRACE-BACK IF ASKED TO */
/*     IF ((NT.GT.0).OR.(NF.LT.0)) CALL SYSTEM ROUTINE FOR TRACEBACK */
/*     DECREMENT MESSAGE COUNT */
    if (nf > 0) {
	--nf;
    }
    erxset_(&nf, &nt);
/*     IF ALL IS WELL, RETURN */
    if (nf >= 0) {
	return 0;
    }
/*     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL, */
/*     THEN EXPLAIN ERXSET USAGE. */
    if (*nchars > 0) {
	s_wsfe(&io___239);
	e_wsfe();
    }
    s_wsfe(&io___240);
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* errchk_ */

/* Subroutine */ int onechk_(integer *nchars, integer *narray)
{
    /* Initialized data */

    static struct {
	char e_1[4];
	integer e_2;
	} equiv_244 = { {'.', '$', ',', '*'}, 0 };

#define nflag (*(integer *)&equiv_244)


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer nf, nt;
    extern /* Subroutine */ int errget_(integer *, integer *), erxset_(
	    integer *, integer *), errprt_(...);


/*     ABSTRACT */
/*         ONECHK IS A COMPANION ROUTINE OF ERRCHK.  IT IS CALLED */
/*         JUST LIKE ERRCHK, AND MESSAGES FROM IT MAY BE SUPPRESSED */
/*         BY AN APPROPRIATE CALL TO ERXSET.  IT DIFFERS FROM ERRCHK */
/*         IN THAT EACH CALL TO ONECHK WILL PRODUCE NO MORE THAN ONE */
/*         PRINTED MESSAGE, REGARDLESS OF HOW MANY TIMES THAT CALL IS */
/*         EXECUTED, AND ONECHK NEVER TERMINATES EXECUTION. */
/*         ITS PURPOSE IS TO PROVIDE ONE-TIME-ONLY INFORMATIVE */
/*         DIAGNOSTICS. */

/*     DESCRIPTION OF ARGUMENTS */
/*         NCHARS - NUMBER OF CHARACTERS IN THE MESSAGE. */
/*                  IF NEGATED, THE MESSAGE WILL BE PRINTED (ONCE) EVEN */
/*                  IF NFATAL HAS BEEN SET TO 0 (SEE ERXSET). */
/*         NARRAY - SAME AS IN ERRCHK */

    /* Parameter adjustments */
    --narray;

    /* Function Body */
    if (narray[1] == nflag) {
	return 0;
    }
    errget_(&nf, &nt);
    if (nf == 0 && *nchars > 0) {
	return 0;
    }
    errprt_(&c__59, "THE FOLLOWING INFORMATIVE DIAGNOSTIC WILL APPEAR ONLY O\
NCE.");
    i__1 = abs(*nchars);
    errprt_(&i__1, &narray[1]);
    if (nf > 0) {
	--nf;
    }
    erxset_(&nf, &nt);
    narray[1] = nflag;
    return 0;
} /* onechk_ */

#undef nflag


/* Subroutine */ int errprt_(integer *nchars, integer *narray)
{
    /* Format strings */
    static char fmt_1[] = "(1x,13a10)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer i__, nch, nwords;

    /* Fortran I/O blocks */
    static cilist io___247 = { 0, 6, 0, fmt_1, 0 };



/*     UTILITY ROUTINE TO SIMPLY PRINT THE HOLLERITH MESSAGE IN NARRAY, */
/*     WHOSE LENGTH IS NCHARS CHARACTERS. */


/*     NOTE - NCH MUST BE THE NUMBER OF HOLLERITH CHARACTERS STORED */
/*     PER WORD.  IF NCH IS CHANGED, FORMAT 1 MUST ALSO BE */
/*     CHANGED CORRESPONDINGLY. */

    /* Parameter adjustments */
    --narray;

    /* Function Body */
    nch = 10;
/*     FOR LINE PRINTERS, USE */
/* L1: */
/*     FOR DATA TERMINALS, USE */
/*   1 FORMAT (1X,7A10) */
    nwords = (*nchars + nch - 1) / nch;
    s_wsfe(&io___247);
    i__1 = nwords;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&narray[i__], (ftnlen)sizeof(integer));
    }
    e_wsfe();
    return 0;
} /* errprt_ */

/* Subroutine */ int erxset_(integer *nfatal, integer *ntrace)
{
    extern /* Subroutine */ int erstgt_(integer *, integer *, integer *);


/*     ABSTRACT */
/*         ERXSET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK. */
/*         ERXSET ASSIGNS THE VALUES OF NFATAL AND NTRACE RESPECTIVELY */
/*         TO NF AND NT IN COMMON BLOCK MLBLK0 THEREBY SPECIFYING THE */
/*         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK. */

/*     DESCRIPTION OF ARGUMENTS */
/*         BOTH ARGUMENTS ARE INPUT ARGUMENTS OF DATA TYPE INTEGER. */
/*         NFATAL - IS A FATAL-ERROR / MESSAGE-LIMIT FLAG. A NEGATIVE */
/*                  VALUE DENOTES THAT DETECTED DIFFICULTIES ARE TO BE */
/*                  TREATED AS FATAL ERRORS.  NONNEGATIVE MEANS NONFATAL. */
/*                  A NONNEGATIVE VALUE IS THE MAXIMUM NUMBER OF NONFATAL */
/*                  WARNING MESSAGES WHICH WILL BE PRINTED BY ERRCHK, */
/*                  AFTER WHICH NONFATAL MESSAGES WILL NOT BE PRINTED. */
/*                  (DEFAULT VALUE IS -1.) */
/*         NTRACE - .GE.1 WILL CAUSE A TRACE-BACK TO BE GIVEN, */
/*                        IF THIS FEATURE IS IMPLEMENTED ON THIS SYSTEM. */
/*                  .LE.0 WILL SUPPRESS ANY TRACE-BACK, EXCEPT FOR */
/*                        CASES WHEN EXECUTION IS TERMINATED. */
/*                  (DEFAULT VALUE IS 0.) */

/*         *NOTE* -- SOME CALLS TO ERRCHK WILL CAUSE UNCONDITIONAL */
/*         TERMINATION OF EXECUTION.  ERXSET HAS NO EFFECT ON SUCH CALLS. */

/*     EXAMPLES */
/*         1. TO PRINT UP TO 100 MESSAGES AS NONFATAL WARNINGS USE */
/*            CALL ERXSET(100,0) */
/*         2. TO SUPPRESS ALL MATHLIB WARNING MESSAGES USE */
/*            CALL ERXSET(0,0) */

    erstgt_(&c__0, nfatal, ntrace);
    return 0;
} /* erxset_ */

/* Subroutine */ int errget_(integer *nfatal, integer *ntrace)
{
    extern /* Subroutine */ int erstgt_(integer *, integer *, integer *);


/*     ABSTRACT */
/*         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK. */
/*         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES */
/*         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE */
/*         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK. */

/*     DESCRIPTION OF ARGUMENTS */
/*         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER. */
/*         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.) */
/*         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.) */

    erstgt_(&c__1, nfatal, ntrace);
    return 0;
} /* errget_ */

/* Subroutine */ int erstgt_(integer *k, integer *nfatal, integer *ntrace)
{
    /* Initialized data */

    static integer lnf = -1;
    static integer lnt = 0;


/*     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS */
/*     THE FLAGS AS LOCAL VARIABLES. */

/*     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN */
/*     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE */
/*     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING */
/*     VALUES IN THE MAIN PROGRAM. */

    if (*k <= 0) {
	lnf = *nfatal;
    }
    if (*k <= 0) {
	lnt = *ntrace;
    }
    if (*k > 0) {
	*nfatal = lnf;
    }
    if (*k > 0) {
	*ntrace = lnt;
    }
    return 0;
} /* erstgt_ */

/* Subroutine */ int ode_(U_fp f, integer *neqn, doublereal *y, doublereal *t,
	 doublereal *tout, doublereal *relerr, doublereal *abserr, integer *
	iflag, doublereal *work, integer *iwork)
{
    /* Initialized data */

    static integer ialpha = 1;
    static integer ibeta = 13;
    static integer isig = 25;
    static integer iv = 38;
    static integer iw = 50;
    static integer ig = 62;
    static integer iphase = 75;
    static integer ipsi = 76;
    static integer ix = 88;
    static integer ih = 89;
    static integer ihold = 90;
    static integer istart = 91;
    static integer itold = 92;
    static integer idelsn = 93;

    extern /* Subroutine */ int de_(U_fp, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , doublereal *, doublereal *, integer *, logical *, integer *, 
	    integer *, integer *);
    static integer ip, iyp, iwt, iyy, iphi;
    static logical nornd, start, phase1;
    static integer iypout;


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
    if (abs(*iflag) < 2 || abs(*iflag) > 6) {
	goto L1;
    }
    start = work[istart] > 0.;
    phase1 = work[iphase] > 0.;
    nornd = iwork[2] != -1;
L1:
    de_((U_fp)f, neqn, &y[1], t, tout, relerr, abserr, iflag, &work[iyy], &
	    work[iwt], &work[ip], &work[iyp], &work[iypout], &work[iphi], &
	    work[ialpha], &work[ibeta], &work[isig], &work[iv], &work[iw], &
	    work[ig], &phase1, &work[ipsi], &work[ix], &work[ih], &work[ihold]
	    , &start, &work[itold], &work[idelsn], &iwork[1], &nornd, &iwork[
	    3], &iwork[4], &iwork[5]);
    work[istart] = -1.;
    if (start) {
	work[istart] = 1.;
    }
    work[iphase] = -1.;
    if (phase1) {
	work[iphase] = 1.;
    }
    iwork[2] = -1;
    if (nornd) {
	iwork[2] = 1;
    }
    return 0;
} /* ode_ */

/* Subroutine */ int de_(S_fp f, integer *neqn, doublereal *y, doublereal *t, 
	doublereal *tout, doublereal *relerr, doublereal *abserr, integer *
	iflag, doublereal *yy, doublereal *wt, doublereal *p, doublereal *yp, 
	doublereal *ypout, doublereal *phi, doublereal *alpha, doublereal *
	beta, doublereal *sig, doublereal *v, doublereal *w, doublereal *g, 
	logical *phase1, doublereal *psi, doublereal *x, doublereal *h__, 
	doublereal *hold, logical *start, doublereal *told, doublereal *
	delsgn, integer *ns, logical *nornd, integer *k, integer *kold, 
	integer *isnold)
{
    /* Initialized data */

    static doublereal fouru = 3.08e-8;
    static integer maxnum = 500;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    integer i_sign(integer *, integer *);
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer l;
    static doublereal del, eps;
    static integer isn, kle4;
    static doublereal tend;
    extern /* Subroutine */ int step1_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, integer *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, integer *, logical *);
    static logical crash, stiff;
    extern /* Subroutine */ int intrp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     doublereal *);
    static doublereal absdel, abseps;
    extern /* Subroutine */ int errchk_(integer *, integer *);
    static doublereal releps;
    static integer nostep;


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

    if (*neqn < 1) {
	goto L10;
    }
    if (*t == *tout) {
	goto L10;
    }
    if (*relerr < (float)0. || *abserr < 0.) {
	goto L10;
    }
    eps = max(*relerr,*abserr);
    if (eps <= 0.) {
	goto L10;
    }
    if (*iflag == 0) {
	goto L10;
    }
    isn = i_sign(&c__1, iflag);
    *iflag = abs(*iflag);
    if (*iflag == 1) {
	goto L20;
    }
    if (*t != *told) {
	goto L10;
    }
    if (*iflag >= 2 && *iflag <= 5) {
	goto L20;
    }
    if (*iflag == 6 && *abserr > 0.) {
	goto L20;
    }
    if (*iflag != 6 && *iflag != 7) {
	goto L10;
    }
    if (*iflag == 6) {
	errchk_(&c_n78, (integer*)"IN ODE   , IFLAG = 6 AND ABSERR NON-POSIT\
IVE ON INPUT -- EXECUTION TERMINATED.");
    }
    if (*iflag == 7) {
	errchk_(&c_n54, (integer*)"IN ODE   , IFLAG = 7 ON INPUT -- EXECUTIO\
N TERMINATED.", (ftnlen)54);
    }
    return 0;
L10:
    *iflag = 7;
    errchk_(&c__36, (integer*)"IN ODE   , INVALID INPUT PARAMETERS.", (ftnlen)
	    36);
    return 0;

/*   ON EACH CALL SET INTERVAL OF INTEGRATION AND COUNTER FOR NUMBER OF */
/*   STEPS.  ADJUST INPUT ERROR TOLERANCES TO DEFINE WEIGHT VECTOR FOR */
/*   SUBROUTINE  STEP1 */

L20:
    del = *tout - *t;
    absdel = abs(del);
    tend = *t + del * 10.;
    if (isn < 0) {
	tend = *tout;
    }
    nostep = 0;
    kle4 = 0;
    stiff = FALSE_;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1) {
	goto L30;
    }
    if (*isnold < 0) {
	goto L30;
    }
    if (*delsgn * del > 0.) {
	goto L50;
    }

/*   ON START AND RESTART ALSO SET WORK VARIABLES X AND YY(*), STORE THE */
/*   DIRECTION OF INTEGRATION AND INITIALIZE THE STEP SIZE */

L30:
    *start = TRUE_;
    *x = *t;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    if ((d__1 = *x - *t, abs(d__1)) < absdel) {
	goto L60;
    }
    intrp_(x, &yy[1], tout, &y[1], &ypout[1], neqn, kold, &phi[phi_offset], &
	    psi[1]);
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   IF CANNOT GO PAST OUTPUT POINT AND SUFFICIENTLY CLOSE, */
/*   EXTRAPOLATE AND RETURN */

L60:
    if (isn > 0 || (d__1 = *tout - *x, abs(d__1)) >= fouru * abs(*x)) {
	goto L80;
    }
    *h__ = *tout - *x;
    (*f)(x, &yy[1], &yp[1]);
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    if (nostep < maxnum) {
	goto L100;
    }
    *iflag = isn << 2;
    if (stiff) {
	*iflag = isn * 5;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    for (l = 1; l <= i__1; ++l) {
	wt[l] = releps * (d__1 = yy[l], abs(d__1)) + abseps;
	if (wt[l] <= 0.) {
	    goto L140;
	}
/* L110: */
    }
    step1_((S_fp)f, neqn, &yy[1], x, h__, &eps, &wt[1], start, hold, k, kold, 
	    &crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], &
	    beta[1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd);

/*   TEST FOR TOLERANCES TOO SMALL */

    if (! crash) {
	goto L130;
    }
    *iflag = isn * 3;
    *relerr = eps * releps;
    *abserr = eps * abseps;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    if (*kold > 4) {
	kle4 = 0;
    }
    if (kle4 >= 50) {
	stiff = TRUE_;
    }
    goto L50;

/*   RELATIVE ERROR CRITERION INAPPROPRIATE */

L140:
    *iflag = isn * 6;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L150: */
	y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    errchk_(&c__97, (integer*)"IN ODE   , PURE RELATIVE ERROR IMPOSSIBLE WHE\
N SOLUTION VANISHES.  USE A NON-ZERO VALUE OF ABSERR", (ftnlen)97);
    return 0;
} /* de_ */

/* Subroutine */ int step1_(S_fp f, integer *neqn, doublereal *y, doublereal *
	x, doublereal *h__, doublereal *eps, doublereal *wt, logical *start, 
	doublereal *hold, integer *k, integer *kold, logical *crash, 
	doublereal *phi, doublereal *p, doublereal *yp, doublereal *psi, 
	doublereal *alpha, doublereal *beta, doublereal *sig, doublereal *v, 
	doublereal *w, doublereal *g, logical *phase1, integer *ns, logical *
	nornd)
{
    /* Initialized data */

    static doublereal twou = 1.54e-8;
    static doublereal fouru = 3.08e-8;
    static doublereal two[13] = { 2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,
	    2048.,4096.,8192. };
    static doublereal gstr[13] = { .5,.0833,.0417,.0264,.0188,.0143,.0114,
	    .00936,.00789,.00679,.00592,.00524,.00468 };

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, l;
    static doublereal r__;
    static integer iq, im1, km1, km2, ip1, kp1, kp2;
    static doublereal erk, err, tau, rho, sum;
    static integer nsm2, nsp1, nsp2;
    static doublereal absh, hnew;
    static integer knew;
    static doublereal xold, erkm1, erkm2, erkp1, temp1, temp2, temp3, temp4, 
	    temp5, temp6, p5eps;
    static integer ifail;
    static doublereal reali, round;
    static integer limit1, limit2;
    static doublereal realns;


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
    d__1 = fouru * abs(*x);
    *h__ = d_sign(&d__1, h__);
    return 0;
L5:
    p5eps = *eps * .5;

/*   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE */

    round = 0.;
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L10: */
/* Computing 2nd power */
	d__1 = y[l] / wt[l];
	round += d__1 * d__1;
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
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + phi_dim1] = yp[l];
	phi[l + (phi_dim1 << 1)] = 0.;
/* L20: */
/* Computing 2nd power */
	d__1 = yp[l] / wt[l];
	sum += d__1 * d__1;
    }
    sum = sqrt(sum);
    absh = abs(*h__);
    if (*eps < sum * 16. * *h__ * *h__) {
	absh = sqrt(*eps / sum) * .25;
    }
/* Computing MAX */
    d__2 = absh, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    *h__ = d_sign(&d__1, h__);
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
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    realns = (doublereal) (*ns);
    alpha[*ns] = 1. / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.;
    if (*k < nsp1) {
	goto L110;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	temp2 = psi[im1];
	psi[im1] = temp1;
	beta[i__] = beta[im1] * psi[im1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	reali = (doublereal) i__;
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
    i__1 = *k;
    for (iq = 1; iq <= i__1; ++iq) {
	temp3 = (doublereal) (iq * (iq + 1));
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
    temp4 = (doublereal) (*k * kp1);
    v[*k] = 1. / temp4;
    nsm2 = *ns - 2;
    if (nsm2 < 1) {
	goto L130;
    }
    i__1 = nsm2;
    for (j = 1; j <= i__1; ++j) {
	i__ = *k - j;
/* L125: */
	v[i__] -= alpha[j + 1] * v[i__ + 1];
    }

/*   UPDATE V(*) AND SET W(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i__1 = limit1;
    for (iq = 1; iq <= i__1; ++iq) {
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
    i__1 = kp1;
    for (i__ = nsp2; i__ <= i__1; ++i__) {
	limit2 = kp2 - i__;
	temp6 = alpha[i__ - 1];
	i__2 = limit2;
	for (iq = 1; iq <= i__2; ++iq) {
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
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	temp1 = beta[i__];
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
/* L205: */
	    phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
	}
/* L210: */
    }

/*   PREDICT SOLUTION AND DIFFERENCES */

L215:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
	phi[l + kp1 * phi_dim1] = 0.;
/* L220: */
	p[l] = 0.;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__ = kp1 - j;
	ip1 = i__ + 1;
	temp2 = g[i__];
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
	    p[l] += temp2 * phi[l + i__ * phi_dim1];
/* L225: */
	    phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
	}
/* L230: */
    }
    if (*nornd) {
	goto L240;
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
	p[l] = y[l] + tau;
/* L235: */
	phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
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
	d__1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
	erkm2 += d__1 * d__1;
L260:
/* Computing 2nd power */
	d__1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
	erkm1 += d__1 * d__1;
L265:
/* Computing 2nd power */
	d__1 = temp4 * temp3;
	erk += d__1 * d__1;
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
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 1. / beta[i__];
	ip1 = i__ + 1;
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
/* L305: */
	    phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[
		    l + ip1 * phi_dim1]);
	}
/* L310: */
    }
    if (*k < 2) {
	goto L320;
    }
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L315: */
	psi[i__ - 1] = psi[i__] - *h__;
    }

/*   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP */
/*   SIZE */

L320:
    ++ifail;
    temp2 = .5;
    if ((i__1 = ifail - 3) < 0) {
	goto L335;
    } else if (i__1 == 0) {
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
    d__1 = fouru * abs(*x);
    *h__ = d_sign(&d__1, h__);
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
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
	y[l] = p[l] + rho;
/* L405: */
	phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
    }
    goto L420;
L410:
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L415: */
	y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
    }
L420:
    (*f)(x, &y[1], &yp[1]);

/*   UPDATE DIFFERENCES FOR NEXT STEP */

    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
/* L425: */
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * 
		phi_dim1];
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
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
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L440: */
/* Computing 2nd power */
	d__1 = phi[l + kp2 * phi_dim1] / wt[l];
	erkp1 += d__1 * d__1;
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
    temp2 = (doublereal) (*k + 1);
    d__1 = p5eps / erk;
    d__2 = 1. / temp2;
    r__ = pow_dd(&d__1, &d__2);
/* Computing MAX */
    d__1 = .5, d__2 = min(.9,r__);
    hnew = absh * max(d__1,d__2);
/* Computing MAX */
    d__2 = hnew, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    hnew = d_sign(&d__1, h__);
L465:
    *h__ = hnew;
    return 0;
/*       ***     END BLOCK 4     *** */
} /* step1_ */

/* Subroutine */ int intrp_(doublereal *x, doublereal *y, doublereal *xout, 
	doublereal *yout, doublereal *ypout, integer *neqn, integer *kold, 
	doublereal *phi, doublereal *psi)
{
    /* Initialized data */

    static doublereal g[13] = { 1. };
    static doublereal rho[13] = { 1. };

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal w[13], hi;
    static integer ki, jm1;
    static doublereal eta;
    static integer kip1;
    static doublereal term, temp1, temp2, temp3, gamma;
    static integer limit1;
    static doublereal psijm1;


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

/*   INITIALIZE W(*) FOR COMPUTING G(*) */

    i__1 = ki;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = (doublereal) i__;
/* L5: */
	w[i__ - 1] = 1. / temp1;
    }
    term = 0.;

/*   COMPUTE G(*) */

    i__1 = ki;
    for (j = 2; j <= i__1; ++j) {
	jm1 = j - 1;
	psijm1 = psi[jm1];
	gamma = (hi + term) / psijm1;
	eta = hi / psijm1;
	limit1 = kip1 - j;
	i__2 = limit1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	    w[i__ - 1] = gamma * w[i__ - 1] - eta * w[i__];
	}
	g[j - 1] = w[0];
	rho[j - 1] = gamma * rho[jm1 - 1];
/* L15: */
	term = psijm1;
    }

/*   INTERPOLATE */

    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
	ypout[l] = 0.;
/* L20: */
	yout[l] = 0.;
    }
    i__1 = ki;
    for (j = 1; j <= i__1; ++j) {
	i__ = kip1 - j;
	temp2 = g[i__ - 1];
	temp3 = rho[i__ - 1];
	i__2 = *neqn;
	for (l = 1; l <= i__2; ++l) {
	    yout[l] += temp2 * phi[l + i__ * phi_dim1];
/* L25: */
	    ypout[l] += temp3 * phi[l + i__ * phi_dim1];
	}
/* L30: */
    }
    i__1 = *neqn;
    for (l = 1; l <= i__1; ++l) {
/* L35: */
	yout[l] = y[l] + hi * yout[l];
    }
    return 0;
} /* intrp_ */

#ifdef __cplusplus
	}
#endif
