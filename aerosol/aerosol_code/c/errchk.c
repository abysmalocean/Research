
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
/* Subroutine */ 
int errchk(int *nchars, int *narray)
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
    int i1;

    /* Builtin functions */
    int s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static int nf, nt;
    extern /* Subroutine */ int errget_(int *, int *), erxset_(int *, int *), errprt_(int *, int *);

    /* Fortran I/O blocks */
    static cilist io_237 = { 0, 6, 0, fmt_5, 0 };
    static cilist io_238 = { 0, 6, 0, fmt_10, 0 };
    static cilist io_239 = { 0, 6, 0, fmt_15, 0 };
    static cilist io_240 = { 0, 6, 0, fmt_20, 0 };

    /* Parameter adjustments */
    --narray;

    /* Function Body */
    errget_(&nf, &nt);
/*     IF ERRCHK WAS CALLED WITH NEGATIVE CHARACTER COUNT, SET FATAL FLAG */
    if (*nchars < 0) 
    {
	nf = -1;
    }
/*     IF MESSAGES ARE TO BE SUPPRESSED, RETURN */
    if (nf == 0) 
    {
	return 0;
    }
/*     IF CHARACTER COUNT IS INVALID, STOP */
    if (*nchars == 0) 
    {
	s_wsfe(&io_237);
	e_wsfe();
    }
    if (*nchars == 0) 
    {
	s_stop("", (ftnlen)0);
    }
/*     PRINT MESSAGE */
    i1 = abs(*nchars);
    errprt_(&i1, &narray[1]);
/*     IF LAST MESSAGE, SAY SO */
    if (nf == 1) 
    {
	s_wsfe(&io_238);
	e_wsfe();
    }
/*     PRINT TRACE-BACK IF ASKED TO */
/*     IF ((NT.GT.0).OR.(NF.LT.0)) CALL SYSTEM ROUTINE FOR TRACEBACK */
/*     DECREMENT MESSAGE COUNT */
    if (nf > 0) 
    {
	--nf;
    }
    erxset_(&nf, &nt);
/*     IF ALL IS WELL, RETURN */
    if (nf >= 0) 
    {
	return 0;
    }
/*     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL, */
/*     THEN EXPLAIN ERXSET USAGE. */
    if (*nchars > 0) 
    {
	s_wsfe(&io_239);
	e_wsfe();
    }
    s_wsfe(&io_240);
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* errchk_ */

