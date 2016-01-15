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

/* Subroutine */ 
int onechk(int *nchars, int *narray)
{
    /* Initialized data */

    static struct 
    {
	char e_1[4];
	int e_2;
    } equiv_244 = { {'.', '$', ',', '*'}, 0 };

#define nflag (*(integer *)&equiv_244)


    /* System generated locals */
    int i1;

    /* Local variables */
    static int nf, nt;
    extern /* Subroutine */ int errget_(int *, int *), erxset_(int *, int *), errprt_(...);

    /* Parameter adjustments */
    --narray;

    /* Function Body */
    if (narray[1] == nflag) 
    {
	return 0;
    }
    errget_(&nf, &nt);
    if (nf == 0 && *nchars > 0) 
    {
	return 0;
    }
    errprt_(&c__59, "THE FOLLOWING INFORMATIVE DIAGNOSTIC WILL APPEAR ONLY ONCE.");
    i1 = abs(*nchars);
    errprt_(&i1, &narray[1]);
    if (nf > 0) 
    {
	--nf;
    }
    erxset_(&nf, &nt);
    narray[1] = nflag;
    return 0;
} /* onechk_ */

