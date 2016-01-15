/* Subroutine */ int errprt(int *nchars, int *narray)
{
    /* Format strings */
    static char fmt_1[] = "(1x,13a10)";

    /* System generated locals */
    int i1;

    /* Builtin functions */
    int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static int i, nch, nwords;

    /* Fortran I/O blocks */
    static cilist io_247 = { 0, 6, 0, fmt_1, 0 };



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
    s_wsfe(&io_247);
    i1 = nwords;
    for (i = 1; i <= i1; ++i) 
    {
	do_fio(&c__1, (char *)&narray[i], (ftnlen)sizeof(integer));
    }
    e_wsfe();
    return 0;
} /* errprt_ */

