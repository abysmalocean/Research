/*     ABSTRACT */
/*         ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK. */
/*         ERRGET ASSIGNS TO NFATAL AND NTRACE RESPECTIVELY THE VALUES */
/*         OF NF AND NT IN COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE */
/*         STATE OF THE OPTIONS WHICH CONTROL THE EXECUTION OF ERRCHK. */

/*     DESCRIPTION OF ARGUMENTS */
/*         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER. */
/*         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.) */
/*         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.) */

/* Subroutine */ int errget(int *nfatal, int *ntrace)
{
    extern /* Subroutine */ int erstgt_(int *, int *, int *);

    erstgt_(&c__1, nfatal, ntrace);
    return 0;
} /* errget_ */

