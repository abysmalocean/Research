/*     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS */
/*     THE FLAGS AS LOCAL VARIABLES. */

/*     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN */
/*     CALLS ON THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE */
/*     PLACED IN A COMMON BLOCK AND PRESET TO THE FOLLOWING */
/*     VALUES IN THE MAIN PROGRAM. */


/* Subroutine */ 
int erstgt(int *k, int *nfatal, int *ntrace)
{
    /* Initialized data */
    static int lnf = -1;
    static int lnt = 0;
    if (*k <= 0) 
    {
	lnf = *nfatal;
    }
    if (*k <= 0) 
    {
	lnt = *ntrace;
    }
    if (*k > 0) 
    {
	*nfatal = lnf;
    }
    if (*k > 0) 
    {
	*ntrace = lnt;
    }
    return 0;
} /* erstgt */

