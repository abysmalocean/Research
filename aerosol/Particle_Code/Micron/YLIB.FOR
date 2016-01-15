      subroutine FNDKEY(keys,   maxkey, tokens, ntoken, tie,    loctie, 
     &                  sorted)
      integer           loctie, maxkey, ntoken
      character*(*)     keys,   tokens, tie
      logical           sorted
      dimension         keys(*), tokens(*)
************************************************************************
* FNDKEY searches an array of tokens for a key word.  The search begins 
* at the end of the array of tokens.
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           itoken
*-----------------------------------------------------------------------
      loctie = 0
      itoken = ntoken
10    if ((loctie .eq. 0) .and. (itoken .gt. 0)) then
        tie = tokens(itoken)
        call OPLOOK(maxkey, keys, 1, sorted, tie, loctie)
        itoken = itoken - 1
        go to 10
      end if
      return
*
*___end of FNDKEY___*
*
      end
      subroutine SPOPEN(nunit,  msglev, root, suffix)
      integer           nunit,  msglev
      character*(*)     root,   suffix
************************************************************************
*
* SPOPEN closes the default files and reopens the file under a 
* different name.
*
* LAST MODIFIED: 2 FEB 1988
************************************************************************
      character*12      filenm
*-----------------------------------------------------------------------
      if (msglev .gt. 0) then
        close(nunit)
        call MAKEFN(root, suffix, filenm)
        open(nunit, file = filenm, status = 'unknown')
      end if
      return
*
*___end of SPOPEN___*
*
      end
      subroutine MAKEFN(root, suffix, filenm)
      character*(*)     root
      character*(*)     suffix
      character*(*)     filenm
************************************************************************
*
* MAKEFN is used to make file names of the form "root"."suffix"
* Called by FILES.
*
* INPUTS:
* root, suffix
*
* OUTPUTS:
* filenm
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      filenm = ' '
      call RMBLK(root)
      call RMBLK(suffix)
      call MERGE(filenm,root)
      call MERGE(filenm,suffix)
      return
*
*___End of MAKEFN___*
*      
      end
      subroutine MERGE (word1,  word2)
************************************************************************
*
* MERGE appends word2 to word1 starting after the first nonblank charac-
* ter at the end of word1.  
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      integer           len1,   len2,   ifind,  lblnk,  irepl,  iword2, 
     &                  nrepl
      character*(*)     word1,  word2
*-----------------------------------------------------------------------
c
c ##find the position of the last blank in word1 (lblnk)
c
      lblnk = 1
      len1  = len(word1)
      len2  = len(word2)
      do 100 ifind = len1, 1, -1
        if (word1(ifind:ifind) .ne. ' ') then
          lblnk = ifind + 1
          go to 1000
        end if
100   continue
c
c ##append as many letters as possible from word2 to the trailing blanks
c ##of word1
c 
1000  nrepl = min0(len1,lblnk + len2 - 1)
      do 200 irepl = lblnk, nrepl
        iword2 = irepl - lblnk + 1
        word1(irepl:irepl) = word2(iword2:iword2)
200   continue
      return
c
*___end of MERGE___*
c
      end
      subroutine RMBLK(word)
      character*(*)     word
************************************************************************
*
* RMBLK removes the blank spaces in a word and places them at
* the end of the word.
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      integer           iilp,   iolp,   lword
      character*1       blank,  let
      parameter        (blank  = ' ')
*-----------------------------------------------------------------------
      lword = len(word)
      do 100 iolp = 1, lword - 1
        if (word(iolp:iolp) .eq. blank) then
          do 110 iilp = iolp, lword - 1
            let             = word(iilp+1:iilp+1)
            word(iilp:iilp) = let
110       continue
          word(lword:lword) = blank
        end if
100   continue
      return
*
*___end of RMBLK___*
*
      end
      subroutine INTCHR(num, word)
      integer           num
      character*(*)     word
************************************************************************
*
* INTCHR finds the character equivalent of an integer.  Nega-
* tive numbers are converted to positive.  If the number is too small 
* for the word, blanks are appended, and if the word is too small to
* hold the number, the rightmost digits are deleted.
*
* LAST MODIFIED:16 NOV 1988
************************************************************************
      integer           ichar,  idigit, lword,  nchar,  ndig,   NDIGIT, 
     &                  numt
*-----------------------------------------------------------------------
      word  = ' '
      numt  = abs(num)
      ndig  = NDIGIT(numt)
      lword = len(word)
      nchar = min0(ndig,lword)
      if (lword .lt. ndig) numt = numt * 10**(lword-ndig)
      do 100 ichar = nchar, 1, -1
        idigit = numt - 10*(numt/10)
        write(word(ichar:ichar), '(i1)') idigit
        numt = numt/10
100   continue
      return
*
*___end of INTCHR___*
*
      end
      integer function NDIGIT(num)
************************************************************************
*
* NDIGIT is an integer function which finds the number of digits in an
* integer.
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      integer num,numpos
      real xdig
c
      if (num .eq. 0) then
        NDIGIT = 1
      else
        numpos = abs(num)
        xdig   = alog10(real(numpos) + 0.1)
        NDIGIT = int(xdig) + 1
      end if
      return
*
*___end of NDIGIT___*
*
      end
      subroutine FLINT (slope,  xlow,   ylow,   x,      answer)
************************************************************************
*
* LNINT performs sinlgle precision linear interpolation.
*
* INPUTS:
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      double precision  answer, slope,  x,      xlow,   ylow
*-----------------------------------------------------------------------
      answer = ylow + (x - xlow) * slope
      return
*
*___end of FLINT___*
*
      end
      subroutine LNINT (xlow,   xhigh,  ylow,   yhigh,  x,      answer, 
     &                  rmin)
************************************************************************
*
* LNINT performs linear interpolation.
*
* INPUTS:
* rmin     - the smallest positive machine number
*
* LAST MODIFIED:16 NOV 1987
************************************************************************
      double precision  answer, rmin,   x,      xhigh,  xlow,   yhigh,
     &                  ylow
*-----------------------------------------------------------------------
      if (dabs(xlow - xhigh) .lt. rmin) then
        answer = 0.5d00*(ylow + yhigh)
      else
        answer = ylow + (x - xlow)*(yhigh - ylow)/(xhigh - xlow)
      end if
      return
*
*___end of LNINT___*
*
      end
      function NEWVEC  (n,      x,      incx,   y,      incy,   d,      
     &                  incd) 
      integer           incd,   incx,   incy,   n
      double precision  d,      x,      y
      logical           NEWVEC
      dimension         d(*), x(*), y(*)
************************************************************************
*
* NEWVEC determines if the vectors x and y are nearly equal.
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           i,      id,     ix,     iy
      logical           new
*-----------------------------------------------------------------------
      new = .false.
      i   = 1
      id  = 1
      ix  = 1
      iy  = 1
10    if ((.not. new) .and. (i .le. n)) then
        new   = (dabs(x(ix) - y(iy)) .gt. d(id))
        i = i + 1
        ix = ix + incx
        iy = iy + incy
        id = id + incd
        go to 10
      end if
      NEWVEC = new
      return
*
*___end of NEWVEC___*
*
      end
      subroutine NZPRN(nprt,   rmin,   x,      incx, nunit) 
************************************************************************
c
c NZPRN prints out the nonzero elements of a vector and reports which
c element is nonzero.
c
c VARIABLES:(type, input/output, array)
c
c LAST MODIFIED: 1 JAN 1988
c
************************************************************************
      integer           incx,   iwrite, ix,     nprt,   nunit
      double precision  rmin,   x
      dimension         x(*)
c-----------------------------------------------------------------------
c
c ##CHECK IF THE INPUTS ARE REASONABLE
c
      if (incx  .lt. 1)    return
      if (nprt  .lt. 1)    return
      if (nunit .lt. 1)    return
      if (rmin  .lt. 0.d+0) return
c
      ix = 1
      do 100 iwrite=1, nprt
        if (dabs(x(ix)) .gt. rmin) 
     &    write(nunit,1) iwrite, x(ix)
        ix = ix + incx
100   continue
      return
1     format(1x,'(',i4,')',2x,d11.4)
c
*___end of NZPRN___*
c
      end
      subroutine FUNIT (nunit)
************************************************************************
c
c FUNIT finds a unit number which can be used to open a new file.
c
c
c LAST MODIFIED: 7 NOV 1987
c
************************************************************************
      integer          iunit,  maxu,   minu,   nunit
      logical          conn
      parameter        (minu = 10)
      parameter        (maxu = 100)
c-----------------------------------------------------------------------
      do 100 iunit = minu, maxu
        inquire(unit = iunit, opened = conn)
        if (.not. conn) then
          nunit = iunit
          return
        end if
100   continue
      return
c
*____end of FUNIT___*
c
      end
      subroutine ISWTCH(int1,   int2)
************************************************************************
c
c ISWTCH switches the value of 2 integers
c
c LAST MODIFIED:10 NOV 1987
c
************************************************************************
      integer           int1,   int2,   inttmp
c-----------------------------------------------------------------------
      inttmp = int1
      int1  = int2
      int2   = inttmp
      return
c
*___end of ISWTCH___*
c
      end
      integer function ISUM(n, ix, incix)
************************************************************************
c
c ISUM is an integer function which computes the sum of an integer
c vector.
c
c INPUTS:
c n      - the number of elements which will be summed
c ivec   - the vector to be summed
c
c LAST MODIFIED:18 NOV 1987
c
************************************************************************
      integer          incix,  itemp,  ix,     jx,     n
      dimension        ix(*)
c-----------------------------------------------------------------------
      if (n .lt. 1) return
      itemp = ix(1)
      do 10 jx = 1 + incix, 1 + (n - 1)*incix, incix
        itemp = itemp + ix(jx)
10    continue
      ISUM = itemp
      return
c
*___end of ISUM___*
c
      end
      subroutine CHLOAD(nload,  const,  chvec,  ich)
************************************************************************
c
c CHLOAD load the character constant 'const' into the vector 'charr'.
c
c LAST MODIFIED:15 DEC 1987
c
************************************************************************
      integer           ich,    iload,  nload
      character*(*)     const,  chvec
      dimension         chvec(*)
c-----------------------------------------------------------------------
      if (nload .eq. 0) return
      do 100 iload = 1, 1 + (nload - 1)*ich, ich
        chvec(iload) = const
100   continue
      return
c
*___end of CHLOAD___*
c
      end
      subroutine TXPOSE(n, ldx, x)
************************************************************************
c
c TXPOSE takes the tranpose of a square matrix stored in the array x.
c
c INPUTS:
c n      - the number of rows (and columns) of the matrix stored in x.
c ldx    - the leading dimension of the array x
c
c LAST MODIFIED:15 DEC 1987
c
************************************************************************
      integer          i,      j,      ldx,    n
      double precision temp,   x
      dimension        x(ldx,*)
c-----------------------------------------------------------------------
      do 100 i = 1, n
        do 110 j = 1, i-1
          temp = x(i,j)
          x(i,j) = x(j,i)
          x(j,i) = temp
110     continue
100   continue
      return
c
*___end of TXPOSE___*
c
      end
      double precision function DSUM(nsum, x, incx)
************************************************************************
c
c DSUM computes the sum of the elements of the vector x
c
c INPUTS:
c nsum   - the number of elements to be summed
c
c LAST MODIFIED:15 DEC 1987
c
************************************************************************
      integer          incx,   isum,   nsum
      double precision temp,   x
      dimension        x(*)
c-----------------------------------------------------------------------
      temp = 0.d+0
      if (nsum .gt. 0) then 
        do 100 isum = 1, 1 + (nsum - 1)*incx, incx
          temp = temp + x(isum)
100     continue
      end if 
      DSUM = temp
      return
c
*___end of DSUM___*
c
      end
      double precision function DNRMDF(n,      x,      incx,   y, incy, 
     &  work)
************************************************************************
c
c DNRMDF computes the euclidean norm of the difference of two vectors 
c x and y
c
c LAST MODIFIED: 4 JAN 1988
c
************************************************************************
      integer          incx,   incy,   n
      double precision DSUM,   work,   x,      y
      dimension        work(*), x(*), y(*)
c
      call     DCOPY (n, x,      incx, work, 1)
      call     DXMY (n, y,    incy, work, 1)
      call     DPSCAL(n, 1,      work, 1,    work, 1)
      DNRMDF = dsqrt(DSUM  (n, work,   1))
c
      return
c
*___end of DNRMDF___*
c
      end
      subroutine DPSCAL(n,      ipower, x,     incx,   y,      incy )
************************************************************************
*
* DPSCAL computes  y:= y * x**ipower
*
* LAST MODIFIED:19 DEC 1987
************************************************************************
      integer           i,      incx,   incy
      integer           ipower, ix,     iy,     n
      double precision  x,      y
      dimension         x(*),   y(*)
*-----------------------------------------------------------------------
      if (n      .le. 0) return
      if (ipower .eq. 0) return
      ix = 1
      iy = 1
      do 100 i = 1, n
        y(iy) = y(iy) * x(ix)**ipower
        ix    = ix + incx
        iy    = iy + incy
100   continue
      return
*
*___End of DPSCAL___*
*
      end
      subroutine DESCAL(n,      escale, x,      incx,   y,     incy )
************************************************************************
c
c DESCAL computes  y:= y * EXP[x*escale]
c
c LAST MODIFIED:19 DEC 1987
c
************************************************************************
      integer           i,      incx,   incy,   ix,     iy,    n
      double precision  escale, x,      y
      dimension         x(*),   y(*)
c
      if (n      .le. 0)     return
      if (escale .eq. 0.d+0) return
      ix = 1
      iy = 1
      do 100 i = 1, n
        y(iy) = y(iy) * dexp( x(ix)*escale )
        ix    = ix + incx
        iy    = iy + incy
100   continue
      return
c
*___end of DESCAL___*
c
      end
      double precision function TRAP( delta, n, x, incx )
************************************************************************
c
c TRAP returns the integral of a function stored as a vector using the
c trapezoidal rule.  Delta corresponds to the spacing of the values at
c which the function is known.
c
c LAST MODIFIED:19 DEC 1987
c
************************************************************************
      integer          i,      iend,   incx,   n
      double precision delta,  x,      temp,   zero
      parameter        (zero = 0.d+0)
      dimension        x(*)
c
      if (n .le. 0)        return
      if (delta .eq. zero) return
      temp = zero
      iend = 1 + (n - 1)*incx
      do 100 i = 1, iend, incx
        temp = temp + x(i)
100   continue
      temp = temp - 0.5d+0*( x(1) + x(iend))
      TRAP = temp * delta
      return
c
*___end of TRAP___*
c
      end
      integer function IDXGT(n, imax,   i,      inci)
************************************************************************
c
c IDXGT returns the first element in the array i which is greater than
c imax.  If none of the elements are greater than imax, 0 is returned.
c
c LAST MODIFIED: 1 JAN 1988
c
************************************************************************
      integer           i,      imax,   inci,   j,      n
      dimension         i(*)
c-----------------------------------------------------------------------
c
      IDXGT = 0
      if (n .le. 0) return
      do 100 j = 1, 1 + (n - 1)* inci, inci
        if ( i(j) .gt. imax) then
          IDXGT = j
          return
        end if
100   continue
c
      return
c
*___end of IDXGT___*
c
      end
      double precision function DROUGH(n, delta, f, work)
************************************************************************
c
c DROUGH evaluates the integral of [f"(x)]^2.  The function values are
c stored as a vector, the second derivative is approximated by finite
c difference, and the integeral is approximated using the trapezoidal
c rule.
c
c LAST MODIFIED: 5 JAN 1988
c 
************************************************************************
      integer           n,      nm2
      double precision  delta,  f,      minus,  scal,   TRAP,  work,   
     &                  zero
      parameter        (zero   = 0.d+0)
      parameter        (minus  = -1.d+0)
      dimension         f(*),   work(*)
c-----------------------------------------------------------------------
      scal    = 1.d+0/delta**2
      work(1) = zero
      work(n) = zero
      nm2     = n - 2
      call      DCOPY (nm2,    f(2),    1,      work(2), 1)
      call      DSCAL (nm2,    2.d+0,   work(2), 1)
      call      DAXPY (nm2,    minus,   f,       1,    work(2), 1)
      call      DAXPY (nm2,    minus,   f(3),    1,    work(2), 1)
      call      DSCAL (nm2,    scal,    work(2), 1)
      call      DPSCAL(nm2,    1,       work(2), 1,    work(2), 1)
      DROUGH  = TRAP  (delta,  n,       work,    1)
      return
c
*___end of DROUGH___*
c
      end
      subroutine INTVPD(n,     tolpad,  x,      incx)
************************************************************************
c
c INTVPD adds padding to the interval (x(i), x(i + 1)) to increase it's 
c length to |x(i + 1) - x(i)| * (1 + tolgap)
c
c LAST MODIFIED: 8 JAN 1988
c
************************************************************************
      integer           i,      incx,   n
      double precision  expad,  tolpad, x
      dimension         x(*)
c
      if (tolpad .eq. 0.d+0) return
      do 100 i = 1, 1 + (n - 1)*incx
        expad    = (x(i + 1) - x(i))*tolpad /2.d+0
        x(i)     = x(i) - expad
        x(i + 1) = x(i + 1) + expad
100   continue
c
      return
c
c
      end
************************************************************************
c
c DABSVC performs x := abs(x)
c
c LAST MODIFIED:11 JAN 1988
c
************************************************************************
      subroutine DABSVC(n, x, incx)
      integer           i,      incx,   n
      double precision  x
      dimension         x(*)
c
      do 100 i = 1, 1 + (n - 1)*incx
        x(i) = dabs(x(i))
100   continue
c
      return
c
c
      end
      logical function IOUTSD(i, j, k)
************************************************************************
c
c IOUTSD is true if ( i < min{j, k} ) or ( i > max{j, k})
c
c LAST MODIFIED: 2 FEB 1988
c
************************************************************************
      integer         i, j, k
      IOUTSD = (i .lt. min0(j,k)) .or. (i .gt. max0(j,k))
      return
c
*___end of IOUTSD___*
c
      end
      logical function DOUTSD(x, y, z)
************************************************************************
c
c DOUTSD is true if (x < min{y, z})  or (x > max{y,z})
c
c LAST MODIFIED: 2 FEB 1988
c
************************************************************************
      double precision x, y, z
      DOUTSD = (x .lt. dmin1(y,z)) .or. (x .gt. dmax1(y,z))
      return
c
*___end of DOUTSD___*
c
      end
      subroutine FNWORD(string, jtoken, first,  mark)
************************************************************************
c
c FNWORD finds the first and last character of the "jtoken"th word
c in string.  If no tokens are found, mark is set to 0, and first is
c unchanged.
c
c VARIABLES:(type, input/output, array)
c 
c string    (C,I)   - the string scanned for tokens
c jtoken    (I,I/O) - (I) the position of the "jtoken"th token is 
c                   - requested, (O) the position of the "jtoken"th
c                   - token that was found.
c first     (I,I/O) - (I) the search for the word begins at the "first"
c                   - character in string. (O) the position of the
c                     first character of the target word.
c mark      (I,O)   - the position in string of the last character of
c                   - the target token
c
c LAST MODIFIED: 7 FEB 1988
c
************************************************************************
      integer           first,  jtoken, last,   mark,   ntoken, tfirst,
     &                  tlast,  tmark
      character*(*)     string
c
      if (first .lt. 0) first = 1
      ntoken = 0
      mark   = first - 1
      tfirst = first
      last   = len(string)
c
c ##SCAN string UNTIL tmark = 0
c
10    if (ntoken .lt. jtoken) then
        tlast  = last
        call OPSCAN(string, tfirst, tlast, tmark)
        if (tfirst .ne. 0) then
          first  = tfirst
          mark   = tmark
          tfirst = tmark  + 1
          ntoken = ntoken + 1
          go to 10
        end if
      end if
c
c ##END SCAN
c
      jtoken = ntoken
      return
*
*___end of FNWORD___*
*
      end
      subroutine RDSTR (ninp,   ioerr,  string,   ntoken, tokens)
************************************************************************
*
* RDSTR finds the next string in the file(UNIT="ninp") containing data
* and finds the tokens in the string of data.  Lines having a '#' as
* the first non-blank character are skipped.
*
* VARIABLES:(type, input/output, array)
*
* ioerr     (I,O) - the iostat varaible associated with the read 
*                   statement.  Set to a nonzero value if a reading
*                   error occurred.
* string      (I,O) - the string of data (from the data file) that was
*                   parsed
* ntoken    (I,I/O)   - (I) the number of tokens requested
*                       (O) the number of tokens found
* tokens    (C,O,A)   - character strings not separated by blanks, '=',
*                       tabs, ',', or ':'.
*
* LAST MODIFIED: 2 FEB 1988
************************************************************************
      integer           first,  ioerr,  last,   mark,   ninp,   ntoken,
     &                  tlast
      character*(*)     string
      character*(*)     tokens
      dimension         tokens(*)
      ioerr  = 0
c
c ##CHECK FOR QUICK RETURN
c
      if ((ntoken .le. 0) .or. (ninp .le. 0)) then
        string = ' '
        ntoken = 0
        return
      end if
c
c
      mark = 0
      last = len(string)
c
c ##FIND THE FIRST LINE WITH TOKENS
c
10    if ((mark .eq. 0)  .or. (string(first:first) .eq. '#')) then
        read(ninp, 1, end = 1000, err = 1000, iostat = ioerr) 
     &       string
        first  = 1
        tlast  = last
        call OPSCAN(string, first, tlast, mark)
        go to 10
      end if
c
c ##PARSE THE LINE OF DATA
c
      call OPTOKN(string, ntoken, tokens)
c
1000  if (ioerr .ne. 0) ntoken = 0
      return
1     format(a80)
*
*___end of RDSTR___*
*
      end
      subroutine CHCOPY(n,      cx,     incx,   cy,     incy)
************************************************************************
*
* CHCOPY
*
* LAST MODIFIED: 1 MAY 1988
************************************************************************
      integer           i,      incx,   incy,   ix,     iy,     n
      character*(*)     cx,     cy
      dimension         cx(*),  cy(*)
*-----------------------------------------------------------------------
      if (n .le. 0) return
      ix = 1
      iy = 1
      do 10 i = 1, n
        cy(iy) = cx(ix)
        ix     = ix + incx
        iy     = iy + incy
10    continue
      return
*
*___end of CHCOPY___*
*
      end
      double precision function RAN1(idum)
************************************************************************
*
* RAN1 is a portable random number generator that uses three linear
* congruential generators.  This random number generator was taken from
* NUMERICAL RECIPIES, Press, Flannery, Teukolsky, Vetterling, 1986,
* pp 186.
*
* LAST MODIFIED: 3 MAR 1988
************************************************************************
      integer           ia,     ia1,    ia2,    ia3,    ic,     idum,   
     &                  iff,    ic1,    ic2,    ic3,    ix,     ix1,    
     &                  ix2,    ix3,    j,      lcg,    lcgerr, m,      
     &                  m1,     m2,     m3,     mshuff
      double precision  r,      rm1,    rm2
      logical           IOUTSD
      parameter        (m1 = 259200, ia1 = 7141, ic1 = 54733)
      parameter        (m2 = 134456, ia2 = 8121, ic2 = 28411)
      parameter        (m3 = 243000, ia3 = 4561, ic3 = 51349)
      parameter        (rm1 = 1.d+0/259200.d+0)
      parameter        (rm2 = 1.d+0/134456.d+0)
      parameter        (mshuff = 97)
      dimension         r(mshuff)
      save              iff,    ix1,    ix2,    ix3,    lcgerr
      save              r
      data              iff, lcgerr /0, 0/
      lcg(ia,ix,ic,m) = mod(ia*ix + ic,m)
c
c ##INITIALIZE OR REINITIALIZE THE SEQUENCE IF idum < 0 OR
c ##IF IT IS THE FIRST CALL
c
      if ((idum .lt. 0) .or. (iff .eq. 0)) then
        iff  = 1
        idum = iabs(idum)
c
c ##SEED THE FIRST ROUTINE
c
        ix1 = mod(ic1 - idum, m1)
        ix1 = lcg(ia1, ix1, ic1, m1) 
c
c ##SEED THE second ROUTINE
c
        ix2 = mod(ix1,m2)
        ix1 = lcg(ia1, ix1, ic1, m1) 
c
c ##SEED THE second ROUTINE
c
        ix3 = mod(ix1, m3)
c
c ##FILL THE TABLE WITH RANDOM NUMBERS
c
        do 10 j = 1, mshuff
          ix1  = lcg(ia1, ix1, ic1, m1)
          ix2  = lcg(ia2, ix2, ic2, m2)
          r(j) = (dble(ix1) + dble(ix2)*rm2)*rm1
10      continue
        idum = 1
      end if 
      ix1 = lcg(ia1, ix1, ic1, m1)
      ix2  = lcg(ia2, ix2, ic2, m2)
c
c ##RANDOMLY CHOOSE AN ELEMENT IN THE TABLE
c
      ix3  = lcg(ia3, ix3, ic3, m3)
      j   = 1 + (mshuff * ix3)/m3
      if (IOUTSD(j, 1, mshuff) .and. (lcgerr .eq. 0)) then
        lcgerr = 1
        write(6,*) ' Error in random number generator.'
      end if
      RAN1 = r(j)
c
c ##RETURN THE ELEMENT IN THE TABLE
c
      r(j) = (dble(ix1) + dble(ix2)*rm2)*rm1
      return
*
*___end of RAN1___*
*
      end
      double precision function GASDEV(idum)
************************************************************************
*
* GASDEV returns normally distributed deviates with zero mean and unit 
* using RAN1 as the source of uniform deviates.  Taken from NUMERICAL 
* RECIPIES, Press, Flannery, Teukolsky, Vetterling, 1986,  pp 186.
* 
* LAST MODIFIED:20 FEB 1988
************************************************************************
      integer           idum,   iset
      double precision  fac,    gset,   one,    r,      RAN1,   two,    
     &                  v1,     v2,     zero
      parameter        (one  = 1.d+0)
      parameter        (two  = 2.d+0)
      parameter        (zero = 0.d+0)
      save              iset,   gset
      data              iset   /0/
*-----------------------------------------------------------------------
      if (iset .eq. 0) then
1       v1 = two*RAN1(idum) - one
        v2 = two*RAN1(idum) - one
        r  = v1**2 + v2**2
        if ((r .ge. one) .or. (r .le. zero)) go to 1
        fac    = dsqrt(-two*dlog(r)/r)
        gset   = v1*fac
        GASDEV = v2*fac
        iset   = 1
      else
        GASDEV = gset
        iset   = 0
      end if
      return
*
*___end of GASDEV___*
*
      end
      subroutine NOISE (n,      nseed,  sdev,   incs,   x,      incx)
      integer           incs,   incx,   n,      nseed
      double precision  sdev,   x
      dimension         sdev(*), x(*)
************************************************************************
*
* NOISE adds noise to the vector x.  The noise is assumed to be normally
* distributed with the standard deviations stored in sdev.
*
* LAST MODIFIED: 6 JAN 1988
************************************************************************
      integer           i,      ix,     is
      double precision  GASDEV
      if (n .le. 0) return
      ix = 1
      is = 1
      do 100 i = 1, n
        if (sdev(is) .ne. 0.d+0) 
     &    x(ix) = x(ix) + sdev(is)*GASDEV(nseed)
        ix = ix + incx
        is = is + incs
100   continue
      return
*
*___end of NOISE___*
*
      end
      subroutine SETMIN(n,      const,  x,      incx)
      integer           incx,   n
      double precision  const,  x
      dimension         x(*)
************************************************************************
*
* SETMIN does the following:    x := max{const, x}
*
* LAST MODIFIED: 6 JAN 1988
************************************************************************
      integer           i,      ix
*-----------------------------------------------------------------------
      if (n .le. 0) return
      ix = 1
      do 100 i = 1, n
        if (x(ix) .lt. const) x(ix) = const
        ix = ix + incx
100   continue
      return
*
*___end of SETMIN___*
*
      end
      subroutine DNLOAD(n,      nseed,  x,      incx)
      integer           n,      nseed,  incx
      double precision  x
      dimension         x(*)
************************************************************************
*
* DNLOAD makes a vector of random variables.
*
* LAST MODIFIED: 6 JAN 1988
************************************************************************
      integer           i,      ix
      double precision  GASDEV
      ix = 1
      do 100 i = 1, n
        x(ix) = GASDEV(nseed)
        ix    = ix + incx
100   continue
      return
*
*___end of DNLOAD___*
*
      end
      integer function IDSAFE(doub)
      double precision  doub
************************************************************************
* IDSAFE converts the double precision number to an integer without
* overflow.
*
* LAST MODIFIED: 6 July 1988
************************************************************************
      double precision  dimax
      parameter        (dimax = 2.147483d+09)
c                       dimax = 2^31 - 1
*-----------------------------------------------------------------------
      if (dabs(doub) .gt. dimax) then
        IDSAFE = nint(dsign(dimax,doub))
      else
        IDSAFE = nint(doub)
      end if
      return
*
*___end of IDSAFE___*
*
      end
      double precision function TV(n, x, incx)
***********************************************************************
*
* TV returns the total variation of the function stored as a vector x.
*
* LAST MODIFIED:11 JAN 1988
***********************************************************************
      integer           i,      incx,   n
      double precision  temp,   x,      zero
      dimension         x(*)
      parameter         (zero = 0.d+0)
*-----------------------------------------------------------------------
      TV   = zero
      if (n .eq. 1) return
      temp = zero
      do 100 i = 1, 1 + (n - 2)*incx
        temp = abs(x(i) - x(i + incx)) + temp
100   continue
      TV = temp
      return
*
*___end of TV___*
*
      end
      INTEGER FUNCTION IDDFMN( N, X, INCX, Y, INCY )
***********************************************************************
      INTEGER                  N, INCX, INCY
      DOUBLE PRECISION         X( * ), Y( *)

C  IDDFMN returns the smallest value of i such that
C
C     abs( x( i ) - y( i) ) = min( abs( x( j ) - y(j) ) )
C                      j
C
***********************************************************************
      INTRINSIC           DABS
      INTEGER             I     , IMAX  , IX , IY
      DOUBLE PRECISION    TEMP,  XMAX

      IF( N.LT.1 )THEN
         IDDFMN = 0
         RETURN
      END IF

      IMAX = 1
      IF( N.GT.1 )THEN
         XMAX = DABS( X( 1 ) - Y( 1 ) )
         IX   = 1
         IY   = 1
         DO 10, I = 2, N
            IX = IX + INCX
            IY = IY + INCY
            TEMP = DABS( X ( IX ) - Y ( IY ) )
            IF( XMAX.GT.TEMP )THEN
               XMAX = TEMP
               IMAX = I
            END IF
   10    CONTINUE
      END IF

      IDDFMN = IMAX
      RETURN

*     End of IDDFMN.

      END
      subroutine IDLOAD(n, ixzero, ixdelt, ix, incix)
      integer           incix,  ix,     ixdelt, ixzero, n
      dimension         ix(*)
***********************************************************************
*
* IDLOAD performs: ix(i) = ixzero + (i - 1) * ixdelt
*
* LAST MODIFIED: 1 JAN 1988
***********************************************************************
      integer           i,      ixval,  j

      if ((ixzero.eq.1) .and. (ixdelt.eq.1) .and. (incix.eq.1)) then
        do 100 i = 1, n
          ix(i) = i
100     continue
      else
        ixval = ixzero
        j     = 1
        do 200 i = 1, n
          ix(j) = ixval
          j     = j     + incix
          ixval = ixval + ixdelt
200     continue
      end if
      return
*
*___end of IDLOAD___*
*
      end
      subroutine KRICHK(fail,   iinst,  maxfev, maxkey, maxtok, mfail,  
     &                  minfev, minst,  mintg,  mpar,   mrecl,  msgscr, 
     &                  nch,    npar,   nscr,   decdiv, dpar,   dspc,   
     &                  par,    tol,    wmach,  work,   cfiles, inme,   
     &                  parnm,  tokens)
      integer           fail,   iinst,  maxfev, maxkey, maxtok, mfail,  
     &                  minfev, minst,  mintg,  mpar,   mrecl,  msgscr, 
     &                  nch,    npar,   nscr
      double precision  decdiv, dpar,   dspc,   par,    tol,    wmach,  
     &                  work
      character*(*)     cfiles, inme,   parnm,  tokens
      dimension         fail(*), nch(*), npar(*), dpar(*), dspc(*), 
     &                  par(*), tol(*), wmach(*), work(*), cfiles(*), 
     &                  inme(*), parnm(*), tokens(*)
************************************************************************
*
* KRICHK checks if the data file micron.kri is current. If micron.kri
* is not current, then it is brought up to date.
*
* LAST MODIFIED: 1 Apr 1989
************************************************************************
      integer           ichk,   ioerr,  ipar,   itmp,   iupdte, joerr,  
     &                  jupdte, minstt, mreclt, ncht,   nkri,   npart
      character*8       ctmp
      logical           NEWVEC
*-----------------------------------------------------------------------
c
c ##CHECK IF micron.kri NEEDS TO BE REMADE
c
      iupdte = 0
      joerr  = 0
      call OPNKRI(ioerr, 2, mrecl, nkri,  cfiles)
      if (ioerr .eq. 0) then
        read(nkri, rec=1, err=2000, iostat=joerr) 
     &    itmp, ncht, ctmp, npart, minstt, mreclt
        if ((minst .ne. minstt) .or. (mrecl .ne. mreclt)) iupdte = -1
      else
        iupdte = -1
      end if
c
c ##CHECK IF A SINGLE INSTRUMENT NEEDS TO BE UPDATED
c
      ichk   = 1
10    if ((iupdte .ge. 0) .and. (ichk .le. minst)) then
        jupdte = 0
        read(nkri, rec=ichk, err=2000, iostat=joerr) 
     &    itmp, ncht, ctmp, npart
        if ((ncht .ne. nch(ichk)) .or. (npart .ne. npar(ichk))) 
     &    jupdte = ichk
        if (npart .gt. 0) then
          read(nkri, rec=ichk+minst, err=2000, iostat=joerr) 
     &    (work(ipar),ipar=1,npart)  
          if (NEWVEC(npart, work, 1, par(ichk), minst, 
     &      dpar(ichk), minst)) jupdte = ichk
        end if
c
c ##CHECK IF TWO UPDATES WERE SPECIFIED
c
        if (jupdte .ne. 0) then
          if (iupdte .eq. 0) then
            iupdte = jupdte
          else
            iupdte = -1
          end if
        end if
        ichk = ichk + 1
        go to 10
      end if
c
c ##UPDATE OR MAKE
c
2000  if (joerr .ne. 0) then
        call WARNG(13, nscr, msgscr, 0, 0.d+0, ' ')
        iupdte = -1
      end if
      if (ioerr .eq. 0) close(nkri)
      if (iupdte .eq. -1) then
        call KERINF(fail,   maxfev, maxtok, mfail,  minfev, minst,  
     &              mintg,  mpar,   mrecl,  msgscr, nch,    npar,   
     &              nscr,   decdiv, dspc,   par,    tol,    wmach,  
     &              work,   cfiles, inme,   parnm,  tokens)
      else if (iupdte .gt. 0) then
        iinst = iupdte
        call UPDATE(fail,   maxfev, mfail,  minfev, minst,  mintg, 
     &              mpar,   mrecl,  msgscr, nch,    npar,   nscr,   
     &              decdiv, dspc,   par,    tol,    wmach,  work,   
     &              cfiles, inme,   parnm)
      end if 
      return
*
*___end of IDLOAD___*
*
      end
      subroutine TOKNUM(ndigit, ntoken, tokens, x,      incx)
      integer           incx,   ndigit, ntoken
      double precision  x
      character*(*)     tokens
      dimension         x(*), tokens(*)
************************************************************************
*
* TOKNUM finds numeric values in the character strings tokens, and loads
* them into the vector x.
*
* LAST MODIFIED: 1 JAN 1989
************************************************************************
      integer           idigit, ie,     itoken, ix
      logical           OPNUMB
*-----------------------------------------------------------------------
      idigit = 1
      itoken = 1
      ix     = 1
10    if ((idigit .le. ndigit) .and. (itoken .le. ntoken)) then
        if (OPNUMB(tokens(itoken))) then
          ie = index(tokens(itoken),'E')
          if (ie .ne. 0) tokens(itoken)(ie:ie) = 'D'
          read(tokens(itoken), '(bn, d16.0)') x(ix)
          idigit = idigit + 1
          ix     = ix     + incx
        end if
        itoken = itoken + 1
        go to 10
      end if
      ndigit = idigit - 1
      return
*
*___end of IDLOAD___*
*
      end
