      subroutine SCRNM (ident,  nuni,   msglev, intgr,  doub,   word)
      double precision  doub
      character*(*)     word
************************************************************************
*
* SCRNM is used to write messages concerning MICRON's current status.
*
* INPUTS:
* ident  - message identifier
* nuni   - the unit number to which the message is written
* msglev - an integer parameter used to decide if the message should
*          be written
* intgr  - an integer which may be included in the message
* doub   - a double precision variable which may be included in the
*          message
* word  - a character which may be included as part of the message
*
* LAST MODIFIED:13 NOV 1987
************************************************************************
      integer           ident,  intgr,  iwrite, merrm,  msgchk, msglev, 
     &                  nuni
      character*17      wfail
      character*16      tword
      parameter        (wfail = ' Inversion failed')
      parameter        (merrm = 50)
      dimension         msgchk(merrm)
      data msgchk      / 30, 30, 30, 30, 40,  30,  40,  30,  20,   20, 
     1                   20,  0, 20, 30, 20,  30,  20,  10,  10,   20,
     2                   20, 30, 30, 30, 30,  30,  30,  30,  30,   20,
     3                   30, 20, 20, 20, 20,  30,   0,   0,   0,    0,
     4                   40, 40, 40,  0, 20,  20,  20,  20,  30,    0/
c
      call WRCHK (ident,  msglev, nuni,   msgchk, merrm,  iwrite)
      if (iwrite .le. 0) return
c
      if (ident  .gt. 10) go to 1000
c
      if      (ident .eq. 1) then
        write (nuni, 101)
      else if (ident .eq. 2)  then
        write (nuni, 102) word
      else if (ident .eq. 3)  then
        write (nuni, 103)
      else if (ident .eq. 4)  then
        write (nuni, 104)
      else if (ident .eq. 5)  then
        write (nuni, 105) word, doub
      else if (ident .eq. 6)  then
        write (nuni, 106) 
      else if (ident .eq. 7) then
        write (nuni, 107) doub
      else if (ident .eq. 8) then
        write (nuni, 108) 
      else if (ident .eq. 9) then
        write (nuni, 109)
      else if (ident .eq. 10) then
        write (nuni, 110) word
      end if
      return
101   format(/,2x,' Rewriting the micron.kri data file.')
102   format(/,2x,' Computing the calibration matrix for the ',a8)
103   format(/,2x,' Determining the inversion interval.        ')
104   format(/,2x,' Checking the inversion interval.        ')
105   format(/,' The ', a5, ' inversion interval diameter found by exami
     &ning the kernel',/,' function support data contained in micron.kri
     & is ', d10.3)
106   format(//, ' MICRON''s labels reset by micron.lbl.')
107   format(/,' The initial smoothing parameter was reset to ',d10.3)
*
108   format(/,
     &' It appears that the size of the inversion interval could have',/,
     &' been reduced for this  inversion.    This  would  have  saved',/, 
     &' computer time and improved solution accuracy.')
*
109   format(//,
     &' Performing an inversion; input file skipped.')
110   format(//' Performing the inversion described by the data file ',
     &a12)
c
1000  if (ident .gt. 20) go to 2000
c
      if      (ident .eq. 11) then
        write (nuni, 111) word
      else if (ident .eq. 12) then
        write (nuni, 112) word
      else if (ident .eq. 13) then
        write (nuni, 113) word
      else if (ident .eq. 14) then
        write (nuni, 114) intgr, word
      else if (ident .eq. 15) then
        write (nuni, 115)
      else if (ident .eq. 16) then
        write (nuni, 116) intgr, word
      else if (ident .eq. 17) then
        write (nuni, 117) 
      else if (ident .eq. 18) then
        write (nuni, 118)
      else if (ident .eq. 19) then
        write (nuni, 119)
      else if (ident .eq. 20) then
        write (nuni, 120) word
      end if
      return
111   format(//' Calculating the information for the data file ',
     &a12)
112   format(//,' Generating the .plt data file for the ',a8)
113   format(/,2x,' Generating data for the ',a8)
114   format(/,2x,' Generating the .plt data file for channel ', i3,/,
     &       ' of the ',a8)
115   format(//,' Generating the micron.kri data file.')
116   format(/, 2x,' Computing the micron.kri data for channel ',
     &         i3,' of the ',a8)
117   format(//,' Copying the data from the micron.kri data file and',/,
     &' placing it in micron.cpy.')
118   format(//,' MICRON is finished.',//)
119   format(//,' MICRON is finished due to error.'//)
120   format(//,' Updating the information for the ',a8,/,
     &           ' in the micron.kri data file.')
c
2000  if (ident .gt. 30) go to 3000
c
      if (ident .eq. 21) then
        write (nuni, 121) word
      else if (ident .eq. 22) then
        write (nuni, 122) word
      else if (ident .eq. 23) then
        write (nuni, 123)
      else if (ident .eq. 24) then
        write (nuni, 124)
      else if (ident .eq. 25) then
        write (nuni, 125)
      else if (ident .eq. 26) then
        write (nuni, 126)
      else if (ident .eq. 27) then
        write (nuni, 127)
      else if (ident .eq. 28) then
        write (nuni, 128)
      else if (ident .eq. 29) then
        write (nuni, 129)
      else if (ident .eq. 30) then
        write (nuni, 130)
      end if
      return
121   format(//,' Generating the simulated data described by the data fi
     &le ',a12)
122   format(/,2x,' Reading ',a12)
123   format(/,2x,' Echoing the input data.')
124   format(/,2x,' Computing the "inverse" of the error array.')
125   format(/,2x,' Determining the inputs to the least squares problem.
     &')
126   format(/,2x,' Searching for the optimal smoothing parameter.')
127   format(/,2x,' Writing the ouput file.')
128   format(/,2x,' Solving the least squares problem.')
129   format(/,2x,' Creating the information matrices.')
130   format(/,2x,' Inversion completed successfully.')
c
3000  if (ident .gt. 40) go to 4000
c
c ##SAVE 33-50 FOR PREOUT
c
      if      (ident .eq. 31) then
        write (nuni, 131)
      else if (ident .eq. 33) then
        write (nuni, 133) wfail
      else if (ident .eq. 34) then
        write (nuni, 134) wfail
      else if (ident .eq. 35) then
        write (nuni, 135) wfail
      else if (ident .eq. 36) then
        write (nuni, 136) 
      end if
      return
131   format(/,2x,' Minimum error assumed to equal zero.')
133   format(/,2x,a17,' while trying to opening I/O data files.')
134   format(/,2x,a17,' while trying to read the input data file.')
135   format(/,2x,a17,' while trying to read the input data file.')
136   format(/,2x,' Information calulations completed successfully.')
c
4000  if (ident .gt. 50) go to 5000
c
      if      (ident .eq. 41) then
        write (nuni, 141) intgr
      else if (ident .eq. 42) then
        write (nuni, 142) 
      else if (ident .eq. 43) then
        write (nuni, 143) intgr
      else if (ident .eq. 45) then
        write (nuni, 145) wfail
      else if (ident .eq. 46) then
        write (nuni, 146) wfail
      else if (ident .eq. 47) then
        write (nuni, 147) wfail
      else if (ident .eq. 48) then
        write (nuni, 148) wfail
      else if (ident .eq. 49) then
        if (intgr .eq. 0) then
          tword = 'CREATING '
        else
          tword = 'RECOVERING'
        end if
        write (nuni, 149) tword, word
      end if
141   format(/,
     &' The solution dimension determined by examining micron.kri',/,
     &' is ',i4,'.')
*
142   format(/,
     &' MICRON was unable to use the data in micron.kri to',/,
     &' determine a reasonable solution dimension.')
*
143   format(
     &' The solution dimension was reset by MICRON to ',i4,'.')
*
145   format(/,2x,a17,' because the solution is too sensitive to',/,
     &' the smoothing parameter.')
*
146   format(/,2x,a17,' because the solution is insensitive to the',/,
     &' smoothing parameter.')
*
147   format(/,2x,a17,' because the number of iterations needed',/,
     &' to find the proper smoothing parameter is greater than the',/,
     &' maximum allowed.')
*
148   format(/,2x,a17,' because numerical  error  caused an',/,
     &' inconsistency while searching for the proper smoothing',/,
     &' parameter.')
*
149   format(//,
     &1x, a10, ' the dropfile ',a10,'.')
c
5000  return
*
*___end of SCRNM___*
*
      end

      subroutine WARNG (ident,  nuni,   msglev, intgr,  doub,   char)
************************************************************************
*
* WARNG writes warning messages.
*
* INPUTS:
* ident  - the warning message identifier
* msglev - an integer used to decide if the warning message needs to be 
*          written
* nuni   - the unit number to which the message is written
*
* LAST MODIFIED: 5 JAN 1988
************************************************************************
      integer           ident,  idoub,  IDSAFE, intgr,  iwrite, merrm,  
     &                  msgchk, msglev, nuni,   warn
      double precision  doub
      character*(*)     char
      character*7       exword
      parameter        (merrm = 40)
      dimension         msgchk(merrm)
      save              warn
      data msgchk      / 30, 30,  30,  30,  30,  30,  30,  30,  30,  30,
     &                   30, 30,  30,  30,  30,  30,  30,  30,  30,  30,
     &                   30, 30,  30,  30,  30,  30,  30,  30,   0,  30,
     &                   30, 30,  30,  30,  30,  30,  90,  90,  90,  90/
*-----------------------------------------------------------------------
      iwrite = 0
      if (ident .eq. 0) then
	warn = 0
      else if (ident .eq. -1) then
        intgr = warn
      else
        warn  = warn + 1
        call WRCHK (ident,  msglev, nuni,   msgchk, merrm,  iwrite)
      end if
      if (iwrite .eq. 0) return
      idoub = IDSAFE(doub)
      write(nuni, 100) ident
c
      if (ident .gt. 10) go to 1000
      if      (ident .eq. 1) then
        write (nuni, 1) intgr
      else if (ident .eq. 2) then
        write (nuni, 2)
      else if (ident .eq. 3) then
        write (nuni, 3) intgr
      else if (ident .eq. 4) then
        if (intgr .eq. 1) then
          exword = 'above  '
        else
          exword = 'below  '
        end if
        write (nuni, 4) char, doub, exword
      else if (ident .eq. 5) then
        write (nuni, 5) char, doub, char
      else if (ident .eq. 6) then
        write (nuni, 6) 
      else if (ident .eq. 7) then
        write (nuni, 7) intgr, doub
      else if (ident .eq. 8) then
        write (nuni, 8) char
      else if (ident .eq. 9) then
        write (nuni, 9) intgr
      else if (ident .eq. 10) then
        write (nuni, 10) char
      end if
c
1     format(
     &' The solution dimension entered in the input data file',/,
     &' was reduced to ',i5,', which is the maximum allowed by',/,
     &' the parameter mdim.')
*
2     format(
     &' An integer was not found after the keyword POINTS.')
*
3     format(
     &' The dimension of the solution vector is less than the',/,
     &' the number of positive data (npdat = ',i4,');   one can',/,
     &' expect the  minimum possible recovered error found by',/,
     &' MICRON will  be larger than the true minimum, and the',/,
     &' solution may be unnecessarrily rough.')
*
4     format(
     &' The ',a5,' inversion diameter determined by examining',/,
     &' the data in micron.kri  and padding the intervals is',/,
     &' equal  to ',d11.3,'.   This  is ',a5,' the  maximum',/,
     &' inversion interval.')
*
5     format(
     &' The ',a5,'  inversion diameter  is  inside  the  inversion',/,
     &' interval  found by examining  the  support data contained',/,
     &' contained in micron.kri (',d10.3,').',/,
     &' For this inversion,  the true  distribution must  be zero',/,
     &' outside the specified ',a5,' interval.')
*
6     format(
     &' A valid key was not found following the keyword = message.')
*
7     format(
     &' There are ',i4,' unsupported gaps in the inversion interval.',/,
     &' "tolgap", the tolerance used in finding these gaps, equals',/,
     &  d10.3)
*
8     format(
     &' The keyword (',a16,') was not recognized.  The input',/,
     &' following the keyword has been ignored.')
*
9     format(
     &' The search  for the target smoothing parameter was',/,
     &' halted because an insufficient number of smoothing',/,
     &' iterations was specified(',i3,').')
*
10    format(
     &' The error model (',a16,         ') requested in the',/,
     &' input  data   file is  not recognized.   A  diagonal',/,
     &' error model will be assumed.')
c
1000  if (ident .gt. 20) go to 2000
c
      if      (ident .eq. 11) then
        write (nuni, 11) doub, char, intgr
      else if (ident .eq. 12) then
        write (nuni, 12) 
      else if (ident .eq. 13) then
        write (nuni, 13) 
      else if (ident .eq. 14) then
        write (nuni, 14) intgr, char
      else if (ident .eq. 15) then
        write (nuni, 15) intgr, char
      else if (ident .eq. 16) then
        write (nuni, 16) char, intgr
      else if (ident .eq. 17) then
        write (nuni, 17)
      else if (ident .eq. 18) then
        write (nuni, 18) doub, char
      else if (ident .eq. 19) then
        write (nuni, 19) intgr, char, idoub
      else if (ident .eq. 20) then
        write (nuni, 20) intgr
      end if
c
11    format(
     &' The sparse error element = ',d10.3,' will be ignored becuase',/,
     &' the ',a6,' index ( = ',i6,') is not valid.')
*
12    format(
     &' The specified inversion interval is too narrow;  MICRON will',/,
     &' use the calculated interval.')
*
13    format(
     &' Unable to read the old micron.kri data file.')
*
14    format(
     &' The kernel function response for channel ',i3,' of the ',a8,/,
     &' is nonzero below the minimum possible inversion diameter.')
*
15    format(' The kernel function response for channel ',i3,' of the ',
     &a8,/,' is nonzero above the maximum possible inversion diameter.')
*
16    format(
     &' The number of smoothing iterations entered in ',a12,' was',/,
     &' reduced to ',i4,' which is the maximum allowed by the MICRON',/,
     &' parameter miter.')
*
17    format(
     &' The number of smoothing iterations has been set to zero.')
*
18    format(
     &' The current smoothing parameter (',d10.3,') seems too ',a5)
*
19    format(' The integration routine DINTG has indicated difficulty wh
     &ile integrating',/,' channel ', i4, ' of the ', a8,/, 
     &'  The error index was set equal to ',i3)
*
20    format(
     &' The error array had ',i4,' eigenvalues equal to zero.')
c
c
2000  if (ident .gt. 30) go to 3000
c
      if      (ident .eq. 21) then
        write (nuni, 21) intgr
      else if (ident .eq. 22) then
        write (nuni, 22) doub
      else if (ident .eq. 23) then
        write (nuni, 23) intgr
      else if (ident .eq. 24) then
        write (nuni, 24) 
      else if (ident .eq. 25) then
        write (nuni, 25) intgr, char
      else if (ident .eq. 26) then
        write (nuni, 26) char
      else if (ident .eq. 27) then
        write (nuni, 27) intgr
      else if (ident .eq. 28) then
        write (nuni, 28) 
      else if (ident .eq. 29) then
        write (nuni, 29) 
      else if (ident .eq. 30) then
        write (nuni, 30) intgr
      end if
c
21    format(
     &' The error  array has ',i4,' negligible  rows.  This is',/,
     &' unrealistic because it implies some to the data have',/,
     &' no error.')
*
22    format(
     &' The QP solver reported a weak minimum for a smoothing',/,
     &' parameter of ',d11.3)
*
23    format(
     &' The column ',i3,' of the error array has been deleted',/,
     &' because all of the elements are negligible.')
*
24    format(' A numerical inconsistency occurred.',/,
     &' The minimum error was reset from a negative number to zero.')
*
25    format(
     &' The rank of the Hessian (',i3,') determined by ',a6, ' is',/,
     &' less than the dimension of the solution vector.')
*
26    format(
     &' No instrument names were found in ',a12,'.')
*
27    format(
     &' The minmization routine returned INFORM = ',i2)
*
28    format(
     &' The solution seems to be collapsing. Re-invert with',/,
     &' the sum of the errors constrained to zero.')
*
29    format(' ')
*
30    format(
     &' The minimum of the smoothing choice functional was not found',/,
     &' to within the specified tolerance after ',i3, ' iterations.')
c
c
3000  if (ident .gt. 40) go to 4000
c
      if      (ident .eq. 31) then
        write (nuni, 31) intgr, char
      else if (ident .eq. 32) then
        if (intgr .eq. 1) then
          exword = 'lower'
        else
          exword = 'upper'
        end if
        write (nuni, 32) char, exword
      else if (ident .eq. 33) then
        write (nuni, 33)
      else if (ident .eq. 34) then
        write (nuni, 34) intgr
      else if (ident .eq. 35) then
        write (nuni, 35) doub
      else if (ident .eq. 36) then
        write (nuni, 36) char, doub
      end if
c
31    format(
     &' The minimization routine BRENT returned INFORM = ',i3,/,
     &' while minimizing ',a4)
*
32    format(
     &' The minimium of the smoothing choice functional ',a4,/,
     &' is at the ',a5,' endpoint of the predetermined interval.')
*
33    format(
     &' The minimum error in  the  trasformed data is unreasonably',/,
     &' large.  Check the least squares solution (smp = 0) to make',/,
     &' sure the specified errors are realistic.')
*
34    format(
     &' The qradratic programming routine has indicated difficulty,',/,
     &' probably due to cycling. (INFORM =',i1,')')
*
35    format(
     &' The user specified smoothing parameter will be used in place',/,
     &' of the smoothing parameter determined by MICRON ( = ',d11.3,')')
*
36    format(
     &' The ',a5,'  inversion diameter  is  substantially outside',/,  
     &' the inversion interval  found by examining  the  support',/,
     &' data contained contained in micron.k (',d10.3,').')
c
c
4000  return
100   format(//,' Warning number ',i3)
*
*___end of WARNG___*
*
      end
