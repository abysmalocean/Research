      subroutine ERROR (ident,  nuni,   msglev, intgr,  doub,   char, 
     &                  fail)
************************************************************************
*
* ERROR is used to write messages and index the fail vector.
*
* INPUTS:
* ident  - the identification number of the message
* nuni   - the unit number the message is written to
* msglev - an integer which helps MICRON decide whether the message
*          should be written
* intgr  - an optional integer which may be written as part of the
*          message
* doub   - an optional double precision variable which may be written
*          as part of the message
* char   - an optional character string which may be written as part
*          of the message
*
* LAST MODIFIED: 5 JAN 1988
************************************************************************
      integer           fail,   ident,  IDSAFE, intgr,  idoub,  iwrite, 
     &                  lenc,   merrm,  msgchk, msglev, nuni
      double precision  doub
      parameter        (merrm = 100)
      character*(*)     char
      dimension         msgchk(merrm)
      data msgchk      /20,  10,  10,  10,  10,  10,  20,  20,  10,  10,
     1                  10,  10,  10,  20,  10,  20,  10,  10,  20,  10,
     2                  20,   0,  20,  20,  20,  20,  10,  10,  10,  10,
     3                  20,   0,   0,  20,  10,   0,  10,  10,   0,   0,
     4                  10,  10,   0,  10,   0,  10,   0,   0,   0,   0,
     5                  10,  10,  10,   0,   0,   0,  10,  10,   0,  10,
     6                   0,  10,  10,  10,  10,  10,   0,   0,   0,   0,
     7                   0,  10,   0,  10,   0,   0,   0,   0,  10,   0,
     8                   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     9                  20,   0,  20,   0,  10,  10,  10,   0,   0,   0/  
c
      call UINDEX(fail)
      call WRCHK(ident, msglev, nuni, msgchk, merrm, iwrite)
      if (iwrite .le. 0) return
      write(nuni, 100) ident
      idoub = IDSAFE(doub)
      lenc  = len(char)
c
      if (ident .gt. 10) go to 1000
c
      if      (ident .eq. 1) then
        write (nuni, 101) char
      else if (ident .eq. 2) then
        write (nuni, 102) char
      else if (ident .eq. 3) then
        write (nuni, 103) char
      else if (ident .eq. 4) then
        write (nuni, 104) char, intgr
      else if (ident .eq. 5) then
        write (nuni, 105) char, intgr
      else if (ident .eq. 6) then
        write (nuni, 106) char
      else if (ident .eq. 7) then
        write (nuni, 107) char
      else if (ident .eq. 8) then
        write (nuni, 108) char, idoub
        write (nuni, 200) intgr
      else if (ident .eq. 9) then
        write (nuni, 109) char, intgr, idoub
      else if (ident .eq. 10) then
        write (nuni, 110) char, intgr
      end if
      go to 9999
c
c
101   format(
     &' The data file ',a12,' generated too many warnings.')
*
102   format(
     &' An error occured while reading instrument information',/,
     &' in the data simulation file near the line:',/,a80)
*
103   format(
     &' The keyword ',a16,' in the micron.jbs data file was not'
     &' recognized.')
*
104   format(
     &' The number of data that need to be simulated in the file',/,
     & a12, ' cannot be larger than MICRON''s parameter mdat = ',i3,'.')
*
105   format(
     &' The number of instruments specified in the file ',a12,/,
     &' cannot be larger than MICRON''s parameter minsiv = ',i3)
*
106   format(
     &' An error occurred while reading test function information',/,
     &' near the line:',/,a80)
*
107   format(
     &' No positive data were found in the input file ',a12,'.')
*
108   format(
     &' Had trouble opening the output file ', a12, ' listed in',/,
     &' position ', i4, ' in micron.jbs.')
*
109   format(
     &' The number of instruments listed in ', a12, ' cannot be less',/,
     &' than 1 or greater than MICRON''s parameter minsiv = ',i3,'.')
*
110   format(
     &' The number of test functions specified in the file ',a12,/,
     &' must be less than MICRON''s parameter msum = ',i3,'.')
c
c
1000  if (ident .gt. 20) go to 2000
c
      if      (ident .eq. 11) then
        write (nuni, 111) intgr
      else if (ident .eq. 12) then
        write (nuni, 112) intgr, char
      else if (ident .eq. 13) then
        write (nuni, 113) char, intgr
      else if (ident .eq. 14) then
        write (nuni, 114) char
      else if (ident .eq. 15) then
        write (nuni, 115) char
      else if (ident .eq. 16) then
        write (nuni, 116) char
      else if (ident .eq. 17) then
        write (nuni, 117) char
      else if (ident .eq. 19) then
        write (nuni, 119) intgr, char, idoub
      else if (ident .eq. 20) then
        write (nuni, 120) intgr, char, doub
      end if
      go to 9999
c
c
111   format(
     &' The test function index specified in micron.sim cannot be',/,
     &' larger than MICRON''s parameter mfun = ',i3,', or less than 1.')
*
112   format(
     &1x,i3,' numbers were found next to the parameter name ',a12,/,
     &' in micron.sim instead of one.')
*
113   format(
     &' The number of data given in the data file ',    a12, ' is',/,
     &' greater than MICRON''s parameter mdat = ',i3,'.')
*
114   format(
     &' An error occured while trying to read instrument information',/,
     &' in the input file near the line:',/,a80)
*
115   format(
     &' The keyword ',a16,' in micron.sim is not valid.')
*
116   format(
     &' The required number of data were not provided for the',a8,'.')
*
117   format( 
     &' The instrument specified by',/,a78,/,
     &' in micron.jbs is not recognized.')
*
119   format(
     &' Unable to  read the  data  for channel',i3,' of the',/,
     &' instrument (',a8,') listed in position', i5,' of',/,
     &' the input data file.')
*
120   format(
     &' The instrument response for channel ',i3,' of the',/,
     &'(',a8,') was negative at diameter = ',d10.3)
c
c
2000  if (ident .gt. 30) go to 3000
c
      if      (ident .eq. 21) then
        write (nuni, 121) char
      else if (ident .eq. 23) then
        write (nuni, 123) intgr
      else if (ident .eq. 24) then
        write (nuni, 124) doub
        write (nuni, 300)
      else if (ident .eq. 25) then
        write (nuni, 125)
        write (nuni, 300)
      else if (ident .eq. 26) then
        write (nuni, 126) intgr
        write (nuni, 300)
      else if (ident .eq. 27) then
        write (nuni, 127) intgr
      else if (ident .eq. 28) then
        write (nuni, 128) char, idoub
        write (nuni, 200) intgr
      else if (ident .eq. 29) then
        write (nuni, 129) intgr
      end if
      go to 9999
c
c
121   format(
     &' A row number was not found for the dense error model near',/,
     &' the line',/,a80)
*
122   format(' ')
*
123   format(
     &' The row number corresponding to dense error model cannot be',/,
     &' larger than MICRON''s parameter mdat = ',i3,'.')
*
124   format(
     &' Unable to meet the target discrepancy within the specified',/,
     &' tolerance = ',d10.3,/,
     &' The sensitivity of the solution to the smoothing parameter',/,
     &' is greater than that allowed by the tolerance tol(3).')
*
125   format(
     &' Unable to meet the target discrepancy due to an insensitivity'
     &,/,' of the solution to the smoothing parameter.')
*
126   format(
     &' Unable to meet the target discrepancy after ',i4,' iter-',/,
     &' ations  which is the maximum allowed by the parameter',/,
     &' miter.')
*
127   format( 
     &1x,i5,' inversions have failed in succession.  This is',/,
     &' the maximum allowed as specified in micron.dat')
*
128   format(
     &' Unable to open the input file ', a12, ' listed in',/,
     &' position ',i5, ' in micron.jbs .')
*
129   format(
     &' The number of elements in a row of the error array cannot',/,
     &' be greater than MICRON''s parameter merr = ',i3,'.')
*
130   format(' ')
*
1130  format(
     &' Check the data file ', a12, ' for an echo of the data',/,
     &' that was read.')
c
3000  if (ident .gt. 40) go to 4000
c
c
      if (ident .eq. 31) then
        write (nuni, 131)
        write (nuni, 300)
      else if (ident .eq. 32) then
        write (nuni, 132) intgr
      else if (ident .eq. 34) then
        write (nuni, 134)
      else if (ident .eq. 35) then
        write (nuni, 135) doub
      else if (ident .eq. 37) then
        write (nuni, 137) intgr
      else if (ident .eq. 38) then
        write (nuni, 138) doub
      else if (ident .eq. 39) then
        write (nuni, 139) 
      else if (ident .eq. 40) then
        write (nuni, 140) 
      end if
      go to 9999
c
c
131   format(
     &' A possible  precision problem was detected.  The',/,
     &' discrepancy increased as the smoothing parameter',/,
     &' decreased.  The  smoother  and more accurate  of',/,
     &' the two solutions was reported.')
*
132   format(' ')
*
133   format(' ')
*
134   format(
     &' The magnitude of the largest eigenvalue of the error',/,
     &' array is nearly zero.')
*
135   format(
     &' The minimum positive number specified (',d10.3,') needs',/,
     &' to be larger than the machine precision.')
*
136   format(' ')
*
137   format(
     &' The singular value decomposition routine exited with an',/,
     &' error.  INFORM was set equal to ',i4)
*
138   format(
     &' The minimum diameter specified in micron.dat (',d10.3,')',/,
     &' should be larger than the smallest positive machine',/,
     &' number.')
*
139   format(' ')
*
140   format(' ')
c
c
4000  if (ident .gt. 50) go to 5000
c
      if      (ident .eq. 41) then
        write(nuni, 141) intgr
      else if (ident .eq. 42) then
        write (nuni, 142) char
        write (nuni, 200) intgr
      else if (ident .eq. 44) then
        write (nuni, 144) intgr, char, doub
      else if (ident .eq. 46) then
        write (nuni, 146)
      else if (ident .eq. 47) then
        write (nuni, 147) 
      else if (ident .eq. 48) then
        write (nuni, 148) 
      else if (ident .eq. 49) then
        write (nuni, 149)
      else if (ident .eq. 50) then
        write (nuni, 150) 
      end if
      go to 9999
c
141   format(
     &' The instrument number listed in micron.dat must be less than',/,
     &' MICRON''s parameter minst = ',i3,'.')
*
142   format(
     &' An error occurred while trying to ',a4,' the micron.kri file.')
*
144   format(
     &' The maximum value of the kernel function for channel ',i5,/,
     &' of the ', a8, ' is ', d10.3,'.  This should be much larger',/,
     &' than the minimum positive number wmach(5).')
*
146   format(' An inconsistency exists in the micron.kri data file.')
*
147   format(' ')
*
148   format(' ')
*
149   format(' ')
*
150   format(' ')
c
c
5000  if (ident .gt. 60) go to 6000
c
      if      (ident .eq. 51) then
        write (nuni, 151) char, intgr
      else if (ident .eq. 52) then
        write (nuni, 152) intgr
      else if (ident .eq. 53) then
        write (nuni, 153) intgr, idoub
      else if (ident .eq. 54) then
        write (nuni, 154) 
      else if (ident .eq. 55) then
        write (nuni, 155) 
      else if (ident .eq. 56) then
        write (nuni, 156) 
      else if (ident .eq. 57) then
        write (nuni, 157) doub
      else if (ident .eq. 58) then
        write (nuni, 158) char
        write (nuni, 200) intgr
      else if (ident .eq. 60) then
        write (nuni, 160)
      end if
      go to 9999
c
151   format(
     &' The number of channels for the ',a8,' listed in micron.dat',/,
     &' must be between 1 and MICRON''s parameter mdat = ',i3,'.')
*
152   format(
     &' The current record length of micron.kri needs to be',/,
     &' increased to at least ',i4)
*
153   format(
     &' No number was found for parameter ',i3,' of the ',a8,/,
     &' in micron.dat.')
*
154   format(' ')
*
155   format(' ')
*
156   format(' ')
*
157   format(
     &' The largest positive number, wmach(7) = ', d10.3,/,
     &' entered in micron.dat is ridiculously small.')
*
158   format(
     &' An error occurred while trying to open the data file ',a12,'.') 
*
160   format(
     &' There were no test function names found in micron.sim.')
c
c
6000  if (ident .gt. 70) go to 7000
c
      if (ident .eq. 61) then
        write (nuni, 161) 
      else if (ident .eq. 62) then
        write (nuni, 162) intgr, char
      else if (ident .eq. 63) then
        write (nuni, 163) char
      else if (ident .eq. 64) then
        write (nuni, 164)
      else if (ident .eq. 65) then
        write (nuni, 165) doub
      else if (ident .eq. 66) then
        write (nuni, 166) doub
      else if (ident .eq. 67) then
        write (nuni, 167) 
      else if (ident .eq. 68) then
        write (nuni, 168) 
      else if (ident .eq. 69) then
        write (nuni, 169) 
      end if
      go to 9999
c
161   format(' ')
*
162   format(
     &' MICRON was only able to find ',i2,' numeric entries after',/,
     &' the keyword ',       a16,' listed in micron.dat.')
*
163   format(
     &' The keyword ',a16,' is not recognized as a valid keyword',/,
     &' for micron.dat.')
*
164   format(
     &' The singular  value decomposition  routine  reported',/,
     &' computational difficulties.  Reinvert the data using',/,
     &' the  method  of  discrepancy and  specify  the first',/,
     &' smoothing parameter.')
*
165   format(
     &' The maximum diameter, dmax = ', d10.3, ' entered in ',/,
     &' micron.dat is ridiculously large.')
*
166   format(
     &' The minimum positive number, wmach(5) = ', d10.3, ' entered',/,
     &' in micron.dat needs to be positive.')
*
167   format(' ')
*
168   format(' ')
*
169   format(' ')
c
c
7000  if (ident .gt. 80) go to 8000
c
      if      (ident .eq. 71) then
        write (nuni, 171) 
      else if (ident .eq. 72) then
        write (nuni, 172) 
        write (nuni, 200) intgr
      else if (ident .eq. 73) then
        write (nuni, 173)
      else if (ident .eq. 74) then
        write (nuni, 174) char
      else if (ident .eq. 75) then
        write (nuni, 175)
      else if (ident .eq. 76) then
        write (nuni, 176) 
      else if (ident .eq. 77) then
        write (nuni, 177) 
      else if (ident .eq. 78) then
        write (nuni, 178) 
      else if (ident .eq. 79) then
        write (nuni, 179) char
      else if (ident .eq. 80) then
        write (nuni, 180) 
      end if
      go to 9999
c
c
171   format(' ')
*
172   format(
     &' An error occurred while trying to open the data file',/,
     &' micron.sim.')  
*
173   format(' ')
*
174   format(
     &' The smoothing optimum, ',a16, ', in micron.dat is not',/,
     &' recogmized.')
*
175   format(' ')
*
176   format(' ')
*
177   format(' ')
*
178   format(' ')
*
179   format(
     &' There were no test functions found in the simulation',/,
     &' data file named ',a12)
*
180   format(' ')
c
c
8000  if (ident .gt. 90) go to 9000
c
      if      (ident .eq. 81) then
        write (nuni, 181) 
      else if (ident .eq. 82) then
        write (nuni, 182) 
        write (nuni, 200) idoub
      else if (ident .eq. 83) then
        write (nuni, 183) 
      else if (ident .eq. 84) then
        write (nuni, 184) 
      else if (ident .eq. 85) then
        write (nuni, 185) 
      else if (ident .eq. 86) then
        write (nuni, 186) 
      else if (ident .eq. 87) then
        write (nuni, 187) 
      else if (ident .eq. 88) then
        write (nuni, 188) 
      else if (ident .eq. 89) then
        write (nuni, 189) 
      else if (ident .eq. 90) then
        write (nuni, 190) 
      end if
      go to 9999
c
181   format(' ')
*
182   format(' ')
*
183   format(' ')
*
184   format(' ')
*
185   format(' ')
*
186   format(' ')
*
187   format(' ')
*
188   format(' ')
*
189   format(' ')
*
190   format(' ')
c
c
9000  if (ident .gt. 100) go to 10000
c
      if      (ident .eq. 91) then
        write (nuni, 191)
      else if (ident .eq. 92) then
        write (nuni, 192)
      else if (ident .eq. 93) then
        write (nuni, 193) intgr
      else if (ident .eq. 94) then
        write (nuni, 194) 
      else if (ident .eq. 95) then
        write (nuni, 195) intgr
      else if (ident .eq. 96) then
        write (nuni, 196)
      else if (ident .eq. 97) then
        write (nuni, 197) char
        write (nuni, 200) intgr
      end if
c
191   format(
     &' An error occurred while trying to find the solution dimension.')
*
192   format(' ')
*
193   format(
     &' The  maximum number of general constraints determined by',/,
     &' MICRON''s FORTRAN parameter [mclin] needs to be increased',/,
     &' to ',i4,' for this inversion.')
*
194   format(' ')
*
195   format(
     &' The minimization routine has returned an error condition.',/,
     &' INFORM = ',i4)
*
196   format(
     &' Numerical difficulty was encountered while decomposing the',/,
     &' Hessian.')
*
197   format(
     &' A FORTRAN error occurred while trying to ',a4,' the drop',/,
     &' file micron.r.')
c
c
10000 return
9999  return
c
100   format(//,' Error number ',i3)
200   format(' The iostat was set equal to ', i5,'.')
300   format(' Check the calculation log for more information.')
      end
