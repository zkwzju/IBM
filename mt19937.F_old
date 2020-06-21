*   genrand() generates one pseudorandom real number (double)
* which is uniformly distributed on [0,1]-interval, for each
* call. sgenrand(seed) set initial values to the working area
* of 624 words. Before genrand(), sgenrand(seed) must be
* called once. (seed is any 32-bit integer except for 0).
* Integer generator is obtained by modifying two lines.
*   Coded by Takuji Nishimura, considering the suggestions by
* Topher Cooper and Marc Rieffel in July-Aug. 1997.
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Library General Public
* License as published by the Free Software Foundation; either
* version 2 of the License, or (at your option) any later
* version.
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU Library General Public License for more details.
* You should have received a copy of the GNU Library General
* Public License along with this library; if not, write to the
* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
* 02111-1307  USA
*
* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
* When you use this, send an email to: matumoto@math.keio.ac.jp
* with an appropriate reference to your work.
*
************************************************************************
* Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
*
*   genrand()      -&gt; double precision function grnd()
*   sgenrand(seed) -&gt; subroutine sgrnd(seed)
*                     integer seed
*
* This program uses the following non-standard intrinsics.
*   ishft(i,n): If n&gt;0, shifts bits in i by n positions to left.
*               If n&lt;0, shifts bits in i by n positions to right.
*   iand (i,j): Performs logical AND on corresponding bits of i and j.
*   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
*   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
*
************************************************************************
* Modifications, in particular the separating of the integer generation
* part of grnd() into a separate function igrnd(), a crude 
* simplification of the int to float conversion in grnd(), wrapping with
* a module and appending some RANDLIB code to generate non-uniform 
* distributions by Robert I A Patterson 2003-4.
*
*************************************************************************



* this main() outputs first 1000 generated numbers
c      program main
*
c      implicit integer(i-n)
c      implicit double precision(a-h,o-z)
*
c      parameter(no=1000)
c      dimension r(0:7)
*
*      call sgrnd(4357)
*                         any nonzero integer can be used as a seed
c      do 1000 j=0,no-1
c        r(mod(j,8))=grnd()
c        if(mod(j,8).eq.7) then
c          write(*,'(8(f8.6,'' ''))') (r(k),k=0,7)
c        else if(j.eq.no-1) then
c          write(*,'(8(f8.6,'' ''))') (r(k),k=0,mod(no-1,8))
c        endif
c 1000 continue
*
c      stop
c      end
************************************************************************
      MODULE MersenneTwister

	CONTAINS

      subroutine sgrnd(seed)
*
      implicit integer(a-z)
*
* Period parameters
      parameter(N     =  624)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
*
*      setting initial seeds to mt[N] using
*      the generator Line 25 of Table 1 in
*      [KNUTH 1981, The Art of Computer Programming
*         Vol. 2 (2nd Ed.), pp102]
*
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue
*
      return
      end subroutine sgrnd
************************************************************************
!r      double precision function grnd()
      integer function igrnd()
*
      implicit integer(a-z)
*
* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
*                                    constant vector a
      parameter(UMASK = -2147483648)
*                                    most significant w-r bits
      parameter(LMASK =  2147483647)
*                                    least significant r bits
* Tempering parameters
      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)
*
      dimension mt(0:N-1)
*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
*                     mti==N+1 means mt[N] is not initialized
*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
*                        mag01(x) = x * MATA for x=0,1
*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
*
      if(mti.ge.N) then
*                       generate N words at one time
        if(mti.eq.N+1) then
*                            if sgrnd() has not been called,
          call sgrnd(4357)
*                              a default initial seed is used
        endif
*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
*
!r      if(y.lt.0) then
!r        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!r      else
!r        grnd=dble(y)/(2.0d0**32-1.0d0)
!r      endif
*
      igrnd = y
      return
	end function igrnd
!r      end function grnd

      double precision function grnd()
	  implicit none
	  integer y

!r    get the random number and drop the sign bit, only works for 32 bit signed integers
        y=IBCLR(igrnd(), BIT_SIZE(y) - 1) 
!r    We should now have a random integer in [0,(2**BIT_SIZE(1)-1)-1]
      grnd = (dble(y)+0.5) / 2147483648.0d0

	end function grnd


C**********************************************************************C
C                                                                      C
C     Following code taken from RANDLIB by Brown and Lovato            C
C     http://odin.mdacc.tmc.edu/anonftp/#RANDLIB                       C
C                                                                      C
C     ?pub domain




      REAL FUNCTION snorm()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
C
C     .. Local Scalars ..
      REAL aa,s,tt,u,ustar,w,y
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL a(32),d(31),h(31),t(31)
C     ..
C     .. External Functions ..
!r      REAL ranf
!r      EXTERNAL ranf
C     ..
C     .. Intrinsic Functions ..
!r      INTRINSIC float,int
C     ..
C     .. Save statement ..
C     JJV added a Save statement for arrays initialized in Data statmts
      SAVE a,d,t,h
C     ..
C     .. Data statements ..
      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991,
     +     .2372021,.2776904,.3186394,.3601299,.4022501,.4450965,
     +     .4887764,.5334097,.5791322,.6260990,.6744898,.7245144,
     +     .7764218,.8305109,.8871466,.9467818,1.009990,1.077516,
     +     1.150349,1.229859,1.318011,1.417797,1.534121,1.675940,
     +     1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     +     .1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     +     .1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     +     .1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     +     .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     +     .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     +     .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     +     .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     +     .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     +     .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     +     .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     +     .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     +     .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     +     .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     +     .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     +     .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     +     .7010474/
C     ..
C     .. Executable Statements ..
C
   10 u = grnd()
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
   20 u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
C
C                                START CENTER
C
   30 ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
C
C                                EXIT   (BOTH CASES)
C
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
C
C                                CENTER CONTINUED
C
   60 u = grnd()
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      ustar = grnd()
   80 IF (ustar.GT.tt) GO TO 50
   90 u = grnd()
      IF (ustar.GE.u) GO TO 70
      ustar = grnd()
      GO TO 40
C
C                                START TAIL
C
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
  130 u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 ustar = grnd()
      IF (ustar.GT.tt) GO TO 50
  170 u = grnd()
      IF (ustar.GE.u) GO TO 150
      u = grnd()
      GO TO 140

      END FUNCTION snorm

      REAL FUNCTION sexpo()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM THE                 C
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  C
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       C
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
C
C     JJV added a Save statement for q (in Data statement)
C     .. Local Scalars ..
      REAL a,q1,u,umin,ustar
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL q(8)
C     ..
C     .. External Functions ..
!r      REAL ranf
!r      EXTERNAL ranf
C     ..
C     .. Equivalences ..
      EQUIVALENCE (q(1),q1)
C     ..
C     .. Save statement ..
      SAVE q
C     ..
C     .. Data statements ..
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
C     ..
C
  210 a = 0.0
      u = grnd()
      GO TO 230

  220 a = a + q1
  230 u = u + u
C     JJV changed the following to reflect the true algorithm and
C     JJV prevent unpredictable behavior if U is initially 0.5.
C      IF (u.LE.1.0) GO TO 220
      IF (u.LT.1.0) GO TO 220
  240 u = u - 1.0
      IF (u.GT.q1) GO TO 260
  250 sexpo = a + u
      RETURN

  260 i = 1
      ustar = grnd()
      umin = ustar
  270 ustar = grnd()
      IF (ustar.LT.umin) umin = ustar
  280 i = i + 1
      IF (u.GT.q(i)) GO TO 270
  290 sexpo = a + umin*q1
      RETURN

      END FUNCTION sexpo

      INTEGER FUNCTION ignpoi(mu)
C**********************************************************************
C
C     INTEGER FUNCTION IGNPOI( MU )
C
C                    GENerate POIsson random deviate
C
C
C                              Function
C
C
C     Generates a single random deviate from a Poisson
C     distribution with mean MU.
C
C
C                              Arguments
C
C
C     MU --> The mean of the Poisson distribution from which
C            a random deviate is to be generated.
C                              REAL MU
C     JJV                    (MU >= 0.0)
C
C     IGNPOI <-- The random deviate.
C                              INTEGER IGNPOI (non-negative)
C
C
C                              Method
C
C
C     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C
C               Ahrens, J.H. and Dieter, U.
C               Computer Generation of Poisson Deviates
C               From Modified Normal Distributions.
C               ACM Trans. Math. Software, 8, 2
C               (June 1982),163-179
C
C**********************************************************************
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C                                                                      C
C     P O I S S O N  DISTRIBUTION                                      C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER GENERATION OF POISSON DEVIATES                C
C               FROM MODIFIED NORMAL DISTRIBUTIONS.                    C
C               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. C
C                                                                      C
C     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  C
C                                                                      C
C**********************************************************************C
C
C      INTEGER FUNCTION IGNPOI(IR,MU)
C
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
C             MU=MEAN MU OF THE POISSON DISTRIBUTION
C     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
C
C
C
C     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR CASE B
C     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
C     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
C
C
C
C     SEPARATION OF CASES A AND B
C
C     .. Scalar Arguments ..
      REAL mu
C     ..
C     .. Local Scalars ..
      REAL a0,a1,a2,a3,a4,a5,a6,a7,b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,
     +     fk,fx,fy,g,muold,muprev,omega,p,p0,px,py,q,s,t,u,v,x,xx
C     JJV I added a variable 'll' here - it is the 'l' for CASE A
      INTEGER j,k,kflag,l,ll,m
C     ..
C     .. Local Arrays ..
      REAL fact(10),pp(35)
C     ..
C     .. External Functions ..
!r      REAL ranf,sexpo,snorm
!r      EXTERNAL ranf,sexpo,snorm
C     ..
C     .. Intrinsic Functions ..
!r      INTRINSIC abs,alog,exp,float,ifix,max0,min0,sign,sqrt
C     ..
C     JJV added this for case: mu unchanged
C     .. Save statement ..
      SAVE s, d, l, ll, omega, c3, c2, c1, c0, c, m, p, q, p0,
     +     a0, a1, a2, a3, a4, a5, a6, a7, fact, muprev, muold
C     ..
C     JJV end addition - I am including vars in Data statements
C     .. Data statements ..
C     JJV changed initial values of MUPREV and MUOLD to -1.0E37
C     JJV if no one calls IGNPOI with MU = -1.0E37 the first time,
C     JJV the code shouldn't break
      DATA muprev,muold/-1.0E37,-1.0E37/
      DATA a0,a1,a2,a3,a4,a5,a6,a7/-.5,.3333333,-.2500068,.2000118,
     +     -.1661269,.1421878,-.1384794,.1250060/
      DATA fact/1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./
C     ..
C     .. Executable Statements ..

      IF (mu.EQ.muprev) GO TO 310
      IF (mu.LT.10.0) GO TO 420
C
C     C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
C
C     JJV This is the case where I changed 'l' to 'll'
C     JJV Here 'll' is set once and used in a comparison once

      muprev = mu
      s = sqrt(mu)
      d = 6.0*mu*mu
C
C             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
C             PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
C             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
C
      ll = ifix(mu-1.1484)
C
C     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
C
  310 g = mu + s*snorm()
      IF (g.LT.0.0) GO TO 320
      ignpoi = ifix(g)
C
C     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
C
      IF (ignpoi.GE.ll) RETURN
C
C     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
C
      fk = float(ignpoi)
      difmuk = mu - fk
      u = grnd()
      IF (d*u.GE.difmuk*difmuk*difmuk) RETURN
C
C     STEP P. PREPARATIONS FOR STEPS Q AND H.
C             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
C             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
C             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
C             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
C             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
C
  320 IF (mu.EQ.muold) GO TO 330
      muold = mu
      omega = .3989423/s
      b1 = .4166667E-1/mu
      b2 = .3*b1*b1
      c3 = .1428571*b1*b2
      c2 = b2 - 15.*c3
      c1 = b1 - 6.*b2 + 45.*c3
      c0 = 1. - b1 + 3.*b2 - 15.*c3
      c = .1069/mu
  330 IF (g.LT.0.0) GO TO 350
C
C             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
C
      kflag = 0
      GO TO 370
C
C     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
C
  340 IF (fy-u*fy.LE.py*exp(px-fx)) RETURN
C
C     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
C             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
C             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
C
  350 e = sexpo()
      u = grnd()
      u = u + u - 1.0
      t = 1.8 + sign(e,u)
      IF (t.LE. (-.6744)) GO TO 350
      ignpoi = ifix(mu+s*t)
      fk = float(ignpoi)
      difmuk = mu - fk
C
C             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
C
      kflag = 1
      GO TO 370
C
C     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
C
  360 IF (c*abs(u).GT.py*exp(px+e)-fy*exp(fx+e)) GO TO 350
      RETURN
C
C     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
C             CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
C
  370 IF (ignpoi.GE.10) GO TO 380
      px = -mu
      py = mu**ignpoi/fact(ignpoi+1)
      GO TO 410
C
C             CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
C             A0-A7 FOR ACCURACY WHEN ADVISABLE
C             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
C
  380 del = .8333333E-1/fk
      del = del - 4.8*del*del*del
      v = difmuk/fk
      IF (abs(v).LE.0.25) GO TO 390
      px = fk*alog(1.0+v) - difmuk - del
      GO TO 400

  390 px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) -
     +     del
  400 py = .3989423/sqrt(fk)
  410 x = (0.5-difmuk)/s
      xx = x*x
      fx = -0.5*xx
      fy = omega* (((c3*xx+c2)*xx+c1)*xx+c0)
      IF (kflag) 340,340,360
C
C     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
C
C     JJV changed MUPREV assignment from 0.0 to initial value
  420 muprev = -1.0E37
      IF (mu.EQ.muold) GO TO 430
C     JJV added argument checker here
      IF (mu.GE.0.0) GO TO 425
      WRITE (*,*) 'MU < 0 in IGNPOI - ABORT'
      WRITE (*,*) 'Value of MU: ',mu
      STOP 'MU < 0 in IGNPOI - ABORT'
C     JJV added line label here
  425 muold = mu
      m = max0(1,ifix(mu))
      l = 0
      p = exp(-mu)
      q = p
      p0 = p
C
C     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
C
  430 u = grnd()
      ignpoi = 0
      IF (u.LE.p0) RETURN
C
C     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
C             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
C             (0.458=PP(9) FOR MU=10)
C
      IF (l.EQ.0) GO TO 450
      j = 1
      IF (u.GT.0.458) j = min0(l,m)
      DO 440 k = j,l
          IF (u.LE.pp(k)) GO TO 480
  440 CONTINUE
      IF (l.EQ.35) GO TO 430
C
C     STEP C. CREATION OF NEW POISSON PROBABILITIES P
C             AND THEIR CUMULATIVES Q=PP(K)
C
  450 l = l + 1
      DO 460 k = l,35
          p = p*mu/float(k)
          q = q + p
          pp(k) = q
          IF (u.LE.q) GO TO 470
  460 CONTINUE
      l = 35
      GO TO 430

  470 l = k
  480 ignpoi = k
      RETURN

      END FUNCTION ignpoi

      END MODULE MersenneTwister
