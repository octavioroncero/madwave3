************************  dwigner  **************************************

      subroutine dwigner(dj,jmax2,m2,mp2,beta,ndim)
      implicit real*8(a-h,o-z)

*     ***********************************************
*     *   Wigner matrizes: d^J_{M,M'} (cos(beta))   *
*     *       kept in dj(2*J) to consider           *
*     *         halfinteger values of J             *
*     *  input:                                     *
*     *     jmax2 is 2 * jmax, jmax being           *
*     *                   the highest  j required   *
*     *     m2 and mp2: are (m * 2) and (mp *2)     *
*     *                   respectively              *
*     *     beta: angle defined between 0 and Pi    *
*     ***********************************************

      dimension dj(0:ndim)

      if(jmax2.gt.ndim)then
           write(6,*)' ** be carefull with dimensions in dwigner **'
           write(6,*)'    2 jmax > ndim',jmax2,ndim
           stop
      endif

      do j=0,jmax2
         dj(j)=0.d0
      enddo

      jmin2 = max0(iabs(m2),iabs(mp2)) 

      y=dcos(beta)

***> jmin=0 ---> M=M'=0 ---> d^J_{M,Mp} = P_J(cos(beta))
********************************************************

      if(jmin2.eq.0)then
         dj(jmin2)=1.d0
         dj(jmin2+2)=y
         do ind=jmin2+4,jmax2,2
            xj=dble(ind)*0.5d0
            indm1=ind-2
            indm2=ind-4
            dj(ind)=y*(2.d0*xj-1.d0)*dj(indm1)-(xj-1.d0)*dj(indm2)
            dj(ind)=dj(ind)/xj
         enddo
      elseif(jmin2.gt.0)then

***> d^J_{M,M'} for J=jmin=max0(iabs(M),iabs(M')) > 0
*a)      Normalization factor

         facnum=0.d0
         do i=1,jmin2
            facnum=facnum+dlog(dble(i))
         enddo

         if(jmin2.eq.iabs(m2))then
            mreal2=mp2
         else
            mreal2=m2
         endif

         facden1=0.d0
         if(jmin2+mreal2.gt.0)then
            do i=2,jmin2+mreal2,2
               facden1=facden1+dlog(dble(i)*0.5d0)
            enddo
         endif

         facden2=0.d0
         if(jmin2-mreal2.gt.0)then
            do i=2,jmin2-mreal2,2
               facden2=facden2+dlog(dble(i)*0.5d0)
            enddo
         endif

         fact=dexp(0.5d0*(facnum-facden1-facden2))

*b) 
         sign=1.d0
         xsin=dsin(beta*0.5d0)
         xcos=dcos(beta*0.5d0)

         if(jmin2.eq.abs(m2))then 
            if(jmin2.eq.m2)then
               if(mod((jmin2-mp2),2).eq.0)then
                  sign=(-1)**(int(dble(jmin2-mp2)*0.5d0+0.5d0))
                  a=xsin**(dble(jmin2-mp2)*0.5d0)
                  b=xcos**(dble(jmin2+mp2)*0.5d0)
               else
                  write(6,*)'  ** danger: complex phase **'
                  stop
               endif
            elseif(jmin2.eq.-m2)then
               a=xsin**(dble(jmin2+mp2)*0.5d0)
               b=xcos**(dble(jmin2-mp2)*0.5d0)
            endif
         elseif(jmin2.eq.abs(mp2))then
            if(jmin2.eq.mp2)then
               a=xsin**(dble(jmin2-m2)*0.5d0)
               b=xcos**(dble(jmin2+m2)*0.5d0)
            elseif(jmin2.eq.-mp2)then
               if(mod((jmin2+m2),2).eq.0)then
                  sign=(-1)**int(dble(jmin2+m2)*0.5d0+0.5d0)
                  a=xsin**(dble(jmin2+m2)*0.5d0)
                  b=xcos**(dble(jmin2-m2)*0.5d0)
               else
                  write(6,*)'  ** danger: complex phase **'
                  stop
               endif
            endif
         endif

*c)

         dj(jmin2)=fact*sign*a*b


***> d^J_{M,M'} for J=jmin+1

         xj=dble(jmin2)*0.5d0
         xjp=xj+1.d0
         xm=dble(m2)*0.5d0
         xmp=dble(mp2)*0.5d0

         facnum=xjp*(2.d0*xj+1.d0)
         facden1=xjp*xjp-xm*xm
         facden2=xjp*xjp-xmp*xmp
         factot=facnum/dsqrt(facden1*facden2)

         facj= y - xm * xmp /(xj*(xj+1.d0))
         
         ind=jmin2

         indp=ind+2
         if(indp.le.jmax2)then
         dj(indp)=factot*facj*dj(ind)

***> d^J_{M,M'} for J>jmin+1

         do indp=jmin2+4,jmax2,2
            ind=indp-2
            indm=indp-4

            xjp=dble(indp)*0.5d0
            xj=xjp-1.d0
            xjm=xjp-2.d0

            facnum=xjp*(2.d0*xj+1.d0)
            facden1=xjp*xjp-xm*xm
            facden2=xjp*xjp-xmp*xmp
            factot=facnum/dsqrt(facden1*facden2)

            facj= y - xm*xmp/(xj*(xj+1.d0))

            fac1jm=xj*xj-xm*xm
            fac2jm=xj*xj-xmp*xmp
            facden=xj*(2.d0*xj+1.d0)
            facjm=dsqrt(fac1jm*fac2jm)/facden

            dj(indp)=factot*( facj*dj(ind) -facjm*dj(indm) )
         enddo


      endif
****************************************

      endif

      return
      end
*************************** gauleg **************************

      SUBROUTINE GAULEG(W,XX,n)
      IMPLICIT real*8 (A-H,O-Z)

      dimension W(n),XX(n)

      x1=-1.d0
      x2=1.d0

      YN=DFLOAT(N)
      EPS=1.D-15
      M=(N+1)/2
      XM=0.5D0*(X1+X2)
      XL=0.5D0*(X2-X1)
      PI=DACOS(-1.D0)
      DO 12 I=1,M
      YI=DFLOAT(I)
      Z=DCOS(PI*(YI-0.25D0)/(YN+0.5D0))
 1    CONTINUE
      P1=1.D0
      P2=0.D0
      DO 11 J=1,N
      P3=P2
      P2=P1
      YJ=DFLOAT(J)
      P1=((2.D0*YJ-1.D0)*Z*P2-(YJ-1.D0)*P3)/YJ
11    CONTINUE
      PP=YN*(Z*P1-P2)/(Z*Z-1.D0)
      Z1=Z
      Z=Z1-P1/PP
      IF(DABS(Z-Z1).GT.EPS)GO TO 1
      XX(I)=XM-XL*Z
      XX(N+1-I)=XM+XL*Z
      W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
12    CONTINUE

      RETURN
      END

********************  FACTORIAL  *********************************

      subroutine factorial

      IMPLICIT real*8(A-H,O-Z)
      common/fct/fact(0:10000)

      FACT(0)=0.D0
      fact(1)=0.d0
      n=10000
      DO 10 I=2,N
         FACT(i)=FACT(i-1)+DLOG(DFLOAT(I))
   10 CONTINUE

      RETURN 
      END
                                                            

*************************  tqli  ****************************

      SUBROUTINE TQLI(D,E,N,NP)
      implicit real*8(a-h,o-z)
****************************************************
***    diagonalization of a tridiagonal matrix   ***
***         from numerical recipies              ***
****************************************************
      DIMENSION D(NP),E(NP)
      IF (N.GT.1) THEN
c overide of silly num recipes convention for off-diagonal
c elements.  It is now assumed off diags are in E(1)..E(N-1)
c        DO 11 I=2,N
c          E(I-1)=E(I)
c11      CONTINUE
        E(N)=0.d0
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=dABS(D(M))+dABS(D(M+1))
            IF (dABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.500)write(6,*) 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.d0*E(L))
            R=dSQRT(G**2+1.d0)
            G=D(M)-D(L)+E(L)/(G+dSIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(dABS(F).GE.dABS(G))THEN
                C=G/F
                R=dSQRT(C**2+1.d0)
                E(I+1)=F*R
                S=1.d0/R
                C=C*S
              ELSE
                S=F/G
                R=dSQRT(S**2+1.d0)
                E(I+1)=G*R
                C=1.d0/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.d0
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END

**************************  SPLIN  **********************************

c     FUNCTION SPLINQ(F,X,IOLD,NX,R,ndim)
      subroutine SPLINQQ(F,X,IOLD,NX,R,ndim,spl)
      IMPLICIT real*8 (A-H,O-Z)
C **********************************************************************
C  USE IOLD TO SAVE TIME IF SPLINT IS GOING TO BE CALLED FOR INCREASING
C  VALUES OF R:
C    - FOR 1ST CALL, SET IOLD TO 2
C    - FOR THE FOLLOWING CALLS, USE THE OUT VALUE OF IOLD
C  IF THE VALUES OF R ARE NOT MONOTONICALLY INCREASING, USE SPLINT OR
C  SPLINQQ OR SPLINP
C----------------------------------------------------------------------
C  C. LEFORESTIER, UNIVERSITE PARIS-SUD, ORSAY, FRANCE.
C  MODIFICATIONS BY N. HALBERSTADT, CNRS, ORSAY, FRANCE.
C **********************************************************************
      DIMENSION U(4)
      DIMENSION F(ndim,2),X(ndim)
      DATA UN/1D0/,TWO/2D0/,THREE/3D0/
C
      IF(R.GE.X(NX)) GO TO 30
      DO 10 IDOL = IOLD, NX
      IF(R.LT.X(IDOL)) GOTO 20
   10 CONTINUE
   20 IOLD = IDOL
      HI=X(IDOL)-X(IDOL-1)
      XR=(R-X(IDOL-1))/HI
      U(1)=XR*XR*(-TWO*XR+THREE)
      U(3)=HI*XR*XR*(XR-UN)
      U(2)=UN-U(1)
      U(4)=HI*XR*((XR-TWO)*XR+UN)
      SPLINQ=U(1)*F(IDOL,1)+U(2)*F(IDOL-1,1)
     &      +U(3)*F(IDOL,2)+U(4)*F(IDOL-1,2)
      spl=splinq
      RETURN
C  >>> WARNING: THIS EXTENSION IS NOT VALID FOR EVERY POTENTIAL <<<
30    RO=X(NX)
      YO=F(NX,1)
      N2X=2*NX
      YP=F(NX,2)
      AIN=YO+YP*RO/6.D0
      C8=-AIN*3.D0*RO**8
      C6=YO*RO**6-C8/RO/RO
      SPLINQ=C6/R**6+C8/R**8
      spl=splinq
         RETURN
      END
      SUBROUTINE SPLSET(F,X,NX,ndim)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C***********************************************************************
C*                                                                     *
C*        THIS ROUTINE SETS THE SPLINE INTERPOLATION ON THE GRID       *
C*    (X(I),I=1,NX) FOR THE FUNCTION (F(I),I=1,NX).                    *
C***********************************************************************
C
C
      PARAMETER (NXMAX=20000)
C
      DIMENSION F(ndim,2),X(ndim)
      DIMENSION HX(NXMAX)
      DIMENSION RLX(NXMAX-1),RMUX(NXMAX-1),XI(NXMAX-1),B(NXMAX-1)
      DIMENSION AB(4),YZ(4),A(4)
C
      DATA UN/1.0D0/,THREE/3.0D0/
C
      IF (NX.GT.NXMAX) GO TO 999
C
      NX2=NX-2
C     COMPUTE THE HX'S
      DO 10 I=2,NX
   10 HX(I-1)=X(I)-X(I-1)
C     COMPUTE LAMBDA'S & MU'S
      DO 40 I=1,NX2
      RLX(I)=HX(I+1)/(HX(I)+HX(I+1))
   40 RMUX(I)=UN-RLX(I)
      MAN=NX-3
      DO 60 I=1,4
      A(I)=X(MAN)
   60 MAN=MAN+1
C
C     SPLINE-FIT DE P(X)
      DO 110 I=1,4
      AB(I)=F(I,1)
  110 YZ(I)=F(NX+I-4,1)
      P0=DLAGRA(X,AB,4,1)
      F(1,2)=P0
      PN=DLAGRA(A,YZ,4,4)
      F(NX,2)=PN
C     CALCUL SECOND MEMBRE
      DO 120 I=1,NX2
  120 B(I)=THREE*RLX(I)/HX(I)*(F(I+1,1)-F(I,1))
     & +THREE*RMUX(I)/HX(I+1)*(F(I+2,1)-F(I+1,1))
      B(1)=B(1)-RLX(1)*P0
      B(NX2)=B(NX2)-RMUX(NX2)*PN
      CALL JORDAN(RMUX,RLX,XI,NX2,B)
      DO 100 I=1,NX2
  100 F(I+1,2)=XI(I)
      RETURN
C
 999  CONTINUE
      PRINT 9000
      PRINT 9999, NXMAX, NX
      STOP
C
 9000 FORMAT(//,2X,20('*'),' STOP IN SPLSET ',20('*'),/)
 9999 FORMAT(2X,'DIMENSION PARAMETER NXMAX = ',I5,' TOO SMALL, '
     &       ,I5,' REQUIRED')
      END
*********************************  SCHR *************************

      SUBROUTINE SCHR(E0,RMIN,RMAX,N,MAXIT,EPS,E2,KV,ITRY,v,p,npunt)
      IMPLICIT real*8(A-H,O-Z) 
      DIMENSION Y(3),p(n),v(n)

      ITRY=0
      H=(RMAX-RMIN)/DFLOAT(N-1)
      H2=H**2
      HV=H2/12.d0
      E=E0
      TEST=-1.d0
      DE=0.d0
      DO 1 I=1,N
1     P(I)=0.d0
C     BOUCLE DES ITERATIONS
   12 DO 171 IT=1,MAXIT
     
      XIT=IT
C     INTEGRATION VERS L-INTERIEUR.PREMIERS PAS
      P(N)=1.d-30
      GN=V(N)-E
      GI=V(N-1)-E
C     E EST-IL TROP GRAND
      IF(GI.GE.0.d0) GO TO 36
      E=V(N-2)
      GO TO 36
  900 PRINT 899
  899 FORMAT('LA TECHNIQUE UTILISEE EST EN DEFAUT')
      ITRY=1
  914 FORMAT(1H ,5(5X,E13.6))
      RETURN
36    APR=(RMAX-H)*DSQRT(GI)-RMAX*DSQRT(GN)
      IF(APR.GT.50.d0) APR=50.d0
      P(N-1)=P(N)*DEXP(-APR)
38    Y(1)=(1.d0-HV*GN)*P(N)
40    Y(2)=(1.d0-HV*GI)*P(N-1)
 
C     INTEGRATION
   44 M=N-2
   46 Y(3)=Y(2)+((Y(2)-Y(1))+H2*GI*P(M+1))
   48 GI=V(M)-E
50    P(M)=Y(3)/(1.d0-HV*GI)
52    IF(DABS(P(M)).LT.1.d+34) GO TO 70
 
C     DEPASSEMENT DE LA LIMITE
      M1=M+1
      PM=P(M1)
179   FORMAT(2X,'PM = ',E16.8/)
   55 DO 56 J=M1,N
   56 P(J)=P(J)/PM
   58 Y(1)=Y(1)/PM
C
C     NOUVEAU DEPART
   60 Y(2)=Y(2)/PM
   62 Y(3)=Y(3)/PM
      GI=V(M+1)-E
      GO TO 46
C     L-INTEGRATION VERS L-INTERIEUR EST-ELLE TERMINEE
   70 IF((DABS(P(M)).LE.DABS(P(M+1))).OR.(M.LE.2))GO TO 90
   81 Y(1)=Y(2)
   82 Y(2)=Y(3)
   84 M=M-1
      GO TO 46
C     L-INTEGRATION VERS L-INTERIEUR EST TERMINEE
   90 PM=P(M)
      MSAVE=M
   92 YIN=Y(2)/PM
   94 DO 96 J=M,N
   96 P(J)=P(J)/PM
C     INTEGRATION VERS L-EXTERIEUR.PREMIERS PAS
100   P(1)=1.d0
  102 Y(1)=0.d0
  104 GI=V(1)-E
106   Y(2)=(1.D0-HV*GI)*P(1)
C
C     INTEGRATION
  108 DO 132 I=2,M
  110 Y(3)=Y(2)+((Y(2)-Y(1))+H2*GI*P(I-1))
  112 GI=V(I)-E
114   P(I)=Y(3)/(1.d0-HV*GI)
116   IF(DABS(P(I)).LT.1.d+34) GO TO 130
C     LA LIMITE A ETE DEPASSEE
  118 I1=I-1
      PM=P(I1)
      DO 120 J=1,I1
  120 P(J)=P(J)/PM
  122 Y(1)=Y(1)/PM
  124 Y(2)=Y(2)/PM
  126 Y(3)=Y(3)/PM
      GI=V(I1)-E
      GO TO 110
  130 Y(1)=Y(2)
  132 Y(2)=Y(3)
C     L-INTEGRATION VERS L-EXTERIEUR EST TERMINEE
C
      PM=P(M)
      IF(PM) 135,149,135
  135 YOUT=Y(1)/PM
  136 YM=Y(3)/PM
  138 DO 140 J=1,M
  140 P(J)=P(J)/PM
C     LES DEUX BRANCHES SONT MAINTENANT RACCORDEES
C     CORRECTION
142   DF=0.d0
  144 DO 146 J=1,N
  146 DF=DF-P(J)**2
148   F=(-YOUT-YIN+2.d0*YM)/H2+(V(M)-E)
      DOLD=DE
      IF(DABS(F).LT.1.d+37) GO TO 150
149   F=9.99999d+29
      DF=-F
      DE=DABS(0.0001d0*E)
      GO TO 152
  150 DE=-F/DF
  152 CONTINUE
  156 FORMAT(I2,5X,E16.8,5X,E16.8,5X,E16.8,5X,E16.8)
  164 EOLD=E
      E=E+DE
      TEST=DMAX1((DABS(DOLD)-DABS(DE)),TEST)
      IF(TEST.LT.0.d0) GO TO 171
      IF(DABS(E-EOLD).LT.EPS)GO TO 172
  171 CONTINUE

      SCHROD=1.d0
      GO TO 173
C     LES ITERATIONS SONT TERMINEES
172   SCHROD=0.d0

C     LES ITERATIONS ONT CONVERGE
C     COMPTAGE DES NOEUDS
  173 KV=0
      NL=N-2
  174 DO 192 J=3,NL
  176 IF(P(J))178,177,177
  177 IF(P(J-1))180,192,192
  178 IF(P(J-1))192,270,184
C     NOEUD POSITIF
  180 IF(P(J+1))192,182,182
  182 IF(P(J-2))190,192,192
C     NOEUD NEGATIF
  184 IF(P(J+1))186,192,192
  186 IF(P(J-2))192,190,190
C     LE NOEUD EST-IL DU A UN SOUS-DEPASSEMENT
  270 IF(P(J+1))280,192,192
  280 IF(P(J-2))192,192,190
  190 KV=KV+1
  192 CONTINUE
      E2=E
C     NORMALISATION

  200 SN=DSQRT(-H*DF)
  202 DO 204 J=1,N
  204 P(J)=P(J)/SN
  250 FORMAT(10X,'SCHR=',I1,/)
C     ECRITURE DE LA SOLUTION
  214 FORMAT('V=',I3,5X,' E=',E16.8/)
  228 FORMAT(10X,E16.8,10X,E16.8)
      RETURN
      END

      SUBROUTINE BESJOT (L,X,F,DF,R)
C
C   THIS PROGRAM GENERATES THE STANDARD AND MODIFIED VERSIONS OF THE SPHERICAL
C   BESSEL-FUNCTIONS OF FIRST(BESJOT)- AND SECOND(BESSEN)-KIND RESPECTIVELY.
C   FOR DEFINITIONS COMPARE: 'NBS-H1NDBOOK OF MATHEMATICAL FUNCTIONS',
C   (ABRAMOWITZ+STEGUN,EDS./N.Y.:1964), SECTIONS 10.1.1 ON PAGE 437 FOR
C   STANDARD VERSIONS AND SS.10.2.2 + 10.2.3 ON P.443 FOR THE MODIFIED ONES.
C   L=INDEX(NATURAL NUMBERS INCLUDING ZERO), X=ARGUMENT(REAL,D.P.), F=OUTPUT.
C
C   THE SIGN OF THE ARGUMENT IS USED TO DETERMINE THE VERSIONS:
C   THE OUTCOMES F(=FIRST-KIND-FUNCTIONS) AND G(=SECOND-KIND-F.) MUST BE
C   DIVIDED (RESP. MULTIPLIED) BY THE L-TH POWER OF THE REDUCTION LOGFAC
C   TO GET THE MODIFIED VERSIONS, USE ARGUMENT WITH NEGATIVE SIGN }
C
C   BY FORMULAS 10.1.31 ON PAGE 439 LOC.CIT. AND 10.2.7 ON P.443 IBID.,
C   SOLUTIONS HAVE BEEN TESTED TO BE CORRECT TO TWELVE PLACES AT LEAST IN THE
C   RANGE COMBINING X=1...441 AND L=0...340 .
C
C   BESJOT IS DIVIDED INTO THREE PARTS, CORRESPONDING TO WETHER X > L, O
C   WHILE X < L, BEEING 0.5*X*X < 2*L OR 0.5*X*X > 2*L RESPECTIVELY .
C
C
C MODIF. POUR X PLUS GRAND QUE L DANS BESJOT  : R=1
C        QQ SOIT X DANS BESSEN : R=1
C   POUR EVITER LES OVERFLOWS OU UNDERFLOWS DANS LE PROG. APPELE POUR 50
C      VERSION JAN. 77
C
C  J.M.LAUNAY, MEUDON, FRANCE
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
       F = 0.D0
       R = 1.D0
      IF(X)  51,50,52
 50   IF(L.EQ.0)  F=1.D0
      DF = 0.D0
      RETURN
C
 51   SINIX = DSINH(-X)
      COSIX = DCOSH(-X)
       W = -1.D0
      GOTO 53
 52   PI = 6.283185307179586D0
      XR = DMOD(X,PI)
      SINIX = DSIN(XR)
      COSIX = DCOS(XR)
       W = +1.D0
 53    Z = 1.D0 / DABS(X)
       A = DBLE(L)
       R = A * Z
      IF(DABS(X)-A) 2,2,1
   2  IF(0.5D0*X*X-2.D0*A) 4,4,3
C
C   FOR THE FOLLOWING VERSION SEE PAGE 439, SECTION 10.1.19 LOC.CIT., AN
C   SECTION 10.1.11 ON PAGE 438 IBIDEM}
C   FOR THE MODIFIED CASE LOOK UP SS.10.2.18 AND 10.2.13 ON PS.444 AND 4
C   THIS VERSION IS USED, IF X > L .
C
  1   F0 = SINIX * Z
      G0 =-COSIX * Z
        R=1.D0
      IF(L) 11,11,12
 11    F = F0
      DF =-G0 - F0*Z
      RETURN
 12   IF(L-1) 13,13,14
 13    F = W * (F0-COSIX) * Z
      DF = F0   - 2.0D0*F*Z
      RETURN
 14   F1 = W * (F0-COSIX) * Z
      IF(L.EQ.2)  GOTO 15
       J = L-2
      DO 10 I=1,J
      F2 = W * (F1*DBLE(2*I+1)*Z - F0  )
      F0 = F1
   10 F1 = F2
   15  F = W * (F1*DBLE(2*L-1)*Z - F0  )
      DF = F1   - (L+1)*F*Z
      RETURN
C
C   FOR THE FOLLOWING VERSION SEE PAGE 453, EXAMPLE 2 LOC.CIT.}
C   THIS VERSION IS USED, IF X < L AND 0.5*X*X > 2*L .
C
  3    N = A + 25.D0 + DSQRT(A)
      B0 = 0.D0
      B1 = 1.D0
      DO 20 J=1,N
      B2 = W * (B1*DBLE(2*(N-J)+3)*Z - B0/R) / R
      B0 = B1
      IF(N-L-J) 22,21,22
 21    F = B2
      GOTO 20
 22   IF(N-L+1-J)  20,23,20
 23   DF = B2
   20 B1 = B2
      DF = W**(L-1) * (DF/B1) * SINIX*Z
       F = W**L * (F/B1) * SINIX*Z
      DF = DF*R - (L+1)*F*Z
      RETURN
C
C   FOR THE FOLLOWING VERSION SEE P1GE 437, FORMULA 10.1.2  LOC.CIT.}
C   FOR THE MODIFIED CASE FORMULA 10.2.5 ON PAGE 443 IS VALID.
C   THIS VERSION IS USED, IF X < L AND 0.5*X*X < 2*L .
C
C
  4    Y = -W * 0.5D0 * X * X
      S0 = DBLE(2*L-1)
      S1 = DBLE(2*L+1)
      P0 = 1.D0
      P1 = 1.D0
      C0 = 1.D0
      C1 = 1.D0
      DO 30 I=1,15
      S0 = S0 + 2.D0
      S1 = S1 + 2.D0
      P0 = Y*P0/(S0*DBLE(I))
      P1 = Y*P1/(S1*DBLE(I))
      C0 = C0 + P0
   30 C1 = C1 + P1
       Q = 1.D0
      IF(L.EQ.1)  GOTO 32
       J = L - 1
      DO 31 I=1,J
   31  Q = Q * A/DBLE(2*I+1)
       F = Q * A/DBLE(2*L+1) * C1
      DF = Q*C0*R - DBLE(L+1)*F*Z
      RETURN
  32   F = C1 / 3.D0
      DF = C0*R - 2.D0*F*Z
      RETURN
C
C
           ENTRY BESSEN(L,X,G,DG,R)
C
C   SPHERICAL BESSEL-(AND MODIFIED BESSEL-) FUNCTIONS OF THE SECOND KIND
C   THIS VERSION IS VALID FOR ALL INDICES AND ARGUMENTS.
C
       G = 0.D0
       A = DBLE(L)
       R = 1.D0
      IF(X) 61,60,62
 60   WRITE(6,300)
300   FORMAT(1H0,'*******  ARGUMENT OF SPHERICAL BESSEL-FUNCTION OF SECO
     1ND KIND SHOULD NOT BE ZERO }')
      RETURN
 61   SINIX = DSINH(-X)
      COSIX = DCOSH(-X)
       W = -1.D0
      GOTO 63
 62   PI = 6.283185307179586D0
      XR = DMOD(X,PI)
      SINIX = DSIN(XR)
      COSIX = DCOS(XR)
       W = +1.D0
 63    Z = 1.D0 / DABS(X)
      G0 =-W * COSIX * Z
      F0 = W * SINIX * Z
      IF(L) 41,41,42
 41    G = G0
      DG = W*F0 - G0*Z
      RETURN
 42        R=1.D0
       IF(L-1) 43,43,44
 43    G = W * (G0-SINIX)*Z
      DG = G0   - 2.0D0*G*Z
      RETURN
 44   G1 = W * (G0-SINIX)*Z
      IF(L.EQ.2)  GOTO 45
       J = L-2
      DO 40 I=1,J
      G2 = W * (G1*DBLE(2*I+1)*Z - G0  )
      G0 = G1
   40 G1 = G2
  45   G = W * (G1*DBLE(2*L-1)*Z - G0  )
      DG = G1   - DBLE(L+1)*G*Z
      RETURN
      END
*******************************  BESSEL  ***********************************

      SUBROUTINE BESPH2 (F,DF,G,DG,AA,ARG,KEY,IBUG)
C#######################################################################
C#    CALCULATES LINEARLY INDEPENDANT SOLUTIONS OF THE EQUATION        #
C#       2    2        2                                               #
C#    ( D / DX  - A / X  + 1 ) Y(X) = 0                                #
C#                       -                                             #
C#    WHERE A = L*(L+1) WITH L INTEGER .GE. 0                          #
C#    + (-) SIGN CORRESPONDS TO AN OPEN (CLOSED) CHANNEL               #
C#    THE SOLUTIONS ARE OBTAINED FROM BESSEL FUNCTIONS OBTAINED        #
C#    IN BESJOT, BESSEN, BESSIK SUBROUTINES                            #
C#    SEE ABRAMOWITZ AND STEGUN (CHAP. 10)                             #
C#    ASYMPTOTIC BEHAVIOUR IS :                                        #
C#    F # SIN (X-L*PI/2) ; G # -COS (X-L*PI/2) FOR OPEN   CHANNELS     #
C#    F # SINH (X); G # EXP (-X)               FOR CLOSED CHANNELS     #
C#---------------------------------------------------------------------#
C#    AA   : L*(L+1)                                                   #
C#    ARG  : ARGUMENT VALUE                                            #
C#           IF POSITIVE THEN ARGUMENT Z IS REAL      (Z = ARG)        #
C#           IF NEGATIVE THEN ARGUMENT Z IS IMAGINARY (Z = -I*ARG)     #
C#           WHERE I = (-1)**0.5                                       #
C#    F,G,DF,DG : REGULAR AND IRREGULAR FUNCTIONS AND THEIR DERIVATIVES#
C#    KEY  : .LT.0 TO SUPPRESS EXPONENTIAL FACTORS IN BESSIK          #
C#    IBUG : .GT. 0 TO PRINT THE OUTPUT                                #
C#---------------------------------------------------------------------#
C#    HAVE BEEN TESTED ON WRONSKIAN RELATION W(F,G) = F*DG-DF*G = 1    #
C#    FOR THE FOLLOWING VALUES OF THE ARGUMENT AND ORDERS :            #
C#    0 =< L  < 100  AND 0 =< Z  < 100   (ERROR IN W LESS THAN 10**-12)#
C#    0 =< L  < 30   AND 0 =< Z  < 100*I (ERROR IN W LESS THAN 10**-10)#
C#    OUTSIDE THIS RANGE CHECK FOR UNDERFLOWS, OVERFLOWS, DIVIDE CHECKS#
C#    AND DEXP CAPACITY.                                               #
C#---------------------------------------------------------------------#
C#    J.M.L. 08/1981 ; UPDATE : 12/1981                                #
C#    J.M.LAUNAY, MEUDON, FRANCE                                       #
C#######################################################################
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FF(2),DFF(2)
      COMMON /TBESP2/ TIME(2,2)
C     DATA PI /3.1415 92653 58979 32384 D0/,TINY /1.D-10/
      pi=dacos(-1.d0)
      tiny=1.d-15
C
      IF (AA .LT. -0.25D0) GO TO 10
      FL = -0.5D0+DSQRT(AA+0.25D0)+TINY
      L = FL
      IF (DABS(FL-L) .LT. 2.D0*TINY) GO TO 20
 10   PRINT 9010,AA
      RETURN
C
 20   X = DABS(ARG)
      SK2 = -1.
      IF (ARG .GT. 0.D0) GO TO 100
      TWOPIM = 2.D0/PI
      LP = L+1
      LM = L-1
      CALL BESSIK (L ,X,BI ,BK ,KEY)
      CALL BESSIK (LM,X,BIM,BKM,KEY)
      CALL BESSIK (LP,X,BIP,BKP,KEY)
      F  = -BI*X
      DF = -( (L*BIM+LP*BIP)/(L+LP)*X+BI)
      G  =  BK*X*TWOPIM
      DG =  (-(L*BKM+LP*BKP)/(L+LP)*X+BK)*TWOPIM
      GO TO 1000
C
 100  CALL BESJOT (L,X,BJ,DBJ,R)
      CALL BESSEN (L,X,BN,DBN,RG)
      SK2 = 1.
      IF (L .EQ. 0)  GO TO 200
      DO 210 I = 1,L
         BJ = BJ/R
         DBJ = DBJ/R
         BN = BN*RG
         DBN = DBN*RG
 210  CONTINUE
 200  F = X*BJ
      G = X*BN
      DF = X*DBJ+BJ
      DG = X*DBN+BN
C
C --- TIME SURFACE INTEGRALS CALCULATION
C
 1000 FF(1) = F
      FF(2) = G
      DFF(1) = DF
      DFF(2) = DG
      VV = AA/(X*X)-SK2
C     DO 1010 I = 1,2
C        DO 1011 J = 1,2
C           TIME(I,J) = -0.125*(FF(I)*DFF(J)+DFF(I)*FF(J)
C    &                         +2.*X*(FF(I)*VV*FF(J)-DFF(I)*DFF(J)) )
 1011    CONTINUE
 1010 CONTINUE
      IF (IBUG .LE. 0) RETURN
C
 2000 WM1 = F*DG-DF*G-1.D0
C      PRINT 9000,AA,ARG,F,DF,G,DG,WM1
      RETURN
C
 9000 FORMAT (' BESPH2 ',1F12.4,1F12.4,1P,4D20.12,1D9.1)
 9010 FORMAT (' ****** BESPH2 ROUTINE; AA = ',F12.2,' IS NOT L*(L+1)'
     &       ,' WITH L INTEGER .GE. 0; REQUEST ABORTED ******')
      END
      SUBROUTINE BESSIK (L,X,BI,BK,KEY)
C#######################################################################
C#    CALCULATION OF MODIFIED SPHERICAL BESSEL FUNCTIONS OF THE THIRD  #
C#    KIND  :                                                          #
C#    (PI/(2*X))**1/2 * I     (X)  BY FORMULA 10.2.5                   #
C#                       L+1/2     FROM THE HANDBOOK (P.443)           #
C#    (PI/(2*X))**1/2 * K     (X)  BY FORMULA 10.2.15                  #
C#                       L+1/2     FROM THE HANDBOOK (P.444)           #
C#---------------------------------------------------------------------#
C#    L     : ORDER (INTEGER)                                          #
C#    X     : ARGUMENT (POSITIVE REAL NUMBER)                          #
C#    BI,BK : BESSEL FUNCTIONS  ;                                      #
C#    KEY  =0: NORMAL, >0: MULTIPLIES BI BY EXP(X), BK BY EXP(-X)      #
C#                     <0:     ''     BI BY EXP(-X), BK BY EXP(X)      #
C#######################################################################
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /LGFAC/ FCT(50000)
      DATA NTIME /0/
      PI=DACOS(-1.D0)
      FACTOR=DLOG(1.D-30)
C
      NTIME  = NTIME+1
      IF (NTIME .EQ. 1) CALL FACLG
      BI = 0.D0
      BK = 0.D0
      IF (L .LT. 0) RETURN
      FL = L
      S2 = DLOG(2.D0)
      SX = DLOG(X)
      HOX = 0.5D0/X
      CI = FL*SX
      SHX2 = DLOG(0.5D0)+2.D0*SX
      ALFA=0.D0
      IF(KEY.LT.0)ALFA=X
      IF(KEY.GT.0)ALFA=-X
      BETA=X-ALFA

C  COMPUTATION OF BK

      SUM = HOX*DEXP(-BETA)
      IF (L .EQ. 0) GO TO 100
         SHOX = DLOG(HOX)
         SHOX1 = 0.D0
         KMIN = 0
         KMAX = L
         SUM = 0.D0
C
         DO 10 K = KMIN,KMAX
            SHOX1 = SHOX+SHOX1
            ST = FCT(L+K+1)-FCT(K+1)-FCT(L-K+1)+SHOX1
            IF (ST .LT. -180.D0) GO TO 100
            SUM = SUM+DEXP(ST-BETA)
 10      CONTINUE
 100  BK = PI*SUM

C  COMPUTATION OF BI

      SUM = 0.D0
      K = 0
      LK = L
 200  CONTINUE
         ST = K*SHX2-FCT(K+1)-( FCT(2*LK+2)-LK*S2-FCT(LK+1) )
         K = K+1
         LK = L+K
         IF(2*LK+2.GT.5000)GO TO 9990
         TERM=DEXP(ST-ALFA+CI)
         SUM=SUM+TERM
         SSUM=0.D0
         IF(SUM.GT.0.D0)SSUM=DLOG(SUM)
      IF(SSUM+ALFA.LT.FACTOR)GO TO 200
      IF(SSUM.EQ.0.D0)GO TO 200
      IF((TERM/SUM).GT.1.D-20)GO TO 200

      BI=SUM

      RETURN
 9990 CONTINUE
      PRINT *,'***> ERROR IN BESSIK: MAX. VALUE OF FACTORIAL=',2*LK+2
      STOP
      END

C***********************************************************************
      SUBROUTINE JORDAN(MU,LAMBDA,X,N,B)
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,L-M,O-Z)
C
      PARAMETER(NX2MAX=20000)
C
      DIMENSION MU(N),LAMBDA(N),X(N),B(N)
      DIMENSION PIV(NX2MAX)
C
      IF(N.GT.NX2MAX) GO TO 999
C
C
C     CALCUL DES PIVOTS
      PIV(1)=2.D0
      DO 10 I=2,N
      PIV(I)=2.D0-LAMBDA(I)*MU(I-1)/PIV(I-1)
   10 B(I)=B(I)-LAMBDA(I)/PIV(I-1)*B(I-1)
C
      X(N)=B(N)/PIV(N)
      I=N-1
   20 X(I)=(B(I)-X(I+1)*MU(I))/PIV(I)
      I=I-1
      IF(I.GT.0) GOTO 20
      RETURN
C
 999  CONTINUE
      PRINT 9000
      PRINT 9999, NX2MAX, N
      STOP
 9000 FORMAT(//,2X,20('*'),' STOP IN JORDAN ',20('*'),/)
 9999 FORMAT(2X,'DIMENSION PARAMETER NX2MAX TOO SMALL , ',I5,
     &       'REQUIRED')
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION DLAGRA(X,Y,MIN,IP)
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MIN),Y(MIN)
      DLAGRA=0.D0
      DO 10 I=1,MIN
      IF(I.EQ.IP) GOTO 10
      YP=Y(I)
      DO 20 J=1,MIN
         IF(J.EQ.IP) GOTO 20
         IF(J.EQ.I) GOTO 20
            YP=YP*(X(IP)-X(J))
   20 CONTINUE
      DO 30 J=1,MIN
         IF(J.EQ.I) GOTO 30
            YP=YP/(X(I)-X(J))
   30 CONTINUE
      DLAGRA=DLAGRA+YP
   10 CONTINUE
      DO 40 I=1,MIN
         IF(I.EQ.IP) GOTO 40
            DLAGRA=DLAGRA+Y(IP)/(X(IP)-X(I))
   40 CONTINUE
      RETURN
      END
      SUBROUTINE FACLG
C#######################################################################
C#    INITIALISATION OF LOGARITHMS OF FACTORIALS ARRAY                 #
C#######################################################################
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /LGFAC/ FCT(50000)
      DATA NTIMES /0/
C
      NTIMES = NTIMES+1
      IF (NTIMES .GT. 1) RETURN
      FCT(1) = 0.D0
      DO 10 I = 1,4999
         AI = I
         FCT(I+1) = FCT(I)+DLOG(AI)
 10   CONTINUE
C
      RETURN
      END
      subroutine sinmom(box,npun,npundim,xmred,hbr,pr,p2r)
      implicit real*8(a-h,o-z)

      dimension pr(npundim),p2r(npundim)

      if(npun.gt.npundim)then
        write(6,*)'  npun= ',npun, ' > npundim = ',npundim,' in fftmom'
        stop
      endif   
         
      pi = dacos(-1.d0)
      dpi = 2.d0*pi
      hbrxm = 0.5d0*hbr*hbr/xmred
      ah = (box)/dble(npun-1)
      box = 2.d0*(dble(npun)+1.d0)*ah
      do ir=1,npundim
         pr(ir)=0.d0
         p2r(ir)=0.d0
      enddo

      do ir=1,npun
         iii=ir
         p = dpi*dble(iii)/box
         pr(ir) = p
         p2r(ir) = p*p*hbrxm 
      enddo

      return
      end
      subroutine noptFFT(nin,nout,ntot)
      integer nin,nout,ntot,i,nmax,imin,idis
      integer nexp2,nexp3,nexp5,nexp7,nexp11
      integer n2,n3,n5,n7,n11
      integer i2,i3,i5,i7,i11

      nmax=10000
      if(nin.gt.nmax.or.ntot.gt.nmax)then
          write(6,*)' nin= ',nin,' , ntot= ',ntot,'  > nmax= ',nmax
          write(6,*)'   in noptFFT of libdyn library'
          stop
      endif

      nout=1
      i=1
      
      nexp2=1
      do n2=0,100
         nexp2=nexp2*2
         if(nexp2.gt.ntot)go to 2
      enddo
 2    continue
      nexp2=n2

      nexp3=1
      do n3=0,100
         nexp3=nexp3*3
         if(nexp3.gt.ntot)go to 3
      enddo
 3    continue
      nexp3=n3

      nexp5=1
      do n5=0,100
         nexp5=nexp5*5
         if(nexp5.gt.ntot)go to 5
      enddo
 5    continue
      nexp5=n5

      nexp7=1
      do n7=0,100
         nexp7=nexp7*7
         if(nexp7.gt.ntot)go to 7
      enddo
 7    continue
      nexp7=n7

      nexp11=1
      do n11=0,100
         nexp11=nexp11*11
         if(nexp11.gt.ntot)go to 11
      enddo
 11   continue
      nexp11=n11

      i=1
      imin=ntot

      do n11=0,nexp11
         i11=11**n11
         do n7=0,nexp7
            i7=7**n7
            do n5=0,nexp5
               i5=5**n5
               do n3=0,nexp3
                  i3=3**n3
                  do n2=0,nexp2
                     i2=2**n2
                     i=i2*i3*i5*i7*i11
                     idis=i-nin                     
                     if(idis.ge.0.and.idis.lt.imin)then
                         nout=i
                         imin=idis
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      if(nout.gt.ntot)nout=ntot

      return
      end
      subroutine fftmom(box,npun,npundim,xmred,hbr,pr,p2r)
      implicit real*8(a-h,o-z)

      dimension pr(npundim),p2r(npundim)

      if(npun.gt.npundim)then
        write(6,*)'  npun= ',npun, ' > npundim = ',npundim,' in fftmom'
        stop
      endif   
         
      pi = dacos(-1.d0)
      dpi = 2.d0*pi
      hbrxm = 0.5d0*hbr*hbr/xmred
      ah = (box)/dble(npun-1)
      box = ah*dble(npun)
      do ir=1,npundim
         pr(ir)=0.d0
         p2r(ir)=0.d0
      enddo

      do ir=1,npun
         if(ir.le.(npun/2))then
            iii=ir-1
         else
            iii=ir-npun-1
         endif
         p = dpi*dble(iii)/box
         pr(ir) = p
         p2r(ir) = p*p*hbrxm 
      enddo

      return
      end
*************************  tresj  **************************

      subroutine tresj(j1,j2,j3,m1,m2,m3,coef)
      implicit real*8 (a-h,o-z)

      common/fct/fact(0:10000)

      coef=0.d0

      if(m1+m2+m3.eq.0.and.trian(j1,j2,j3).ne.0.d0)then

          b=fact(j1+m1)+fact(j1-m1)+fact(j2+m2)
     &     +fact(j2-m2)+fact(j3+m3)+fact(j3-m3)
          b=0.5d0*b

          k1=j3-j2+m1
          k2=j3-j1-m2
          k3=j1+j2-j3
          k4=j1-m1
          k5=j2+m2

          kmin=max0( -k1 , -k2 , 0 )
          kmax=min0( k3 , k4 , k5 )

             isign=-1
             if(mod(kmin,2).eq.0)isign=1

             do 10 k=kmin,kmax 
               a=fact(k)+fact(k3-k)+fact(k4-k)
     &          +fact(k5-k)+fact(k1+k)+fact(k2+k)
               coef=coef+isign*dexp(b-a)
               isign=isign*(-1) 
 10          continue

             isign=-1
             if(mod(j1-j2-m3,2).eq.0)isign=1
             coef=coef*trian(j1,j2,j3)*isign

      endif
   
      return
      end
*******************************  NPLEGM  ********************************

      subroutine nplegm(pm,lmax,m,y,ndim)
      implicit real*8(a-h,o-z)
*     ******************************************
*     **  normalized                          **
*     **     Associated Legedre functions     **
*     **                                      **
*     ** Recursion formula (Num.Rec. pg.182)  **
*     **   P_m^m (y) =(2m-1)!! (1-y^2)^(m/2)  **
*     **         P_{m-1}^m (y)=0              **
*     **   (l-m)P_l^m = y(2l-1)P_{l-1}^m      **
*     **              - (l+m-1)P_{l-2}^m      **
*     **                                      **
*     **  Input:                              **
*     **      lmax: maximum l value desired   **
*     **      m   : m value                   **
*     **      y   : argument (cos(theta)      **
*     **  Output:                             **
*     **      pm  : array with the values..   **
*     **            pm(0)= Y_0^m              **
*     **            pm(1)= Y_1^m              **
*     ****************************************** 

      dimension pm(0:ndim)

      if(lmax.gt.ndim)then
           write(6,*)' ** be carefull with dimensions in nplegm **'
           write(6,*)'     lmax > ndim',lmax,ndim
           stop
      endif
      if(dabs(y).gt.1.d0)then
           write(6,*)' ** Bad argument in plegm **'
           stop
      elseif(m.lt.0)then
           write(6,*)' ** m.lt.0 in plegm **'
           stop
      endif 

      do l=0,lmax
         pm(l)=0.d0
      enddo
      faclog=0.d0

***>> Compute P_m^m and P_{m+1}^m

      if(m.eq.0)then
         pm(m)=1.d0
         pm(m+1)=y
      else
         faclog=0.d0
         do i=1,2*m-1,2
            faclog=faclog+dlog(dble(i))
         enddo
          xxx=dsin(dacos(y))
          xxx=xxx**(dble(m))
          xxx=xxx*(-1.d0)**(m)

c         xxx=(1-y*y)
c         xxx=xxx**(0.5d0*dble(m))
         pm(m)=xxx
         pm(m+1)=y*dble(2*m+1)*pm(m)
      endif

***>> Compute P_l^m  l=m+2,...,lmax

      do l=m+2,lmax
         xlm=dble(l-m)
         pm(l)=y*dble(2*l-1)*pm(l-1)-dble(l+m-1)*pm(l-2)
         pm(l)=pm(l)/xlm
      enddo

***>> normalization

      do l=0,lmax
        facden=0.d0
        if(l+m.gt.0)then
           do i=1,l+m
              facden=facden+dlog(dble(i))
           enddo
        endif
        facnum=0.d0
        if(l-m.gt.0)then
           do i=1,l-m
              facnum=facnum+dlog(dble(i))
           enddo
        endif
        factor=dexp(0.5d0*(facnum-facden)+faclog)
        if(m.eq.0)factor=1.d0
        factor=factor*dsqrt(dble(l)+0.5d0)
        pm(l)=pm(l)*factor
      enddo

      return
      end

*************************  coefread  ***************************

      subroutine coefread(ifile,nmin,nmax,ntot,nbaslie
     &              ,cdis,jd,ld,nv1,nv2,eee)
      implicit real*8(a-h,o-z)
      character*40 fcoef

      dimension cdis(nbaslie,nmin:nmax)
      dimension jd(nbaslie),ld(nbaslie),nv1(nbaslie),nv2(nbaslie)

c      open(ifile,file=fcoef,status='old')
       read(ifile,*)numvec
       if(numvec.lt.nmax)then
          write(6,*)'  ** In coefread: no. of eigenvectors= ',numvec
          write(6,*)'      while desired state are between= ',nmin,nmax
          stop 
       endif

       read(ifile,*)nnn,ntot
       if(ntot.gt.nbaslie)then
          write(6,*)'  ** In coefread: no. of basis fucntions= ',ntot
          write(6,*)'     larger than dimension nbaslie= ',nbaslie
          stop 
       endif
      
* Reading basis set quantum numbers

       do i=1,ntot
          read(ifile,*)nv1(i),ld(i),nv2(i),iii,jd(i)
       enddo

* Reading coefficients

       do ibound=1,nmax
           read(ifile,*)kkkk,eee
           write(6,*)'    ',ibound,kkkk,' bound state energy = ',eee
           call flush(6)
           do i=1,ntot
              read(ifile,*)ccc
              if(ibound.ge.nmin)cdis(i,ibound)=ccc
           enddo
           if(ibound.ge.nmin)then
              write(6,"('First and last coef. of wvf ='
     &          ,e15.7,2x,e15.7)")cdis(1,ibound),cdis(ntot,ibound)
           endif

* Distributions of bound  state

         if(ibound.ge.nmin)then
            write(6,*)'  v distribution '
            m1=nv1(1)
            m2=nv1(ntot)
            do im= m1,m2
               xxx=0.d0
               do i=1,ntot
               if(nv1(i).eq.im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
               enddo
               write(6,*)'      v=',im,' -->',xxx
            enddo
            write(6,*)'  j distribution '
            m1=jd(1)
            m2=jd(ntot)
            do im= m1,m2
               xxx=0.d0
               do i=1,ntot
               if(jd(i).eq.im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
               enddo
               write(6,*)'      j=',im,' -->',xxx
            enddo
            write(6,*)'  n distribution '
            m1=nv2(1)
            m2=nv2(ntot)
            do im= m1,m2
               xxx=0.d0
               do i=1,ntot
               if(nv2(i).eq.im)xxx=xxx+cdis(i,ibound)*cdis(i,ibound)
               enddo
               write(6,*)'      n=',im,' -->',xxx
            enddo
 
          endif
       enddo

c      close(ifile)

       return
       end
******************************  rfunread  **********************************

       subroutine rfunread(ifile,npunr,nvini,nvmax,fread,nuno,
     &                    fspl,npun,rmis,ah,f,x)
       implicit real*8(a-h,o-z)
       dimension fread(npunr,nvini:nvmax),fspl(npun,nvini:nvmax)
       dimension f(npunr,2),x(npunr)
  
       CONVL=.52917726D0
       
       write(6,*)ifile,npunr,nvini,nvmax,nuno,npun,rmis,ah
       call flush(6) 
       read(ifile,*)nnn,nuno,nvinir,nvmaxr
       write(6,*)nnn,nuno,nvinir,nvmaxr
       call flush(6)
       if(nuno.eq.0)then
          if(nnn.gt.npunr)then
             write(6,*)' *** Attention in nfunread ifile= ',ifile,' ***'
             write(6,*)'   number of points in file ',ifile,' = ',nnn
             write(6,*)'   larger than npunlie= ',npunr
             stop
          endif
          if(nvini.gt.nvinir.or.nvmaxr.gt.nvmax)then
             write(6,*)' *** Attention in nfunread ***'
             write(6,*)'   number of functions in file',ifile,'= '
     &                   ,nvinir,nvmaxr
           write(6,*)' are not compatible with dimensions= ',nvini,nvmax
             stop
          endif
         
          do ir=1,nnn
             read(ifile,*)x(ir),(fread(ir,iv),iv=nvinir,nvmaxr)
             x(ir)=x(ir)*convl
          enddo                  
       
          do iv=nvinir,nvmaxr
             xnorm=0.d0
             do ir=1,nnn
                f(ir,1)=fread(ir,iv)/dsqrt(convl)
                f(ir,2)=0.d0
             enddo

             call splset(f,x,nnn,npunr)
             iold=2
             do ir=1,npun
                r=rmis+dble(ir-1)*ah
                if(r.lt.x(1))then
                   fspl(ir,iv)=0.d0
                elseif(r.gt.x(nnn))then
                   fspl(ir,iv)=0.d0
                else
                   call splinqq(f,x,iold,nnn,r,npunr,spl)
c                  fspl(ir,iv)=splinq(f,x,iold,nnn,r,npunr)
                   fspl(ir,iv)=spl
                endif
                xnorm=xnorm+fspl(ir,iv)*fspl(ir,iv)
             enddo
             write(6,*)iv,xnorm*ah
          
          enddo
      endif

      return
      end
c********************  trian   ****************************

      function trian(j1,j2,j3)
      implicit real*8 (a-h,o-z)

      common/fct/fact(0:10000)

      trian=0.d0

      if(j3.ge.iabs(j1-j2).and.j3.le.(j1+j2))then
         cc = fact(j1+j2-j3)+fact(j1-j2+j3)
     &       +fact(-j1+j2+j3)-fact(j1+j2+j3+1)
         trian = dexp(cc/2.d0)
      endif

      return
      end 

************************************************************************

      double precision function tresjd(j1d,j2d,j3d,m1d,m2d,m3d)
      implicit real*8(a-h,o-z)

      fj1=dble(j1d)*0.5d0
      fj2=dble(j2d)*0.5d0
      fj3=dble(j3d)*0.5d0
      fm1=dble(m1d)*0.5d0
      fm2=dble(m2d)*0.5d0
      fm3=dble(m3d)*0.5d0
     
      tresjd=f3j(fj1,fj2,fj3,fm1,fm2,fm3)

      return
      end function tresjd
!---------------------------------------------------------
      double precision function seisjd(j1d,j2d,j3d,j4d,j5d,j6d)
      implicit real*8(a-h,o-z)

      fj1=dble(j1d)*0.5d0
      fj2=dble(j2d)*0.5d0
      fj3=dble(j3d)*0.5d0
      fj4=dble(j4d)*0.5d0
      fj5=dble(j5d)*0.5d0
      fj6=dble(j6d)*0.5d0
     
      seisjd=f6j(fj1,fj2,fj3,fj4,fj5,fj6)

      return
      end function seisjd

c-----------------------------------------------------------------------
      double precision function f3j (fj1,fj2,fj3, fm1,fm2,fm3)
c#######################################################################
c#    calculates 3j coefficients from racah formula                    #
c#    (messiah: t2, p 910; formula 21) .                               #
c#    clebsch-gordan coefficients are given by (p. 908, formula 12) :  #
c#                         j -j +m                |j    j     j|       #
c#    <j j m m |j m> = (-1) 1  2   (2*j+1)**(0.5) | 1    2     |       #
c#      1 2 1 2                                   |m    m    -m|       #
c#                                                | 1    2     |       #
c#######################################################################
      implicit double precision (a-h,o-z)
      integer t,tmin,tmax
      parameter (nfctmx=5001)
      data tiny,zero,one /0.01d0,0.d0,1.d0/ ,ntimes /1/
      common /cfaclog/ fct(nfctmx)
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      cc = zero
      if (fj3 . gt. (fj1+fj2+tiny))      go to 100
      if (dabs(fj1-fj2) .gt. (fj3+tiny)) go to 100
      if (dabs(fm1+fm2+fm3) .gt. tiny)   go to 100
      if (dabs(fm1) .gt. (fj1+tiny))     go to 100
      if (dabs(fm2) .gt. (fj2+tiny))     go to 100
      if (dabs(fm3) .gt. (fj3+tiny))     go to 100
      fk1 = fj3-fj2+fm1
      fk2 = fj3-fj1-fm2
      fk3 = fj1-fm1
      fk4 = fj2+fm2
      fk5 = fj1+fj2-fj3
      fk1m = fk1-tiny
      fk2m = fk2-tiny
      fk1p = fk1+tiny
      fk2p = fk2+tiny
      if (fk1m .lt. zero) k1 = fk1m
      if (fk1p .gt. zero) k1 = fk1p
      if (fk2m .lt. zero) k2 = fk2m
      if (fk2p .gt. zero) k2 = fk2p
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = fj1+fj2-fj3+one+tiny
      n2 = fj2+fj3-fj1+one+tiny
      n3 = fj3+fj1-fj2+one+tiny
      n4 = fj1+fm1+one+tiny
      n5 = fj2+fm2+one+tiny
      n6 = fj3+fm3+one+tiny
      n7 = fj1-fm1+one+tiny
      n8 = fj2-fm2+one+tiny
      n9 = fj3-fm3+one+tiny
      n10 = fj1+fj2+fj3+2.d0+tiny
      x = fct(n1)+fct(n2)+fct(n3)+fct(n4)+fct(n5)+fct(n6)
     &   +fct(n7)+fct(n8)+fct(n9)-fct(n10)
      x = 0.5d0*x
      do 10  t = tmin,tmax
	 phase = one
	 if (mod(t,2) .ne. 0) phase = -one
	 cc = cc+phase*dexp(-fct(t+1)   -fct(k1+t+1)-fct(k2+t+1)
     &                      -fct(k3-t+1)-fct(k4-t+1)-fct(k5-t+1)+x)
 10   continue
      fsp = dabs(fj1-fj2-fm3)+tiny
      ns = fsp
      if (mod(ns,2) .gt. 0) cc = -cc
 100  f3j = cc
      return
      end
c-----------------------------------------------------------------------
      double precision function f6j (fj1,fj2,fj3,fl1,fl2,fl3)
c#######################################################################
c#    calculation of 6j-coefficients                                   #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=5001)
      common /cfaclog/ fct(nfctmx)
      data tiny /.01/ ,ntimes /1/
c
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      d = fdelta (fj1,fj2,fj3)
      d = d*fdelta (fj1,fl2,fl3)
      d = d*fdelta (fl1,fj2,fl3)
      d = d*fdelta (fl1,fl2,fj3)
      f6j = 0.d0
      if (dabs(d) .eq. 0.d0) return
c
      fk1 = fj1+fj2+fj3
      fk2 = fj1+fl2+fl3
      fk3 = fl1+fj2+fl3
      fk4 = fl1+fl2+fj3
      fk5 = fj1+fj2+fl1+fl2
      fk6 = fj2+fj3+fl2+fl3
      fk7 = fj3+fj1+fl3+fl1
      fmin = dmin1 (fk5,fk6,fk7)
      fmax = dmax1 (fk1,fk2,fk3,fk4)
      min = fmin+tiny
      max = fmax+tiny
      k1 = fk1+tiny
      k2 = fk2+tiny
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      k6 = fk6+tiny
      k7 = fk7+tiny
      if (min-max) 1000,3,3
 3    if (max) 1000,4,4
 4    if (min) 1000,5,90
 5    k1 = -k1
      k2 = -k2
      k3 = -k3
      k4 = -k4
      bot = fct(k1+1)+fct(k2+1)+fct(k3+1)+fct(k4+1)+fct(k5+1)+fct(k6+1)
     &     +fct(k7+1)
      bot = dexp(bot)
      f6j = d/bot
      return
c
 90   f6j = 0.
      do 100 i = max,min
	 boite = ((-1.)**i)
	 iz = i+1
	 m1 = i-k1
	 m2 = i-k2
	 m3 = i-k3
	 m4 = i-k4
	 m5 = k5-i
	 m6 = k6-i
	 m7 = k7-i
	 dot = fct(iz+1)
	 bot = fct(m1+1)+fct(m2+1)+fct(m3+1)+fct(m4+1)+fct(m5+1)
     &        +fct(m6+1)+fct(m7+1)
	 b1 = dot-bot
	 boite = boite*dexp(b1)
	 f6j = f6j+boite
 100  continue
      f6j = f6j*d
 1000 return
c
      end
c-----------------------------------------------------------------------
      double precision function fdelta (fl1,fl2,fl3)
      implicit double precision (a-h,o-z)
      parameter (nfctmx=5001)
      common /cfaclog/ fct(nfctmx)
      data   eps /.01/
c
      ia=fl1+fl2+fl3 +eps
      a=2.*(fl1+fl2+fl3)+1.
      ib=a +eps
      ib=ib/2
      if(ib-ia)1,6,1
    6 continue
      ik1=fl1+fl2-fl3+eps
      ik2=fl2+fl3-fl1+eps
      ik3=fl3+fl1-fl2+eps
      kk=fl1+fl2+fl3+1+eps
      if(ik1)1,2,2
    2 if(ik2)1,3,3
    3 if(ik3)1,4,4
    4 d1=fct(kk+1)
      d2=fct(ik1+1)+fct(ik2+1)+fct(ik3+1)
      d3 = (d2 - d1) / 2.d0
      fdelta = dexp (d3)
      go to 5
    1 fdelta=0.
    5 return
c
      end
      
      subroutine faclog
c#######################################################################
c#    initialisation of logarithms of factorials array                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=5001)
      common /cfaclog/ fct(nfctmx)
      data ntimes /0/
c
      ntimes = ntimes+1
      if (ntimes .gt. 1) return
      fct(1) = 0.d0
      do 10 i = 1,nfctmx-1
	 ai = i
	 fct(i+1) = fct(i)+dlog(ai)
 10   continue
c
      return
      end

***********************************************************************

      subroutine bndbcele(Eval,fun,potmatrix,xmu
     &      ,rmis,rfin,nvini,nvmax,j,npun,nelec)
      implicit real*8(a-h,o-z)

*
* calculate the bound state of a diatomic system
* for several coupled electronic states (nelec)
*   for a particular value of the angular momentum j
*
* uses a DVR representation of particles in a box
*
*  uses a.u.
* xmu: reduced mass
*
      parameter(npunaux=1024,nelecaux=10,ntotaux=nelecaux*npunaux)

      dimension Eval(nvini:nvmax),fun(npun,nelec,nvini:nvmax)      
      dimension potmatrix(npun,nelec,nelec)

      dimension nr(ntotaux),ne(ntotaux),Hmat(ntotaux,ntotaux)
      dimension eigen(ntotaux),T(ntotaux,ntotaux)
      dimension wwork(5*ntotaux)
      dimension ind(ntotaux*5)

* checking dimensions

      if(npun.gt.npunaux)then
         write(6,*)' npun= ',npun,' > npunaux= ',npunaux
         write(6,*)'  change it in bndbcele '
         stop
      endif

      if(nelec.gt.nelecaux)then
         write(6,*)' nelec= ',nelec,' > nelecaux= ',nelecaux
         write(6,*)'  change it in bndbcele '
         stop
      endif

* forming basis 

      ii=0
      do ie=1,nelec
      do ir=1,npun
         ii=ii+1
         nr(ii)=ir
         ne(ii)=ie
      enddo
      enddo
      ntot=nelec*npun
      ah=(rfin-rmis)/dble(npun-1)

      xl=0.5d0*dble(j*(j+1))/xmu

      pi=dacos(-1.d0)
      xkinfac=pi*0.5d0/ ( (rfin-rmis)*xmu)
      hbr=1.d0

* forming H matrix

      ii=0
      do i=1,ntot
      do ip=i,ntot
         ii=ii+1
         Hmat(i,ip)=0.d0
         ie=ne(i)
         iep=ne(ip)
         ir=nr(i)
         irp=nr(ip)

* rotational term
         if(ir.eq.irp.and.ie.eq.iep)then
            r=rmis+dble(ir-1)*ah

            Hmat(i,ip)= Hmat(i,ip)+ xl/(r*r)
         endif
* radial kinetic term
         if(ie.eq.iep)then
            cint=0.d0
            menos=(ir-irp)
            mas=ir+irp
            sign=dble( (-1)**(menos))
            if(ir.eq.irp)then
               x1=pi*pi/6.d0
               x2=1.d0/( dble(mas)*dble(mas) )
               cint=x1+x2
            else
               x1=1.d0/( dble(menos)*dble(menos) )
               x2=1.d0/( dble(mas)*dble(mas) )
               cint=x1+x2
            endif
            cint=sign*cint/(ah*ah)

            Tmat=hbr*hbr*cint/xmu
            
            Hmat(i,ip)=Hmat(i,ip)+Tmat
         endif
* He+Hso potential terms

         if(ir.eq.irp)then
           Hmat(i,ip)=Hmat(i,ip)+potmatrix(ir,ie,iep)
         endif

      enddo
      enddo

* diagonalization
 
c       call dspev('v','l',ntot,Hmat,eigen,T,ntotaux,wwork,inf)

        call diagon(hmat,ntot,ntotaux,T,eigen)

* keeping desired eigenstates

      iiv=0
      do iv=nvini,nvmax
         iiv=iiv+1
         Eval(iv)=eigen(iiv)
         do ii=1,ntot
            ir=nr(ii)
            ie=ne(ii)
            fun(ir,ie,iv)=T(ii,iiv)/dsqrt(ah)
         enddo
      enddo

      return
      end

********************************* l2mat ***************************************

      subroutine l2mat(bfl2mat,facmass,iommin,iommax,j,Jtot
     &                                   ,iomdim0,iomdim1)
      implicit real*8(a-h,o-z)
      dimension bfl2mat(iomdim0:iomdim1,iomdim0:iomdim1)

*   l^2 matrix obtained in a body-fixed representation
*                  for fixed Jtot and j 
*            and a limited number of Omega projection's
*                   using parity adapted functions
*                           facmass = hbr*hbr*0.5d0/xmred
      
      do iom=iommin,iommax
      do jom=iommin,iommax
         bfl2mat(iom,jom)=0.d0
      enddo
      enddo

      xj=dble(j)

      if(iommin.lt.iomdim0.or.iommax.gt.iomdim1)then
         if(idproc.eq.0)write(6,*)'   problem with Omega:',
     &         ' j,iommin,iommax= '  ,j,iommin,iommax,'  in l2mat'
         stop
      endif

      do iom=iommin,iommax
      do jom=iommin,iommax
         if(iom.eq.jom)then
            x1=xj*(xj+1.d0)
            x2=dble(Jtot*(Jtot+1)-2*iom*jom)
            xk=(x1+x2)*facmass
            bfl2mat(iom,jom)=xk
         elseif(iabs(iom-jom).eq.1)then
            x1=dble(Jtot*(Jtot+1)-iom*jom)
            x2=xj*(xj+1.d0)-dble(iom*jom)
            x12=dsqrt(x1*x2)
            xk=-x12*facmass
            if(iom*jom.eq.0)xk=xk*dsqrt(2.d0)
            bfl2mat(iom,jom)=xk
         endif
      enddo
      enddo
      
      return
      end 

********************** rgauscolini ***********************************
      subroutine rgauscolini(rgaus,r,rcol,alpha,xk,factor,l)
      real*8 rgaus,r,rcol,alpha,xk,factor,pepe,f,g,df,dg,arg
      integer l,key

      pepe=dble(l*(l+1))
      arg=(r-rcol)*xk
c      key=-1
c      CALL BESPH2(F,DF,G,DG,PEPE,ARG,KEY,0)
c      g=-g*dsqrt(2.d0)
      g=dcos(arg)*dsqrt(2.d0)

      rgaus=g*factor*dexp(-(r-rcol)*(r-rcol)/alpha/alpha)

      return
      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine diagon(Hmat,n,ndim,T,eigen)
!     diagonalize Hmat  --> providing eigen (eigenvalues) and T (eigenvectors)
!            which are ordered in increasing energy
      implicit none
      integer n,ndim,nrot
      real*8 :: Hmat(ndim,ndim),eigen(ndim),T(ndim,ndim)

      call jacobi(Hmat,n,ndim,eigen,T,nrot)

      call eigsrt(eigen,T,n,ndim)

      return
      end subroutine diagon
      
!***************************************************
      subroutine jacobi(a,n,np,d,v,nrot)
      implicit none
      integer,parameter :: NMAX=50000
      INTEGER :: n,np,nrot
      REAL*8 :: a(np,np),d(np),v(np,np)
      INTEGER :: i,ip,iq,j
      REAL*8 :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*dabs(a(ip,iq))
            if((i.gt.4).and.(dabs(d(ip))+
     *g.eq.dabs(d(ip))).and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
              a(ip,iq)=0.d0
            else if(dabs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(dabs(h)+g.eq.dabs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/dsqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      stop 'too many iterations in jacdiagon'
      return
      END subroutine jacobi
C  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
!*********************************************************************
      subroutine eigsrt(d,v,n,np)
!ordering eigenvalues and eigenvectors
      implicit none
      INTEGER:: n,np
      REAL*8 :: d(np),v(np,np)
      INTEGER :: i,j,k
      REAL*8 :: p

      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
C  (C) Copr. 1986-92 Numerical Recipes Software *1n#!-013.
      END subroutine eigsrt


