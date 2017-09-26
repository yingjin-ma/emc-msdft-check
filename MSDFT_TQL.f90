      SUBROUTINE TQL(MD,N,Z,D,E)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION   D(MD),E(MD),Z(MD,MD)
      EPS=1.0D-15
      NITER=50
      CALL TRED2(MD,N,Z,D,E)
      DO 10 I=2,N
  10  E(I-1)=E(I)
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
      DO 20 L=1,N
      J=0
      H=EPS*(DABS(D(L))+DABS(E(L)))
      LP1=L+1
      IF (B-H) 30,40,40
  30  B=H
  40  DO 50 M=L,N
      IF (DABS(E(M))-B) 60,60,50
  50  CONTINUE
  60  IF (M-L) 70,80,70
  70  IF (J-NITER) 90,100,90
  90  J=J+1
      P=(D(LP1)-D(L))/(2*E(L))
      R=Dsqrt(P*P+1)
      IF (P) 110,111,111
  110 H=D(L)-E(L)/(P-R)
      GOTO 130
  111 H=D(L)-E(L)/(P+R)
  130 DO 140 I=L,N
  140 D(I)=D(I)-H
      F=F+H
      P=D(M)
      C=1.0D0
      S=0.0D0
      MM1=M-1
      IF (MM1-L) 270,280,280
  280 DO 120 LMIP=L,MM1
      I=L+MM1-LMIP
      IP1=I+1
      G=C*E(I)
      H=C*P
      IF (DABS(P)-DABS(E(I))) 160,170,170
  170 C=E(I)/P
      R=Dsqrt(C*C+1.0D0)
      E(IP1)=S*P*R
      S=C/R
      C=1.0D0/R
      GOTO 180
  160 C=P/E(I)
      R=Dsqrt(C*C+1)
      E(IP1)=S*E(I)*R
      S=1/R
      C=C/R
  180 P=C*D(I)-S*G
      D(IP1)=H+S*(C*G+S*D(I))
      DO 190 K=1,N
      H=Z(K,IP1)
      Z(K,IP1)=S*Z(K,I)+C*H
  190 Z(K,I)=C*Z(K,I)-S*H
  120 CONTINUE
  270 E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L))-B) 80,80,70
  80  D(L)=D(L)+F
  20  CONTINUE
      DO 112 I=1,N
      IP1=I+1
      K=I
      P=D(I)
      IF (N-I) 230,230,300
  300 DO 210 J=IP1,N
      IF (D(J)-P) 220,210,210
  220 K=J
      P=D(J)
  210 CONTINUE
  230 IF (K-I) 240,112,240
  240 D(K)=D(I)
      D(I)=P
      DO 260 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
  260 Z(J,K)=P
  112 CONTINUE
      RETURN
  100 STOP '  FAIL'
      END
      subroutine EIGTQL(mdimx, ndimx, dzxmat, deigen, dworkx)
      implicit double precision (a-h,o-z)
      parameter(dtorle = 1.0d-12, mitera = 50)
      dimension dzxmat(mdimx, mdimx), deigen( mdimx), dworkx( mdimx)
      call TRED2(mdimx, ndimx, dzxmat, deigen, dworkx)
      do 10 jdimi=2,ndimx
  10  dworkx(jdimi-1)=dworkx(jdimi)
      dworkx(ndimx)=0.0D0
      dfactf=0.0D0
      dfactb=0.0D0
      do 20 jdiml=1,ndimx
         nitera=0
         dfacth=dtorle*(dabs(deigen(jdiml))+dabs(dworkx(jdiml)))
         jdiml1=jdiml+1
         if (dfactb-dfacth) 30,40,40
  30     dfactb=dfacth
  40     do 50 jdimm=jdiml,ndimx
            if (dabs(dworkx(jdimm))-dfactb) 60,60,50
  50     continue
  60     if (jdimm-jdiml) 70,80,70
  70     if (nitera-mitera) 90,100,90
  90     nitera = nitera + 1
         dfactp=(deigen(jdiml1)-deigen(jdiml))/(2*dworkx(jdiml))
         dfactr=dsqrt(dfactp*dfactp+1.0d+0)
         if (dfactp) 110,111,111
  110    dfacth = deigen(jdiml)-dworkx(jdiml)/(dfactp-dfactr)
         goto 130
  111    dfacth = deigen(jdiml)-dworkx(jdiml)/(dfactp+dfactr)
  130    do 140 jdimi = jdiml, ndimx
  140    deigen(jdimi)=deigen(jdimi)-dfacth
         dfactf = dfactf + dfacth
         dfactp =deigen(jdimm)
         dfactc =1.0D0
         dfacts =0.0D0
         jdimm1 = jdimm - 1
         if (jdimm1-jdiml) 270,280,280
  280    do 120 jdimx=jdiml, jdimm1
            jdimi = jdiml + jdimm1 - jdimx
            jdimi1 = jdimi + 1
            dfactg = dfactc * dworkx(jdimi)
            dfacth = dfactc * dfactp
            if (dabs(dfactp)-dabs(dworkx(jdimi))) 160,170,170
  170       dfactc = dworkx(jdimi)/dfactp
            dfactr = dsqrt(dfactc*dfactc + 1.0D0)
            dworkx(jdimi1)=dfacts*dfactp*dfactr
            dfacts = dfactc / dfactr
            dfactc = 1.0d+0 / dfactr
            goto 180
  160       dfactc = dfactp/dworkx(jdimi)
            dfactr = dsqrt(dfactc*dfactc+1.0d+0)
            dworkx(jdimi1)=dfacts*dworkx(jdimi)*dfactr
            dfacts = 1.0d+0 / dfactr
            dfactc = dfactc / dfactr
  180       dfactp = dfactc * deigen(jdimi) - dfacts * dfactg
            deigen(jdimi1)=dfacth+dfacts*(dfactc*dfactg+dfacts*deigen(jdimi))
            do 190 jdimk=1,ndimx
               dfacth=dzxmat(jdimk, jdimi1)
               dzxmat(jdimk, jdimi1)=dfacts*dzxmat(jdimk,jdimi)+dfactc*dfacth
               dzxmat(jdimk, jdimi )=dfactc*dzxmat(jdimk,jdimi)-dfacts*dfacth
  190       continue
  120    continue
  270    dworkx(jdiml)=dfacts*dfactp
         deigen(jdiml)=dfactc*dfactp
         if (dabs(dworkx(jdiml))-dfactb) 80,80,70
   80    deigen(jdiml)=deigen(jdiml)+dfactf
   20 continue
      do 112 jdimi=1,ndimx
         jdimi1=jdimi+1
         jdimk =jdimi
         dfactp=deigen(jdimi)
         if (ndimx-jdimi) 230,230,300
  300    do 210 jdimj=jdimi1,ndimx
            if (deigen(jdimj)-dfactp) 220,210,210
  220    jdimk = jdimj
         dfactp=deigen(jdimj)
  210    continue
  230    if (jdimk - jdimi) 240,112,240
  240    deigen(jdimk)=deigen(jdimi)
         deigen(jdimi)=dfactp
         do 260 jdimj=1,ndimx
            dfactp=dzxmat(jdimj, jdimi)
            dzxmat(jdimj, jdimi)=dzxmat(jdimj, jdimk)
            dzxmat(jdimj, jdimk)=dfactp
  260    continue
  112 continue
      return
  100 write(*,9001)
 9001 format(' fail in EIGTQL ')
      stop
      end
      subroutine TRED2(mdimx, ndimx, dzxmat, deigen, dworkx)
      implicit double precision (a-h,o-z)
      parameter( dtorle = 1.0d-30)
       dimension dzxmat(mdimx,mdimx), deigen(mdimx), dworkx(mdimx)
      do 20 jdimx=2,ndimx
         jdimi=ndimx+2-jdimx
         jdimi1=jdimi-1
         jdimi2=jdimi-2
         jdiml =jdimi2
         dfactf = dzxmat(jdimi, jdimi1)
         dfactg =0.0D0
         if (jdiml) 30,30,40
  40     do 50 jdimk = 1, jdiml
  50     dfactg = dfactg + dzxmat(jdimi, jdimk)*dzxmat(jdimi, jdimk)
  30     dfacth = dfactg + dfactf * dfactf
         if (dfactg-dtorle) 60,60,70
  60     dworkx(jdimi)=dfactf
         dfacth=0.0D0
         goto 180
  70     jdiml = jdiml + 1
         if (dfactf) 80,90,90
  90     dworkx(jdimi)=-dsqrt(dfacth)
         dfactg=dworkx(jdimi)
         goto 100
  80     dworkx(jdimi)=dsqrt(dfacth)
         dfactg = dworkx(jdimi)
 100     dfacth = dfacth - dfactf*dfactg
         dzxmat(jdimi, jdimi1)=dfactf-dfactg
         dfactf=0.0D0
         do 110 jdimj = 1, jdiml
            dzxmat(jdimj, jdimi)=dzxmat(jdimi, jdimj)/dfacth
            dfactg=0.0D0
            do 201 jdimk = 1, jdimj
  201       dfactg = dfactg + dzxmat(jdimj,jdimk)*dzxmat(jdimi,jdimk)
            jdimj1 = jdimj + 1
            if (jdimj1 - jdiml) 130,130,140
  130       do 120 jdimk = jdimj1, jdiml
  120       dfactg = dfactg + dzxmat(jdimk,jdimj)*dzxmat(jdimi,jdimk)
  140       dworkx(jdimj)=dfactg / dfacth
            dfactf = dfactf + dfactg * dzxmat(jdimj, jdimi)
  110    continue
         dhh = dfactf / ( dfacth + dfacth)
         do 160 jdimj = 1, jdiml
            dfactf = dzxmat(jdimi, jdimj)
            dworkx(jdimj) = dworkx(jdimj)- dhh * dfactf
            dfactg = dworkx(jdimj)
            do 170 jdimk = 1, jdimj
  170       dzxmat(jdimj, jdimk)=dzxmat(jdimj,jdimk)-dfactf*dworkx(jdimk)-dfactg*dzxmat(jdimi, jdimk)
  160    continue
  180    deigen(jdimi)=dfacth
  20  continue
      deigen(1)=0.0D0
      dworkx(1)=0.0D0
      do 190 jdimi=1,ndimx
         jdiml=jdimi-1
         if (deigen(jdimi)) 202,210,202
  202    if (jdiml) 210,210,220
  220    do 230 jdimj=1,jdiml
            dfactg=0.0D0
            DO 240 jdimk = 1, jdiml
  240       dfactg=dfactg+dzxmat(jdimi,jdimk)*dzxmat(jdimk, jdimj)
            do 250 jdimk = 1, jdiml
  250       dzxmat(jdimk, jdimj)=dzxmat(jdimk,jdimj)-dfactg*dzxmat(jdimk, jdimi)
  230    continue
  210    deigen(jdimi)=dzxmat(jdimi, jdimi)
         dzxmat(jdimi, jdimi)=1.0D0
         if (jdiml) 260,260,270
  270    do 280 jdimj = 1, jdiml
            dzxmat(jdimi, jdimj)=0.0D0
  280       dzxmat(jdimj, jdimi)=0.0D0
  260    continue
  190 continue
      return
      end

