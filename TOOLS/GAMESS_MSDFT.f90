      SUBROUTINE MSDFT_GAMESS(NAO,NATOM,BASELABLE,NORDER,IORDER)
      IMPLICIT NONE
      INTEGER NORDER
      INTEGER NAO,NATOM
      INTEGER NMO,I1,I2,ILINE1,ILINE2
      INTEGER IORDER(NORDER)
      CHARACTER*40 CF1,CF2
      REAL*8 GAMESSMO(NAO,NAO)
      REAL*8 MSDFTMO(NAO,NAO)
      CHARACTER*30 BASELABLE

      CF1='GAMESS.dat'     
      CF2='MSDFT.dat'
      GAMESSMO=0.0D0
      MSDFTMO=0.0D0 

      WRITE(*,'(A20,I10)')'atoms:             ',NATOM
      WRITE(*,'(A20,I10)')'orbial:            ',NAO
      WRITE(*,'(A20,A20)')'basis set:         ',BASELABLE
      WRITE(*,'(A20,A12,A12)')'output filename:   ',CF1,CF2

      ILINE1=0
      OPEN(23,FILE=CF1)
      DO WHILE(.TRUE.)
        READ(23,*,END=101)
        ILINE1=ILINE1+1
      ENDDO
101   CLOSE(23)

      ILINE2=0
      OPEN(24,FILE=CF2)
      DO WHILE(.TRUE.)
        READ(24,*,END=102)
        ILINE2=ILINE2+1
      ENDDO
102   CLOSE(24)

      IF(ILINE1==0 .AND. ILINE2==0) GOTO 9999
      IF(ILINE1>0 .AND. ILINE2>0) GOTO 9999
      IF(ILINE1>0 .AND. ILINE2==0) THEN
        write(*,*)"ILINE1 : ",ILINE1       
        CALL READGAMESSMO(NAO,GAMESSMO,ILINE1)
        CALL GAMESSTOMSDFT(NAO,GAMESSMO,MSDFTMO,BASELABLE,NATOM)
      ENDIF
      IF(ILINE2>0 .AND. ILINE1==0) THEN
        CALL READMSDFTMO(NAO,MSDFTMO,ILINE2)
        CALL MSDFTTOGAMESS(NAO,MSDFTMO,GAMESSMO,BASELABLE,NATOM)
      ENDIF

      CALL CHANGE(NAO,MSDFTMO,NORDER,IORDER)
      CALL CHANGE(NAO,GAMESSMO,NORDER,IORDER)
      CALL WRITEGAMESSMO(NAO,GAMESSMO)
      CALL WRITEMSDFTMO(NAO,MSDFTMO)

9999  CONTINUE
      WRITE(*,*)'GAMESS FILE:',ILINE1
      WRITE(*,*)'MSDFT  FILE:',ILINE2
      END

      SUBROUTINE GAMESSTOMSDFT(N,MO1,MO2,BASELABLE,NATOM)
      IMPLICIT NONE
      CHARACTER*30 BASELABLE
      INTEGER N,I,J,K,L
      INTEGER NATOM
      INTEGER ATOMCHG(NATOM)
      REAL*8 MO1(N,N),MO2(N,N)
      CHARACTER*30 CC,FILENAME,CTMP,C1,CTYPE

      FILENAME=trim(BASELABLE)//'.bas'
      MO2=MO1
      OPEN(24,FILE='GEOM.xyz')
      DO I=1,NATOM
        READ(24,*,END=101)CTYPE
        IF(CTYPE=='H'  .or. CTYPE=='h')  ATOMCHG(I)=1
        IF(CTYPE=='He' .or. CTYPE=='he') ATOMCHG(I)=2
        IF(CTYPE=='Li' .or. CTYPE=='li') ATOMCHG(I)=3
        IF(CTYPE=='Be' .or. CTYPE=='be') ATOMCHG(I)=4
        IF(CTYPE=='B'  .or. CTYPE=='b')  ATOMCHG(I)=5
        IF(CTYPE=='C'  .or. CTYPE=='c')  ATOMCHG(I)=6
        IF(CTYPE=='N'  .or. CTYPE=='n')  ATOMCHG(I)=7
        IF(CTYPE=='O'  .or. CTYPE=='o')  ATOMCHG(I)=8
        IF(CTYPE=='F'  .or. CTYPE=='f')  ATOMCHG(I)=9
        IF(CTYPE=='Ne' .or. CTYPE=='ne') ATOMCHG(I)=10
        IF(CTYPE=='Na'  .or. CTYPE=='na')  ATOMCHG(I)=11
        IF(CTYPE=='Mg' .or. CTYPE=='mg') ATOMCHG(I)=12
        IF(CTYPE=='Al' .or. CTYPE=='al') ATOMCHG(I)=13
        IF(CTYPE=='Si' .or. CTYPE=='si') ATOMCHG(I)=14
        IF(CTYPE=='P'  .or. CTYPE=='p')  ATOMCHG(I)=15
        IF(CTYPE=='S'  .or. CTYPE=='s')  ATOMCHG(I)=16
        IF(CTYPE=='Cl'  .or. CTYPE=='cl')  ATOMCHG(I)=17
        IF(CTYPE=='Ar'  .or. CTYPE=='ar')  ATOMCHG(I)=18
        IF(CTYPE=='K'  .or. CTYPE=='k')  ATOMCHG(I)=19
        IF(CTYPE=='Ca' .or. CTYPE=='ca') ATOMCHG(I)=20
        IF(CTYPE=='Fe' .or. CTYPE=='fe') ATOMCHG(I)=26
      ENDDO
101   CLOSE(24)

      L=0
      DO I=1,NATOM
       WRITE(CTMP,"(I0.2,A1,A27)")ATOMCHG(I),'-',BASELABLE
       OPEN(23,FILE=FILENAME)
       DO WHILE(.TRUE.)
        READ(23,*,END=301)CC
        IF(CC==CTMP) THEN
           DO WHILE(.TRUE.)
             READ(23,*)CC,C1
             IF(C1 == "END") GOTO 120
             IF(C1 == "s") L=L+1
             IF(C1 == "p") L=L+3
             IF(C1 == "d") THEN
               L=L+6
               CALL ORDER_D1(N,MO2,L-6,L)
             ENDIF
             IF(C1 == "f") THEN
               L=L+10
               CALL ORDER_F1(N,MO2,L-10,L)
             ENDIF
             IF(C1 == 'g') THEN
               L=L+15
               CALL ORDER_G1(N,MO2,L-15,L)
             ENDIF
           ENDDO
         ENDIF
       ENDDO
120    CONTINUE
301    CLOSE(23)
      ENDDO
      END

      SUBROUTINE ORDER_D1(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER ID1(6),ID2(6)
      INTEGER I,J,K
      DO I=1,6
        ID1(I)=I
      ENDDO
      ID2(1)=1
      ID2(2)=4
      ID2(3)=5
      ID2(4)=2
      ID2(5)=6
      ID2(6)=3
      TMO=MO
      DO I=1,6
        J=N1+ID1(I)
        K=N1+ID2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END

      SUBROUTINE ORDER_F1(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER IF1(10),IF2(10)
      INTEGER I,J,K
      DO I=1,10
        IF1(I)=I
      ENDDO
      IF2(1)=1
      IF2(2)=4
      IF2(3)=5
      IF2(4)=6
      IF2(5)=10
      IF2(6)=8
      IF2(7)=2
      IF2(8)=7
      IF2(9)=9
      IF2(10)=3
      TMO=MO
      DO I=1,10
        J=N1+IF1(I)
        K=N1+IF2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END

      SUBROUTINE ORDER_G1(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER ID1(15),ID2(15)
      INTEGER I,J,K
      DO I=1,15
        ID1(I)=I
      ENDDO
      ID2(1)=1
      ID2(2)=4
      ID2(3)=5
      ID2(4)=10
      ID2(5)=13
      ID2(6)=11
      ID2(7)=6
      ID2(8)=14
      ID2(9)=15
      ID2(10)=8
      ID2(11)=2
      ID2(12)=7
      ID2(13)=12
      ID2(14)=9
      ID2(15)=3
      TMO=MO
      DO I=1,15
        J=N1+ID1(I)
        K=N1+ID2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END


      SUBROUTINE MSDFTTOGAMESS(N,MO1,MO2,BASELABLE,NATOM)
      IMPLICIT NONE
      CHARACTER*30 BASELABLE
      INTEGER N,I,J,K,L
      INTEGER NATOM
      INTEGER ATOMCHG(NATOM)
      REAL*8 MO1(N,N),MO2(N,N)
      CHARACTER*30 CC,FILENAME,CTMP,C1,CTYPE

      FILENAME=trim(BASELABLE)//'.bas'
      MO2=MO1
      OPEN(24,FILE='GEOM.xyz')
      DO I=1,NATOM
        READ(24,*,END=101)CTYPE
        IF(CTYPE=='H'  .or. CTYPE=='h')  ATOMCHG(I)=1
        IF(CTYPE=='He' .or. CTYPE=='he') ATOMCHG(I)=2
        IF(CTYPE=='Li' .or. CTYPE=='li') ATOMCHG(I)=3
        IF(CTYPE=='Be' .or. CTYPE=='be') ATOMCHG(I)=4
        IF(CTYPE=='B'  .or. CTYPE=='b')  ATOMCHG(I)=5
        IF(CTYPE=='C'  .or. CTYPE=='c')  ATOMCHG(I)=6
        IF(CTYPE=='N'  .or. CTYPE=='n')  ATOMCHG(I)=7
        IF(CTYPE=='O'  .or. CTYPE=='o')  ATOMCHG(I)=8
        IF(CTYPE=='F'  .or. CTYPE=='f')  ATOMCHG(I)=9
        IF(CTYPE=='Ne' .or. CTYPE=='ne') ATOMCHG(I)=10
        IF(CTYPE=='Na'  .or. CTYPE=='na')  ATOMCHG(I)=11
        IF(CTYPE=='Mg' .or. CTYPE=='mg') ATOMCHG(I)=12
        IF(CTYPE=='Al' .or. CTYPE=='al') ATOMCHG(I)=13
        IF(CTYPE=='Si' .or. CTYPE=='si') ATOMCHG(I)=14
        IF(CTYPE=='P'  .or. CTYPE=='p')  ATOMCHG(I)=15
        IF(CTYPE=='S'  .or. CTYPE=='s')  ATOMCHG(I)=16
        IF(CTYPE=='Cl'  .or. CTYPE=='cl')  ATOMCHG(I)=17
        IF(CTYPE=='Ar'  .or. CTYPE=='ar')  ATOMCHG(I)=18
        IF(CTYPE=='K'  .or. CTYPE=='k')  ATOMCHG(I)=19
        IF(CTYPE=='Ca' .or. CTYPE=='ca') ATOMCHG(I)=20
        IF(CTYPE=='Fe' .or. CTYPE=='fe') ATOMCHG(I)=26
      ENDDO
101   CLOSE(24)

      L=0
      DO I=1,NATOM
       WRITE(CTMP,"(I0.2,A1,A27)")ATOMCHG(I),'-',BASELABLE
       OPEN(23,FILE=FILENAME)
       DO WHILE(.TRUE.)
        READ(23,*,END=301)CC
        IF(CC==CTMP) THEN
           DO WHILE(.TRUE.)
             READ(23,*)CC,C1
             IF(C1 == "END") GOTO 120
             IF(C1 == "s") L=L+1
             IF(C1 == "p") L=L+3
             IF(C1 == "d") THEN
               L=L+6
               CALL ORDER_D2(N,MO2,L-6,L)
             ENDIF
             IF(C1 == "f") THEN
               L=L+10
               CALL ORDER_F2(N,MO2,L-10,L)
             ENDIF
             IF(C1 == "g") THEN
               L=L+15
               CALL ORDER_G2(N,MO2,L-15,L)
             ENDIF
           ENDDO
         ENDIF
       ENDDO
120    CONTINUE
301    CLOSE(23)
      ENDDO
      END

      SUBROUTINE ORDER_D2(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER ID1(6),ID2(6)
      INTEGER I,J,K
      DO I=1,6
        ID1(I)=I
      ENDDO
      ID2(1)=1
      ID2(2)=4
      ID2(3)=6
      ID2(4)=2
      ID2(5)=3
      ID2(6)=5
      TMO=MO
      DO I=1,6
        J=N1+ID1(I)
        K=N1+ID2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END

      SUBROUTINE ORDER_F2(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER IF1(10),IF2(10)
      INTEGER I,J,K
      DO I=1,10
        IF1(I)=I
      ENDDO
      IF2(1)=1
      IF2(2)=7
      IF2(3)=10
      IF2(4)=2
      IF2(5)=3
      IF2(6)=4
      IF2(7)=8
      IF2(8)=6
      IF2(9)=9
      IF2(10)=5
      TMO=MO
      DO I=1,10
        J=N1+IF1(I)
        K=N1+IF2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END

      SUBROUTINE ORDER_G2(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER ID1(15),ID2(15)
      INTEGER I,J,K
      DO I=1,15
        ID1(I)=I
      ENDDO
      ID2(1)=1
      ID2(2)=11
      ID2(3)=15
      ID2(4)=2
      ID2(5)=3
      ID2(6)=7
      ID2(7)=12
      ID2(8)=10
      ID2(9)=14
      ID2(10)=4
      ID2(11)=6
      ID2(12)=13
      ID2(13)=5
      ID2(14)=8
      ID2(15)=9
      TMO=MO
      DO I=1,15
        J=N1+ID1(I)
        K=N1+ID2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END


      SUBROUTINE READGAMESSMO(N,MO,ILINE)
      IMPLICIT NONE
      INTEGER N,ILINE
      REAL*8 MO(N,N)
      INTEGER I,J,K,iMAX,iMIN,IC
      CHARACTER*40 CC
      OPEN(41,FILE='GAMESS.dat')
      DO K=1,ILINE
        READ(41,*)CC
        IF(CC=='$VEC') THEN
          DO 120 J = 1,N
          write(6,*)"MAYJ J :",J 
          iMAX = 0
100       iMIN = iMAX+1
          iMAX = iMAX+5
          IF (iMAX .GT. N) iMAX = N
          READ (41,140) IC,IC,(MO(I,J),I = iMIN,iMAX)
          IF (iMAX .LT. N) GO TO 100
          write(6,*)"MAYJ J done "
120       CONTINUE          
        GOTO 130
        ENDIF
      ENDDO
130   CONTINUE
      CLOSE(41)
140   FORMAT(I2,I3,1P,5E15.8)
      END

      SUBROUTINE WRITEGAMESSMO(N,MO)
      IMPLICIT NONE
      INTEGER N,ILINE
      REAL*8 MO(N,N)
      INTEGER I,J,K,iMAX,iMIN,IC
      CHARACTER*40 CC
      OPEN(41,FILE='GAMESS.dat')
      WRITE(41,*)'$VEC'
      DO J = 1,N
        IC=1
        iMAX = 0
100     iMIN = iMAX+1
        iMAX = iMAX+5
        IF (iMAX .GT. N) iMAX = N
        WRITE (41,140) J,IC,(MO(I,J),I = iMIN,iMAX)
        IF (iMAX .LT. N) THEN
          IC=IC+1
          GO TO 100
        ENDIF
      ENDDO
      WRITE(41,*)'$END'
      CLOSE(41)
140   FORMAT(I2,I3,1P,5E15.8)
      END

      SUBROUTINE READMSDFTMO(N,MO,ILINE)      
      IMPLICIT NONE
      INTEGER N,ILINE,I,J
      REAL*8 MO(N,N)
      OPEN(42,FILE='MSDFT.dat')
      DO I=1,N
        READ(42,*)(MO(J,I),J=1,N)
      ENDDO
      CLOSE(42)
      END

      SUBROUTINE WRITEMSDFTMO(N,MO)
      IMPLICIT NONE
      INTEGER N,I,J
      REAL*8 MO(N,N)
      OPEN(42,FILE='MSDFT.dat')
      DO I=1,N
        WRITE(42,*)(MO(J,I),J=1,N)
      ENDDO
      CLOSE(42)
      END

      SUBROUTINE CHANGE(N,MO,NORDER,IORDER)
      IMPLICIT NONE
      INTEGER N,NORDER,I,J,K,iMIN
      INTEGER IORDER(NORDER),iINDEX(N)
      REAL*8 MO(N,N),TMO(N,N)
      DO I=1,N
        iINDEX(I)=I
      ENDDO
      iMIN=IORDER(1)
      DO I=2,NORDER
        J=IORDER(I)
        IF(J<iMIN) iMIN=J
      ENDDO
      iINDEX(iMIN:(iMIN+NORDER-1))=IORDER(1:NORDER)
      DO I=1,N
        J=iINDEX(I)
        TMO(:,I)=MO(:,J)
      ENDDO
      MO=TMO
      END

