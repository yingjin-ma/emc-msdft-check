      SUBROUTINE MOLPRO_GAMESS(NAO,NATOM,BASELABLE,NORDER,IORDER)
      IMPLICIT NONE
      INTEGER NORDER
      INTEGER NAO,NATOM
      INTEGER NMO,I1,I2,ILINE1,ILINE2
      INTEGER IORDER(NORDER)
      CHARACTER*40 CF1,CF2
      REAL*8 GAMESSMO(NAO,NAO)
      REAL*8 MOLPROMO(NAO,NAO)
      CHARACTER*30 BASELABLE

      CF1='MOLPRO.dat'     
      CF2='GAMESS.dat'
      GAMESSMO=0.0D0
      MOLPROMO=0.0D0
      WRITE(*,'(A20,I10)')'atoms:             ',NATOM
      WRITE(*,'(A20,I10)')'orbial:            ',NAO
      WRITE(*,'(A20,A20)')'basis set:         ',BASELABLE
      WRITE(*,'(A20,A12,A12)')'input  filename:   ',CF1
      WRITE(*,'(A20,A12,A12)')'output filename:   ',CF2

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

      IF(ILINE1==0) GOTO 9999
      IF(ILINE1>0 .AND. ILINE2>0) GOTO 9999
      IF(ILINE1>0 .AND. ILINE2==0) THEN
        CALL READMOLPROMO(NAO,MOLPROMO,ILINE1)
        CALL MOLPROTOGAMESS(NAO,MOLPROMO,GAMESSMO,BASELABLE,NATOM)
      ENDIF
      CALL WRITEGAMESSMO(NAO,GAMESSMO)

9999  CONTINUE
      WRITE(*,*)'MOLPRO FILE:',ILINE1
      WRITE(*,*)'GAMESS FILE:',ILINE2
      END

      SUBROUTINE MOLPROTOGAMESS(N,MO1,MO2,BASELABLE,NATOM)
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
        IF(CTYPE=='Ca' .or. CTYPE=='Ca') ATOMCHG(I)=20
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
               CALL ORDER_D3(N,MO2,L-6,L)
             ENDIF
             IF(C1 == "f") THEN
               L=L+10
               CALL ORDER_F3(N,MO2,L-10,L)
             ENDIF
           ENDDO
         ENDIF
       ENDDO
120    CONTINUE
301    CLOSE(23)
      ENDDO
print*,L
      END

      SUBROUTINE ORDER_D3(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER ID1(6),ID2(6)
      INTEGER I,J,K
      DO I=1,6
        ID1(I)=I
      ENDDO
      ID2(1)=1
      ID2(2)=2
      ID2(3)=3
      ID2(4)=4
      ID2(5)=5
      ID2(6)=6
      TMO=MO
      DO I=1,6
        J=N1+ID1(I)
        K=N1+ID2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END

      SUBROUTINE ORDER_F3(N,MO,N1,N2)
      IMPLICIT NONE
      INTEGER N,N1,N2
      REAL*8 MO(N,N),TMO(N,N)
      INTEGER IF1(10),IF2(10)
      INTEGER I,J,K
      DO I=1,10
        IF1(I)=I
      ENDDO
      IF2(1)=1
      IF2(2)=2
      IF2(3)=3
      IF2(4)=5
      IF2(5)=6
      IF2(6)=4
      IF2(7)=9
      IF2(8)=7
      IF2(9)=8
      IF2(10)=10
      TMO=MO
      DO I=1,10
        J=N1+IF1(I)
        K=N1+IF2(I)
        TMO(J,:)=MO(K,:)
      ENDDO
      MO=TMO
      END


      SUBROUTINE READMOLPROMO(N,MO,ILINE)
      IMPLICIT NONE
      INTEGER N,ILINE
      REAL*8 MO(N,N)
      INTEGER I,J,K,IC
      CHARACTER*40 CC
      OPEN(41,FILE='MOLPRO.dat')
      DO K=1,ILINE
        READ(41,*)CC
        IF(CC=='[MO]') THEN
          DO J = 1,N
            DO I=1,4
              READ(41,*)
            ENDDO
            DO I=1,N
              READ(41,*)IC,MO(I,J)
            ENDDO
          ENDDO
          GOTO 130
        ENDIF
      ENDDO
130   CONTINUE
      CLOSE(41)
      END

