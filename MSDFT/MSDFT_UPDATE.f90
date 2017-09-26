!C*MODULE BLWCAS  *DECK UPDATE
      SUBROUTINE UPDATE(Norb,ENG,T,U,TRANS)
      INTEGER I,J,K,L,M1,M2,M3,M4
      REAL*8 TRANS(Norb,Norb),ENG
      REAL*8 T(Norb,Norb)
      REAL*8 U(Norb,Norb,Norb,Norb)
      REAL*8 MO(Norb,Norb)
      REAL*8 TM(Norb,Norb)
      REAL*8 TD
      INTEGER ITMP(Norb)
!C READ MOs
      OPEN(23,FILE='INT_MO')
      DO I=1,Norb
        READ(23,*)(MO(J,I),J=1,Norb)
      END DO
      CLOSE(23)
      MO=MATMUL(MO,TRANS)
      CALL transform_RAM(NORB,TRANS,T,U)

!C WRITE NEW MOs
      OPEN(23,FILE='INT_MO')
      DO I=1,Norb
        WRITE(23,*)(MO(J,I),J=1,Norb)
      END DO
      CLOSE(23)

      END

