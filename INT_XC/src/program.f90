      Program main
      implicit none
      Integer,Parameter :: Natom=6
      Integer,Parameter :: NORB=86
      Integer,Parameter :: NROOTS=2
      Integer iROOT,iSTATE,NOCC
      Integer ON_KAPPA,ON_ZETA
      Real*8 WEIGHTS(2)
      Real*8 GEOM(Natom,3)
      Integer ATOMCHG(Natom),iERRO
      Real*8 T(NORB,NORB)
      Real*8 U(NORB,NORB,NORB,NORB)
      Real*8 Pa(NORB,NORB,NROOTS),Pb(NORB,NORB,NROOTS)
      Real*8 Energy_XC(NROOTS)
      Real*8 COEFF(NROOTS,NROOTS)
      character     :: baselable*30

      ON_KAPPA=1
      ON_ZETA=1
      baselable='aug-cc-pVDZ'
      COEFF=0.0D0
      COEFF(1,1)=DSQRT(1.0D0)
      COEFF(2,1)=DSQRT(0.0D0)
      COEFF(1,2)=DSQRT(0.1D0)
      COEFF(2,2)=-DSQRT(0.9D0)
      iROOT=1
      iSTATE=1
      WEIGHTS(1)=1.0D0
      WEIGHTS(2)=0.0D0
      NOCC=0
      CALL RD_GEOM(NATOM,GEOM,ATOMCHG,iERRO)
      CALL INTXC_ALLOCATE(NATOM)
      CALL BASIS_INT(NATOM,NORB,GEOM,ATOMCHG,T,U,baselable)
      CALL DFT_Fc(NATOM,NORB,NROOTS,Pa,Pb,Energy_XC, &
                  GEOM,ATOMCHG,1,COEFF,iROOT,iSTATE,WEIGHTS,NOCC, &
                  ON_KAPPA,ON_ZETA)
      CALL INTXC_DEALLOCATE()
      END

      SUBROUTINE RD_GEOM(NATOM,GEOM,ATOMCHG,iERRO)
      IMPLICIT NONE
      INTEGER NATOM,iERRO,I
      INTEGER ATOMCHG(NATOM)
      REAL*8  GEOM(NATOM,3)
      CHARACTER*10 CTYPE(NATOM)
      GEOM=0.0D0
      ATOMCHG=0
      OPEN(23,FILE='GEOM.xyz')
      DO I=1,NATOM
        READ(23,*),CTYPE(I),GEOM(I,1),GEOM(I,2),GEOM(I,3)
      ENDDO
      CLOSE(23)
      DO I=1,NATOM
        IF(CTYPE(I)=='H'  .or. CTYPE(I)=='h')  ATOMCHG(I)=1
        IF(CTYPE(I)=='He' .or. CTYPE(I)=='he') ATOMCHG(I)=2
        IF(CTYPE(I)=='Li' .or. CTYPE(I)=='li') ATOMCHG(I)=3
        IF(CTYPE(I)=='Be' .or. CTYPE(I)=='be') ATOMCHG(I)=4
        IF(CTYPE(I)=='B'  .or. CTYPE(I)=='b')  ATOMCHG(I)=5
        IF(CTYPE(I)=='C'  .or. CTYPE(I)=='c')  ATOMCHG(I)=6
        IF(CTYPE(I)=='N'  .or. CTYPE(I)=='n')  ATOMCHG(I)=7
        IF(CTYPE(I)=='O'  .or. CTYPE(I)=='o')  ATOMCHG(I)=8
        IF(CTYPE(I)=='F'  .or. CTYPE(I)=='f')  ATOMCHG(I)=9
        IF(CTYPE(I)=='Ne' .or. CTYPE(I)=='ne') ATOMCHG(I)=10
      ENDDO
      END SUBROUTINE RD_GEOM

      SUBROUTINE CONF_PURE(COEFF,NDET,iROOT,iSTATE,iOUT)
      INTEGER NDET,iROOT,iSTATE
      INTEGER iOUT(iSTATE)
      INTEGER SS(iROOT),TT(iROOT)
      REAL*8 COEFF(NDET,iROOT),DSUM1,DSUM2
      INTEGER I,J,K1,K2,N1(NDET),N2(NDET),iCOUNT,NUM
      INTEGER iS,iT
      N1=0
      N2=0
      iCOUNT=0
      OPEN(23,FILE='couple_index.tmp')
      DO WHILE(.TRUE.)
        READ(23,*,END=601)I,J
        iCOUNT=iCOUNT+1
        N1(iCOUNT)=I
        N2(iCOUNT)=J
      ENDDO
601   CLOSE(23)
      IF(iCOUNT==0) THEN
        WRITE(*,*)'ERROR: CAN NOT FIND couple_index.tmp FILE'
        GOTO 999
      ENDIF
      NUM=iCOUNT

      iS=0
      iT=0
      DO I=1,iROOT
        DSUM1=0.0D0
        DO J=1,NUM
          K1=N1(J)
          K2=N2(J)
          IF(DABS(COEFF(K2,I))>1.0D-10) THEN
            DSUM1=DSUM1+COEFF(K1,I)/COEFF(K2,I)
          ENDIF
        ENDDO
        IF(DSUM1>0.5) THEN
          iS=iS+1
          SS(iS)=I
        ELSEIF(DSUM1<-0.5) THEN
          iT=iT+1
          TT(iT)=I
        ELSE
          WRITE(*,*)'ERROR IN FIND iROOT'
        ENDIF
      ENDDO

      IF(iS<iSTATE) THEN
        WRITE(*,*)'NOT ENOUGH ROOTS FOR iSTATES'
        GOTO 999
      ENDIF

      iOUT(:)=SS(1:iSTATE)
999   CONTINUE
      END

