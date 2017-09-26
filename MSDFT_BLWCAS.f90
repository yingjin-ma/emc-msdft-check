!C*MODULE BLWCAS  *DECK BLWCAS
      SUBROUTINE BLWCAS(iPARA,DIM,GEOM,ATOMCHG,iERRO)
      INTEGER I,J,K,L,ICIR,DIM,iERRO
      INTEGER iPARA(DIM),NSTATE
      INTEGER IMAX,IINT,IDFT,iTOLER,iROOT
      INTEGER NORB,NELE,NACT,NDET,iMETHOD,NATOM
      INTEGER ATOMCHG(iPARA(7))
      REAL*8 GEOM(iPARA(7),3)
      REAL*8 ENG,DSUM,TD,Converged
      REAL*8,ALLOCATABLE :: READCI(:,:),RDM1(:,:),RDM2(:,:,:,:)
      REAL*8,ALLOCATABLE :: WEIGHTS(:)
      REAL*8,ALLOCATABLE :: T(:,:),U(:,:,:,:)
      REAL*8,ALLOCATABLE :: TRANSU(:,:)
      REAL*8,ALLOCATABLE :: DFT_Fc(:,:)

!C Initialize of parameters
!C##################################################################
!C need to improve 
!C In this version the closed or frozen orbitals are not avalibale
!C So the active orbitals include all the occupied orbitals
!C This means that NACT=FROZEN+CLOSED+ACTIVE
!C And DETs==RDMDETs
      NACT=iPARA(3)+iPARA(4)+iPARA(5)
      CALL SYSTEM('cp DETs RDMDETs')
!C##################################################################
      IMAX=iPARA(61)
      NSTATE=iPARA(31)
      iROOT=iPARA(33)
      IINT=iPARA(91)
      NELE=iPARA(11)
      NORB=iPARA(1)
      NATOM=iPARA(7)
      NDET=iPARA(21)
      IDFT=iPARA(92)
      iMETHOD=iPARA(95)
      iERRO=0
      iTOLER=iPARA(96)
      IF(iTOLER<=0) THEN
        iTOLER=15
      ENDIF
      IF(iROOT==0) THEN
        iROOT=NSTATE
      ENDIF
      Converged=1.0D0
      WRITE(8406,*)''
      WRITE(8406,*)' > STEP 2.1: Setting Initial Parameters'
      WRITE(8406,*)'IN:  PUT_INT.f90'
      CALL PUT_INT(iPARA,DIM)
      WRITE(8406,*)'OUT: PUT_INT.f90'
      WRITE(8406,*)''
      ALLOCATE(T(NORB,NORB))
      ALLOCATE(U(NORB,NORB,NORB,NORB))
!C Initialization: RDMDETs, DETs, INT_MO, FCIDUMP, ONETWOINT
      IF(IINT==0) THEN
        CALL CHECKFILE(0,iERRO)
        IF(iERRO==-1) GOTO 999
        WRITE(8406,*)''
        WRITE(8406,*)'NO INITIAL GUESS ARE FOUNDED'
        WRITE(8406,*)'IN:  MSDFT_INT.f90'
        CALL CAS_INT(NORB,NATOM,GEOM,ATOMCHG,T,U,ENG)
        WRITE(8406,*)'OUT: MSDFT_INT.f90'
        WRITE(8406,*)''
      ELSE
        CALL FOUND_INIT(NORB,NATOM,GEOM,ATOMCHG)
        CALL CHECKFILE(1,iERRO)
        IF(iERRO==-1) GOTO 999
        WRITE(8406,*)''
        WRITE(8406,*)'INITIAL GUESS FOUNDED'
        WRITE(8406,*)''
      END IF
      ALLOCATE(DFT_Fc(NORB,NORB))
      ALLOCATE(TRANSU(NORB,NORB))
      ALLOCATE(WEIGHTS(NSTATE))
      ALLOCATE(READCI(NDET,NSTATE))
      ALLOCATE(RDM1(NACT,NACT))
      ALLOCATE(RDM2(NACT,NACT,NACT,NACT))
!C##################################################################
!C need to improve 
!C In this version the weights for states are not avalible
!C Equal weights are used
      WEIGHTS=1.0D0
      IF(iPARA(32)==1) THEN
        OPEN(51,FILE='WEIGHTS')
        READ(51,*)(WEIGHTS(I),I=1,NSTATE)
        CLOSE(51)
      ENDIF
      DSUM=0.0D0
      DO I=1,NSTATE
        DSUM=DSUM+WEIGHTS(I)
      END DO
      DO I=1,NSTATE
        WEIGHTS(I)=WEIGHTS(I)/DSUM
      END DO
!C##################################################################
      IF(IINT==1) THEN
      OPEN(25,FILE='ONETWOINT',form='unformatted')
      DO I=1,Norb
      DO J=1,I
      DO K=1,Norb
      DO L=1,K
        READ(25)TD
        U(I,J,K,L)=TD
        U(J,I,K,L)=TD
        U(I,J,L,K)=TD
        U(J,I,L,K)=TD
      END DO
      END DO
      END DO
      END DO
      DO I=1,Norb
      DO J=1,Norb
        READ(25)TD
        T(i,j)=TD
      END DO
      END DO
      READ(25)ENG
      CLOSE(25)
      ENDIF
      WRITE(8406,*)''
      WRITE(8406,*)' > STEP 2.2: CAS(MSDFT)SCF Program'
      WRITE(8406,*)''
      CALL INPUT()
      CALL t_deallocate()
      DO ICIR=1,IMAX
        WRITE(8406,'(A25,I5,A20)')'-----  Iter.',ICIR,'-----          '
        WRITE(8406,*)''
        WRITE(8406,*)'>> STEP 2.2.1 Entering the MSDFT CI calculation '
        WRITE(8406,*)''
        DFT_Fc=0.0D0
        RDM1=0.0D0
        RDM2=0.0D0
        WRITE(8406,*)'Compute 1-RDM, 2-RDM and CI energy '
        WRITE(8406,*)''
        WRITE(8406,*)'IN:  MSDFT_1-2RDM.f90'
        CALL OneAndTwoRDM(iMETHOD,NATOM,NORB,NACT,NDET,RDM1(:,:), &
        RDM2(:,:,:,:),T,U,ENG,IDFT,NSTATE,WEIGHTS(:),iTOLER,iROOT, &
        GEOM,ATOMCHG)
        WRITE(8406,*)'OUT: MSDFT_1-2RDM.f90'
        WRITE(8406,*)''
        WRITE(8406,*)' ------ Normal Termination with STEP 2.2.1 ------'
        WRITE(8406,*)''
        OPEN(24,file='oneRDM.0.0')
        WRITE(24,*)NELE
        DO I=1,Nact
        DO J=1,Nact
          WRITE(24,*)I-1,J-1,RDM1(I,J)
        END DO
        END DO
        CLOSE(24)
        OPEN(24,file='twoRDM.0.0')
        WRITE(24,*)NELE
        DO I=1,Nact
        DO J=1,Nact
          DO K=1,Nact
          DO L=1,Nact
            WRITE(24,*)I-1,K-1,L-1,J-1,RDM2(I,J,K,L)
          END DO
          END DO
        END DO
        END DO
        CLOSE(24)
        WRITE(*,*)' ------ Normal Termination with STEP 2.2.1 ------'
        IF(iPARA(94)==1) GOTO 900
!C STEP 4: 
!C IN:  oneRDM, twoRDM, FCIDUMP, ONETWOINT
!C OUT: TRANSU
        WRITE(8406,*)''
        WRITE(8406,*)'>> STEP 2.2.2 Entering the orbital optimization '
        WRITE(8406,*)''
        WRITE(8406,*)'IN:  mcscf'
!********************************************************************
!Read DFT_Fc from the disk: FCIJ.tmp
        IF(IDFT==1) THEN
          OPEN(35,FILE='FCIJ.tmp')
          READ(35,*)DFT_Fc
          CLOSE(35)
        ENDIF
!********************************************************************
        CALL SCF_ORBOPT(T,U,NORB,ENG,TRANSU,DFT_Fc)
        WRITE(8406,*)'OUT: mcscf'
        WRITE(8406,*)''
        OPEN(41,FILE='SACASENG.tmp')
        READ(41,*)Converged
        CLOSE(41)
        WRITE(*,*)'Converged = ',Converged
        WRITE(8406,*)'orbital Convergence = ',Converged
        IF(DABS(Converged)<=1.0D-5) GOTO 900

!C STEP 5: 
!C IN:  TRANSU, INT_MO, ONETWOINT
!C OUT: INT_MO, ONETWOINT
        IF(IMAX>1) THEN
          WRITE(8406,*)'IN:  MSDFT_UPDATE.f90'
          CALL UPDATE(NORB,ENG,T,U,TRANSU)
          WRITE(8406,*)'OUT: MSDFT_UPDATE.f90'
        ENDIF
        WRITE(8406,*)''
        WRITE(8406,*)' ------ Normal Termination with STEP 2.2.2 ------'
        WRITE(8406,*)''
      END DO
900   CONTINUE
!C WRITE ONETWOINT
      OPEN(24,FILE='ONETWOINT',form='unformatted')
      DO I=1,Norb
      DO J=1,I
      DO K=1,Norb
      DO L=1,K
        TD=U(I,J,K,L)
        WRITE(24)TD
      END DO
      END DO
      END DO
      END DO
      DO I=1,Norb
      DO J=1,Norb
        TD=T(I,J)
        WRITE(24)TD
      END DO
      END DO
      WRITE(24)ENG
      CLOSE(24)

      DEALLOCATE(WEIGHTS)
      DEALLOCATE(READCI)
      DEALLOCATE(RDM1)
      DEALLOCATE(RDM2)
      DEALLOCATE(T)
      DEALLOCATE(U)
      DEALLOCATE(TRANSU,DFT_Fc)
999   CONTINUE
      END  

