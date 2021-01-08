      MODULE GLOBLE_COMM
      IMPLICIT NONE
      INTEGER COMM_NATOM
      INTEGER COMM_IMULT
      INTEGER COMM_ICHARGE
      INTEGER COMM_KAPPA
      INTEGER COMM_ZETA
      INTEGER COMM_NORB
      INTEGER,ALLOCATABLE :: COMM_ATOMCHG(:)
      REAL*8,ALLOCATABLE :: COMM_GEOM(:,:)
      CHARACTER COMM_FUNCTIONAL*30
      END MODULE GLOBLE_COMM

      SUBROUTINE SET_COMM(NATOM,SPIN,ATOMCHG,GEOM,ON_KAPPA,ON_ZETA,NORB,CHARGE)
      USE GLOBLE_COMM
      IMPLICIT NONE
      INTEGER NATOM,SPIN,ON_KAPPA,ON_ZETA,NORB,CHARGE
      INTEGER ATOMCHG(NATOM)
      REAL*8 GEOM(NATOM,3)
      COMM_NATOM=NATOM
      COMM_IMULT=SPIN
      COMM_ICHARGE=CHARGE
      COMM_KAPPA=ON_KAPPA
      COMM_ZETA=ON_ZETA
      COMM_NORB=NORB
      ALLOCATE(COMM_ATOMCHG(COMM_NATOM))
      ALLOCATE(COMM_GEOM(COMM_NATOM,3))
      COMM_ATOMCHG=ATOMCHG
      COMM_GEOM=GEOM      
      END

      SUBROUTINE CLOSE_COMM()
      USE GLOBLE_COMM
      IMPLICIT NONE
      DEALLOCATE(COMM_ATOMCHG)
      DEALLOCATE(COMM_GEOM)
      END

      SUBROUTINE METHOD100(NORB,NDETS,XHF,EXC)
      IMPLICIT NONE
      INTEGER NORB,NA,NB,NDETS  ! input
      REAL*8 XHF(NDETS,NDETS) ! input
      REAL*8 EXC(NDETS,NDETS)  ! input and output
      !worktmp: integertmp, realtmp and charactertmp
      INTEGER DET(NORB,2,NDETS)
      INTEGER INDEX1(NORB,2)
      INTEGER INDEX2(NORB,2)
      INTEGER T1(NORB,2),T2(NORB,2)
      INTEGER I,J,K,L,N(2),ICOUNT
      REAL*8 COEFF(NORB,NORB)
      REAL*8 TEXC1,TEXC2
      !******************************************
      ! READ DETs and INT_MO from disk
      DET=0
      OPEN(23,file='DETs')
      READ(23,*)NA,NB
      DO I=1,NDETS
        READ(23,*)(DET(J,1,I),J=1,NA)
        READ(23,*)(DET(J,2,I),J=1,NB)
      ENDDO
      CLOSE(23)
      COEFF=0.0D0
      OPEN(23,FILE='INT_MO')
      DO I=1,NORB
        READ(23,*)(COEFF(J,I),J=1,NORB)
      ENDDO
      CLOSE(23)
      !*****************************************
      DO I=1,NDETS
        INDEX1=0
        DO K=1,NA
          INDEX1(DET(K,1,I),1)=1
        ENDDO
        DO K=1,NB
          INDEX1(DET(K,2,I),2)=1
        ENDDO
        DO J=I+1,NDETS
          INDEX2=0
          DO K=1,NA
            INDEX2(DET(K,1,J),1)=1
          ENDDO
          DO K=1,NB
            INDEX2(DET(K,2,J),2)=1
          ENDDO
          L=COUNT(ABS(INDEX1(:,1)-INDEX2(:,1)) .EQ. 1)
          IF((COUNT((INDEX1(:,1)-INDEX1(:,2)) .EQ. 0)==NORB).AND. &
             (COUNT((INDEX2(:,1)-INDEX2(:,2)) .EQ. 0)==NORB).AND. &
             (L==2)) THEN
            ICOUNT=0
            DO K=1,NORB
              IF(ABS(INDEX1(K,1)-INDEX2(K,1))==1) THEN
                ICOUNT=ICOUNT+1
                N(ICOUNT)=K
                IF(ICOUNT==2) GOTO 100
              ENDIF
            ENDDO
100         CONTINUE
            T1(:,:)=INDEX2(:,:)
            T2(:,:)=INDEX2(:,:)
            T1(N(1),1)=1
            T1(N(2),1)=0
            T1(N(1),2)=0
            T1(N(2),2)=1
            T2(N(1),1)=1
            T2(N(2),1)=1
            T2(N(1),2)=0
            T2(N(2),2)=0
            CALL FINDEDET(NORB,NDETS,NA,NB,T1,DET,K)
            CALL XC100(T1,T2,COEFF,NORB,EXC(I,J),TEXC1,TEXC2)
! General Case:
!            EXC(I,J)=EXC(I,J)*((EXC(I,I)+EXC(J,J))/(XHF(I,I)+XHF(J,J)))*(XHF(K,K)/EXC(K,K))
! Diatom Case:
!            EXC(I,J)=(EXC(I,I)+EXC(J,J))/2.0D0+TEXC2-2.0D0*TEXC1
            CALL TDF34(XHF,EXC,DET,NA,NB,NDETS,NORB,I,J,T1,TEXC1,TEXC2)
            EXC(I,J)=TEXC2
          ELSEIF((COUNT((INDEX1(:,1)-INDEX2(:,2)) .EQ. 0)==NORB).AND. &
             (COUNT((INDEX2(:,1)-INDEX1(:,2)) .EQ. 0)==NORB).AND. &
             (L==2)) THEN
            ICOUNT=0
            DO K=1,NORB
              IF(ABS(INDEX1(K,1)-INDEX2(K,1))==1) THEN
                ICOUNT=ICOUNT+1
                N(ICOUNT)=K
                IF(ICOUNT==2) GOTO 200
              ENDIF
            ENDDO
200         CONTINUE
            T1(:,:)=INDEX1(:,:)
            T2(:,:)=INDEX1(:,:)
            T2(N(1),1)=1
            T2(N(2),1)=1
            T2(N(1),2)=0
            T2(N(2),2)=0
            CALL XC100(T1,T2,COEFF,NORB,EXC(I,J),TEXC1,TEXC2)
            EXC(I,J)=EXC(I,J)
          ELSE
            EXC(I,J)=XHF(I,J)*(EXC(I,I)+EXC(J,J))/(XHF(I,I)+XHF(J,J))
          ENDIF
          EXC(J,I)=EXC(I,J)
        ENDDO
      ENDDO
      END

      SUBROUTINE XC100(DET1,DET2,COEFF,NORB,EXC,EXCAB,EXCAA)
      USE GLOBLE_COMM
      IMPLICIT NONE
      INTEGER NORB
      INTEGER DET1(NORB,2),DET2(NORB,2)
      REAL*8 COEFF(NORB,NORB)
      REAL*8 EXC
      INTEGER I,J,K,NA,NB
      REAL*8 MOA(NORB,NORB),MOB(NORB,NORB)
      REAL*8 PA(NORB,NORB),PB(NORB,NORB)
      REAL*8 EXCAB,EXCAA
      NA=0
      MOA=0
      DO I=1,NORB
        IF(DET1(I,1)==1) THEN
          NA=NA+1
          MOA(:,NA)=COEFF(:,I)
        ENDIF
      ENDDO
      NB=0
      MOB=0
      DO I=1,NORB
        IF(DET1(I,2)==1) THEN
          NB=NB+1
          MOB(:,NB)=COEFF(:,I)
        ENDIF
      ENDDO
      CALL RHOAB(NORB,NORB,NA,NB,MOA,MOB,PA,PB)
      CALL ENGINEUP(COMM_NATOM,1,COMM_ICHARGE, &
                    COMM_FUNCTIONAL,&
                    PA,PB,PA,PB,COMM_NORB,COMM_GEOM, &
                    COMM_ATOMCHG,&
                    EXCAB,1,COMM_KAPPA,COMM_ZETA)
      NA=0
      MOA=0
      DO I=1,NORB
        IF(DET2(I,1)==1) THEN
          NA=NA+1
          MOA(:,NA)=COEFF(:,I)
        ENDIF
      ENDDO
      NB=0
      MOB=0
      DO I=1,NORB
        IF(DET2(I,2)==1) THEN
          NB=NB+1
          MOB(:,NB)=COEFF(:,I)
        ENDIF
      ENDDO
      CALL RHOAB(NORB,NORB,NA,NB,MOA,MOB,PA,PB)
      CALL ENGINEUP(COMM_NATOM,3,COMM_ICHARGE, &
                    COMM_FUNCTIONAL,&
                    PA,PB,PA,PB,COMM_NORB,COMM_GEOM, &
                    COMM_ATOMCHG,&
                    EXCAA,1,COMM_KAPPA,COMM_ZETA)
      EXC=EXCAB-EXCAA
      END

      SUBROUTINE FINDEDET(NORB,NDETS,NA,NB,T1,DET,NUM)
      IMPLICIT NONE
      INTEGER NORB,NDETS,NUM,NA,NB
      INTEGER T1(NORB,2)
      INTEGER DET(NORB,2,NDETS)
      INTEGER I,J,K,L
      INTEGER INDEX1(NORB,2)
      DO I=1,NDETS
        INDEX1=0
        DO K=1,NA
          INDEX1(DET(K,1,I),1)=1
        ENDDO
        DO K=1,NB
          INDEX1(DET(K,2,I),2)=1
        ENDDO
        L=0
        DO J=1,NORB
          DO K=1,2
            L=L+ABS(INDEX1(J,K)-T1(J,K))
          ENDDO
        ENDDO
        IF(L==0) THEN
          NUM=I
          GOTO 300
        ENDIF
300   CONTINUE
      ENDDO
      END

      SUBROUTINE TDF34(XHF,EXC,DET,NA,NB,N,NORB,II,JJ,DETAB,EXCAB,EXCAA)
      IMPLICIT NONE
      INTEGER NA,NB,N,NORB,II,JJ
      REAL*8 XHF(N,N),EXC(N,N)
      INTEGER DET(NORB,2,N)
      INTEGER DETAB(NORB,2),TMPAB(NORB,2)
      REAL*8 EXCAA,EXCAB,EXCIJ
      INTEGER I,J,K,L1,IA,IB,ISUM
      REAL*8 TMP1,TMP2,TMP3,TMP4,TMP5

      IA=0
      IB=0
      DO I=1,NORB
        IF(DETAB(I,1)==1) THEN
          IA=IA+1
          TMPAB(IA,1)=I
        ENDIF  
        IF(DETAB(I,2)==1) THEN
          IB=IB+1
          TMPAB(IB,2)=I
        ENDIF
      ENDDO
      DO I=1,N
        ISUM=0
        DO J=1,NORB
          DO K=1,2
            ISUM=ISUM+ABS(TMPAB(J,K)-DET(J,K,I))
          ENDDO
        ENDDO
        IF(ISUM==0) THEN
          L1=I
          GOTO 100
        ENDIF
      ENDDO
100   CONTINUE

L1=3
print*,II,JJ,L1
print*,EXCAA,EXCAB
!      TMP1=DSQRT((XHF(II,II)-XHF(JJ,JJ))**2+4.0D0*(XHF(II,JJ))**2)
! with dV12=0
      TMP1=2.0D0*XHF(II,JJ)
      TMP2=XHF(II,II)+XHF(JJ,JJ)-2.0D0*XHF(L1,L1)
      TMP3=(TMP1+2.0D0*XHF(II,JJ))/TMP2
      TMP4=XHF(II,II)+XHF(JJ,JJ)+EXC(II,II)+EXC(JJ,JJ)-2.0D0*(XHF(L1,L1)+EXCAB)
!      TMP5=(TMP3*TMP4-2.0D0*(XHF(II,JJ)+EXCAB-EXCAA))**2-(XHF(II,II)+EXC(II,II)-XHF(JJ,JJ)-EXC(JJ,JJ))**2
      TMP5=(TMP3*TMP4-2.0D0*(XHF(II,JJ)+EXCAB-EXCAA))**2
      EXCIJ=DSQRT(TMP5)/2.0D0-XHF(II,JJ)
      EXCAA=EXCIJ
print*,tmp1,tmp2,tmp3,tmp4,tmp5
print*,EXCIJ
      END
