      SUBROUTINE MSDFT_Fc(NORB,NCONF,COEFF,DFTFCA,DFTFCB,FF,EXC, &
                 Energy_XC,IMETHOD,iROOT,iSTATE,WEIGHTS)
      IMPLICIT NONE
      INTEGER NORB,NCONF,IMETHOD,iROOT,iSTATE
      REAL*8 WEIGHTS(iSTATE)
      REAL*8 COEFF(NCONF,NCONF),EXC(NCONF,NCONF)
      REAL*8 FF(NCONF,NCONF),Energy_XC
      REAL*8 DFTFCA(NORB,NORB,NCONF),DFTFCB(NORB,NORB,NCONF)
      REAL*8 FC(NORB,NORB,NCONF,NCONF)
      REAL*8 TMP(NORB,NORB),MO(NORB,NORB)
      REAL*8 MTR1(NORB,NORB),DSUM
      INTEGER I,J,K,L,IS,I1,I2
      integer       :: conf_index(iSTATE)
      DO I=1,NCONF
        FC(:,:,I,I)=DFTFCA(:,:,I)+DFTFCB(:,:,I)
      ENDDO

      DO I=1,NCONF
        DO J=I+1,NCONF
          DO K=1,NORB
          DO L=1,NORB
            FC(K,L,I,J)=FF(I,J)*(FC(K,L,I,I)+FC(K,L,J,J))/2.0D0
            FC(K,L,J,I)=FC(K,L,I,J)
          ENDDO
          ENDDO
        ENDDO
      ENDDO


      TMP=0.0D0
      DSUM=0.0D0

      conf_index=0
      IF(ABS(iROOT-iSTATE)>0) THEN
        call CONF_PURE(Coeff(:,1:iROOT),NCONF,iROOT,iSTATE,conf_index)
      ENDIF
      DO I2=1,iSTATE
        IF(ABS(iROOT-iSTATE)>0) THEN
          I1=conf_index(I2)
        ELSE
          I1=I2
        ENDIF
        DO I=1,NCONF
          DO J=1,NCONF
            TMP=TMP+COEFF(I,I1)*COEFF(J,I1)*FC(:,:,I,J)*WEIGHTS(I2)
            DSUM=DSUM+COEFF(I,I1)*COEFF(J,I1)*EXC(I,J)*WEIGHTS(I2) 
          ENDDO
        ENDDO
      ENDDO

      IF(IMETHOD==1) THEN
        TMP=0.0D0
        TMP=TMP+DFTFCA(:,:,1)+DFTFCB(:,:,1)
      ENDIF
      TMP=TMP*0.5D0
      OPEN(23,FILE='INT_MO')
      DO I=1,NORB
        READ(23,*)(MO(J,I),J=1,NORB)
      END DO
      CLOSE(23)


      MTR1=MATMUL(TRANSPOSE(MO),TMP)
      TMP=MATMUL(MTR1,MO)

! Print (Fc)ij
      OPEN(34,FILE='FCIJ.tmp')
      WRITE(34,*)TMP
      CLOSE(34)          
! Print EXC
!      OPEN(35,FILE='Corr_ENG.tmp')
!      WRITE(35,*)DSUM,Energy_XC
!      CLOSE(35)
      END
