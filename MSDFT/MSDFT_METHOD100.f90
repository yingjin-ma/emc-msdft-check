      SUBROUTINE METHOD100(NORB,NA,NB)
      INTEGER NORB,NACT,NA,NB
      INTEGER INDEXA(NA),INDEXB(NB)
      INTEGER I,J,K
      REAL*8 COEFF(NORB,NORB),EAA,EABBA

      DO I=1,NA
        INDEXA(I)=I
      ENDDO
      DO I=1,NB
        INDEXB(I)=I
      ENDDO
      INDEXB(NB)=NB+1
!C READ Initial MOs
      OPEN(23,FILE='INT_MO')
      DO I=1,NORB
        READ(23,*)(COEFF(J,I),J=1,NORB)
      ENDDO
      CLOSE(23)
!C WRITE MOABBA.tmp
      OPEN(24,FILE='MOABBA.tmp')
      DO J=1,NA
        K=INDEXA(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      DO J=1,NB
        K=INDEXB(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      DO J=1,NB
        K=INDEXB(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      DO J=1,NA
        K=INDEXA(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      CLOSE(24)
!C WRITE MOAA.tmp
      OPEN(24,FILE='MOAA.tmp')
      DO J=1,NA
        K=INDEXA(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      K=INDEXB(NB)
      WRITE(24,*)COEFF(:,K)
      DO J=1,NB-1
        K=INDEXB(J)
        WRITE(24,*)COEFF(:,K)
      ENDDO
      CLOSE(24)
!C RUN EABBA
      CALL SYSTEM('$MOVB/runabba BLWCIABBA.inp')
!C RUN EAA
      CALL SYSTEM('$MOVB/runaa BLWCIAA.inp')
!C WRITE V12.tmp
      OPEN(23,FILE='EABBA.tmp')
      READ(23,*)EABBA
      CLOSE(23)
      OPEN(23,FILE='EAA.tmp')
      READ(23,*)EAA
      CLOSE(23)
      OPEN(24,FILE='V12.tmp')
      WRITE(24,*)EABBA-EAA
      CLOSE(24)
      END
      SUBROUTINE METHOD101(HH11,HH22,HH12,EXC1, &
                           EXC2,EXC12)
      IMPLICIT NONE
      REAL*8 HH11,HH22,HH12,EXC1,EXC2,EXC12
      REAL*8 A,B,C,DT

      DT=DSQRT((HH11-HH22)**2+4.0D0*(HH12**2))
      A=4.0D0
      B=8.0D0*HH12
      IF(HH11>HH22) THEN
        C=4.0D0*(HH12**2)+(HH11-HH22+EXC1-EXC2)**2-(DT+EXC1-EXC2)**2
      ELSEIF(HH11<HH22) THEN
        C=4.0D0*(HH12**2)+(HH11-HH22+EXC1-EXC2)**2-(DT+EXC2-EXC1)**2
      ELSE
        C=0.0D0
      ENDIF
      DT=B**2-4.0D0*A*C
      IF(DT<0.0D0) THEN
         DT=0.0D0
      ENDIF
      IF(HH12>0.0D0) THEN
        EXC12=0.5D0*(-B+DSQRT(DT))
      ELSEIF(HH12<0.0D0) THEN
        EXC12=0.5D0*(-B-DSQRT(DT))
      ELSE
        EXC12=0.0D0
      ENDIF
      END
