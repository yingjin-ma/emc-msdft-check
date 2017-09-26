      SUBROUTINE CALDEXCR(HH,HD,XHFJ,XHFK,NDET,EXC,iMETHOD,NORB,NA, &
                          NB,Energy_XC,FACTOR)
      INTEGER NDET,I,J,K,L,M,N,iMETHOD,NORB,NA,NB
      REAL*8 HD(NDET,NDET),DSIGN
      REAL*8 EXC(NDET,NDET),HH(NDET,NDET)
      REAL*8 XHF(NDET,NDET),XDAT(NDET,NDET)
      REAL*8 EXCR(NDET,NDET),SS(NDET,NDET) 
      REAL*8 DFT(NDET),HF(NDET),DEC(NDET)
      REAL*8 XHFJ(NDET,NDET),XHFK(NDET,NDET)
      REAL*8 F1,Energy_XC(NDET),FACTOR(NDET,NDET)
      INTEGER NE(NDET,NDET)

      DSIGN=1.0D0
      EXC=0.0D0
      FACTOR=0.0D0
      DO I=1,NDET
        EXC(I,I)=Energy_XC(I)
        FACTOR(I,I)=1.0D0
      ENDDO
      XHF=HD+XHFJ+XHFK
      F1=0.0D0
      DO I=1,NDET
        DO J=I+1,NDET
          IF    (iMETHOD==0) THEN
            F1=2.0D0*HH(I,J)/(HH(I,I)+HH(J,J))*DSIGN
          ELSEIF(iMETHOD==1) THEN
            F1=2.0D0*XHF(I,J)/(XHF(I,I)+XHF(J,J))*DSIGN
          ELSEIF(iMETHOD==3) THEN
            F1=2.0D0*XHFJ(I,J)/(XHFJ(I,I)+XHFJ(J,J))*DSIGN
          ELSEIF(iMETHOD==2) THEN
            F1=2.0D0*(XHFJ(I,J)+XHFK(I,J))/(XHFJ(I,I)+XHFJ(J,J) &
               +XHFK(I,I)+XHFK(J,J))*DSIGN
          ELSEIF(iMETHOD==10) THEN
            F1=0.0D0
          ELSEIF(iMETHOD==101) THEN
            CALL METHOD101(HH(I,I),HH(J,J),HH(I,J),EXC(I,I), &
                 EXC(I,J),EXC(I,J))
            FACTOR(I,J)=(2.0D0*EXC(I,J))/(EXC(I,I)+EXC(J,J))
          ELSEIF(iMETHOD==100) THEN
            CALL METHOD100(NORB,NA,NB)
            OPEN(23,FILE='V12.tmp')
            READ(23,*)F1
            CLOSE(23)
            EXC(I,J)=F1-HH(I,J)
            FACTOR(I,J)=(2.0D0*EXC(I,J))/(EXC(I,I)+EXC(J,J))
            WRITE(8406,*)'correlation V12:'
            WRITE(8406,*)I,J,F1,EXC(I,J)
          ELSE
            WRITE(8406,'(A30,I3)')'Error for computed DFTEXCR',iMETHOD
          ENDIF
          IF(iMETHOD<100) THEN
            WRITE(8406,*)'F1=',F1
            EXC(I,J)=F1*(EXC(I,I)+EXC(J,J))/2.0D0
            FACTOR(I,J)=F1
          ENDIF 
          EXC(J,I)=EXC(I,J)
          FACTOR(J,I)=FACTOR(I,J)
        ENDDO
      ENDDO
      END
