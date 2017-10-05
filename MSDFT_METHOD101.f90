      SUBROUTINE METHOD101(HH11,HH22,HH12,EXC1, &
                           EXC2,EXC12)
      IMPLICIT NONE
!     dE(MSDFT)=dE(WFT)+dE(Ec1-Ec2)
      REAL*8 HH11,HH22,HH12,EXC1,EXC2,EXC12
      REAL*8 DH,DE
      REAL*8 A,B

      DH=DABS(HH11-HH22)
      DE=DABS(EXC1-EXC2)
      A=DSQRT(DH**2+4.0D0*HH12)+DE-DH
      B=0.5D0*DSQRT(A**2-DE**2)
      IF(HH12>0.0D0) THEN
        EXC12=B-HH12
      ELSEIF(HH12<0.0D0) THEN
        EXC12=-B-HH12
      ELSE
        EXC12=0.0D0
      ENDIF
      END
