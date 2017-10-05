      SUBROUTINE XC_ZETA(RHOA,RHOB)
      IMPLICIT NONE
      REAL*8 RHOA,RHOB ! input and output
      REAL*8 RHO,FF,GG,AA,BB,RS,PI,ZETA
      AA=0.194D0
      BB=0.525D0
      PI=3.141592653589793D0
      ZETA=RHOA-RHOB
      RHO=RHOA+RHOB
      IF(RHO>1.0D-20) THEN
        RS=(3.0D0/(4.0D0*PI*RHO))**(1.0D0/3.0D0)
        GG=0.5D0*(1.0D0+2.0D0*AA*RS)/(1.0D0+BB*RS+AA*BB*RS*RS)**2
        FF=DSQRT(DABS(1.0D0-4.0D0*GG*(1.0D0-ZETA**2)))
      ELSE
        FF=1.0D0
      ENDIF
      RHOA=(RHO/2.0D0)*(1.0D0+FF)
      RHOB=(RHO/2.0D0)*(1.0D0-FF)  
      END
