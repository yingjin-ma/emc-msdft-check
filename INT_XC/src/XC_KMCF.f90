      SUBROUTINE XC_KMCF(IGRID,NAO,RHOA,RHOB,TRHOA,TRHOB,KAPA)
      USE GRID_INFO
      IMPLICIT NONE
      INTEGER IGRID,NAO ! input
      REAL*8 RHOA(NAO,NAO),RHOB(NAO,NAO) ! input
       REAL*8 TRHOA(NAO,NAO),TRHOB(NAO,NAO) ! input
      REAL*8 KAPA ! output
      INTEGER I,J,K
      REAL*8 PI,BMN(6,5),RHO,TRHO,RS,KK,XX

      ! setting parameters
      PI=3.141592653589793D0 
      BMN(1,1)=-0.2207193D1
      BMN(2,1)= 0.1128469D1
      BMN(3,1)=-0.2475593D0
      BMN(4,1)= 0.8616560D-1
      BMN(5,1)=-0.6500077D-2
      BMN(6,1)=-0.2491486D-2
      BMN(1,2)= 0.6807648D1
      BMN(2,2)=-0.2535669D1
      BMN(3,2)= 0.4243142D0
      BMN(4,2)=-0.1715714D0
      BMN(5,2)= 0.1714085D-1
      BMN(6,2)= 0.5321373D-2
      BMN(1,3)=-0.6386316D1
      BMN(2,3)= 0.2432821D1
      BMN(3,3)=-0.2565175D0
      BMN(4,3)= 0.1067547D0
      BMN(5,3)=-0.1462187D-1
      BMN(6,3)=-0.3704699D-2
      BMN(1,4)= 0.2860522D1
      BMN(2,4)=-0.1064058D1
      BMN(3,4)= 0.8294749D-1
      BMN(4,4)=-0.2392882D-1
      BMN(5,4)= 0.4423830D-2
      BMN(6,4)= 0.9700054D-3
      BMN(1,5)=-0.7466076D-1
      BMN(2,5)= 0.3843687D-1
      BMN(3,5)=-0.3184296D-2
      BMN(4,5)= 0.2579856D-2
      BMN(5,5)=-0.4427570D-3
      BMN(6,5)=-0.9518308D-4

      ! compute 2*phi*phi at r
      TRHO=0.0D0
      I=IGRID
      DO J=1,NAO
        DO K =1,NAO
         TRHO=TRHO+(TRHOA(J,K)+TRHOB(J,K)) &
            *GRIDS(I)%VAL0(J)*GRIDS(I)%VAL0(K)*GRIDS(I)%WEIGHT
        ENDDO
      ENDDO
      ! compute rho at r
      RHO=0.0D0
      I=IGRID
      DO J=1,NAO
        DO K =1,NAO
         RHO=RHO+(RHOA(J,K)+RHOB(J,K)) &
            *GRIDS(I)%VAL0(J)*GRIDS(I)%VAL0(K)*GRIDS(I)%WEIGHT
        ENDDO
      ENDDO
      RHO =dmax1(0.0D0,RHO )
      TRHO=dmax1(0.0D0,TRHO)
      IF(RHO>1.0D-20 .AND. TRHO>RHO) THEN
      ! compute k,rs and x at r
      KK=(TRHO/RHO)**(1.0D0/3.0D0)
      RS=(3.0D0/(4.0D0*PI*RHO))**(1.0D0/3.0D0)
      IF(RS>25.0D0) RS=25.0D0
      XX=LOG(RS)
      ! compute facter ka at r
      KAPA=0.0D0
      DO I=1,6
        DO J=1,5
          KAPA=KAPA+(BMN(I,J)*(XX**(I-1))*(KK**(J-1)))
        ENDDO
      ENDDO
      KAPA=1.0D0/KAPA
      ELSE
      KAPA=1.0D0
      ENDIF
      IF(RHO==TRHO) KAPA=1.0D0
      END
      
