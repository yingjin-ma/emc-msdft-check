      SUBROUTINE RHOAB(NAO,NMO,NA,NB,MOA,MOB,RHOA,RHOB)
      INTEGER NAO,NMO,NA,NB
      REAL*8 MOA(NAO,NMO),MOB(NAO,NMO)
      REAL*8 RHOA(NAO,NAO),RHOB(NAO,NAO)
      REAL*8 TRHO(NAO,NAO)
      INTEGER I,J,K,L
      RHOA=0.0D0
      DO I=1,NA
        DO J=1,NAO
        DO K=1,NAO
          TRHO(J,K)=MOA(J,I)*MOA(K,I)
        ENDDO
        ENDDO
        RHOA=RHOA+TRHO
      ENDDO
      RHOB=0.0D0
      DO I=1,NB
        DO J=1,NAO
        DO K=1,NAO
          TRHO(J,K)=MOB(J,I)*MOB(K,I)
        ENDDO
        ENDDO
        RHOB=RHOB+TRHO
      ENDDO 
      END SUBROUTINE
