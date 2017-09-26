      SUBROUTINE MCCONF(DIM,NDET,NORB,IERRO)
      INTEGER DIM,NDET,NORB,IERRO
      INTEGER DETs(NDET,2,NORB)
      INTEGER CONFs(DIM)
      INTEGER I,J,K,NA,NB

      IERRO=0
      OPEN(23,FILE='DETs')
      READ(23,*)NA,NB
      DO I=1,NDET
        READ(23,*)(DETs(I,1,J),J=1,NA)
        READ(23,*)(DETs(I,2,J),J=1,NB)
      ENDDO
      CLOSE(23)

      OPEN(24,FILE='CONFs')
      READ(24,*)(CONFs(I),I=1,DIM)
      CLOSE(24)

      OPEN(25,FILE='DETs')
      WRITE(25,*)NA,NB
      DO K=1,DIM
        I=CONFs(K)
        WRITE(25,*)(DETs(I,1,J),J=1,NA)
        WRITE(25,*)(DETs(I,2,J),J=1,NB)
      ENDDO
      CLOSE(25)
      END

      

