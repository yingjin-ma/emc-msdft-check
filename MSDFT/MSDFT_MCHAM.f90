        SUBROUTINE MCHAM(HD,SD,XDAT,XHF,EXCR,EN,DIM,NDET,ICASE,HH,SS)
        INTEGER DIM,NDET,ICASE
        REAL*8 HD(DIM),SD(DIM),XDAT(DIM),XHF(DIM),EXCR(DIM),EN
        REAL*8 HH(NDET,NDET),SS(NDET,NDET)
        real*8 OVERLAP(DIM),EC,HF(NDET,NDET)
        real*8 F11,F12,F21,F22
        real*8 TXHF(NDET,NDET),TXDAT(NDET,NDET),TEXCR(NDET,NDET)
        integer II,i,j,k,l
        character cc

        OVERLAP=SD

        II=0
        DO I = 1,NDET
          DO J = 1,I
            II = II+1
            TXHF(J,I)=XHF(II)
            TXHF(I,J)=TXHF(J,I)
            TXDAT(J,I)=XDAT(II)
            TXDAT(I,J)=TXDAT(J,I)
            TEXCR(J,I)=EXCR(II)
            TEXCR(I,J)=TEXCR(J,I)
            SS(J,I) = SD(II)
            SS(I,J) = SS(J,I)
            IF( I .EQ. J ) THEN
              HH(J,I) = HD(II)+XDAT(II)/2.0D0+EXCR(II)
              HH(J,I) = HH(J,I)*OVERLAP(II)
              HF(J,I) = (HD(II)+XHF(II)/2.00D0)*OVERLAP(II)
            ENDIF
          ENDDO
        ENDDO

        DO I = 1,NDET
           SD(I) = 1.D0/DSQRT(SS(I,I))
        ENDDO

        II = 0
        DO I = 1,NDET
          DO J = 1,I
            II = II+1
            IF( I .NE. J) THEN
              IF(ICASE==10) THEN
                HH(J,I) = (HD(II) + XHF(II)/2.00D0)*OVERLAP(II)
              ELSE IF(ICASE==20) THEN                 
                EC=(HH(J,J)-HF(J,J)+HH(I,I)-HF(I,I))/2.0D0
                HH(J,I)=(HD(II)+XHF(II)/2.00D0)*OVERLAP(II)+SS(J,I)*EC
              ELSE 
                F11=( XHF(II)*(TXDAT(I,I)+2*TEXCR(I,I))/ TXHF(I,I)+& 
     &        XHF(II)*(TXDAT(J,J)+2*TEXCR(J,J))/ TXHF(J,J))/2.0D0
                F12=(XDAT(II)*(TXDAT(I,I)+2*TEXCR(I,I))/TXDAT(I,I)+&
     &        XDAT(II)*(TXDAT(J,J)+2*TEXCR(J,J))/TXDAT(J,J))/2.0D0
                F22=0.75D0*F11/XHF(II)+0.25D0*F12/XDAT(II)
                HH(j,I)=(HD(II)+F22*XHF(II)/2.0D0)*OVERLAP(II)
              END IF
              HH(I,J) = HH(J,I)
            ENDIF
           ENDDO
        ENDDO
        DO I = 1,NDET
          DO J = 1,NDET
            SS(I,J) = SS(I,J)*SD(I)*SD(J)
            HH(I,J) = HH(I,J)*SD(I)*SD(J)
            HH(I,J) = HH(I,J)+SS(I,J)*EN
          ENDDO
        ENDDO

        END
        
