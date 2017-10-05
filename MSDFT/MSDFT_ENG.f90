        SUBROUTINE MSDFTENG(HH,Coeff,Dim)
        INTEGER DIM,I,J
        REAL*8 HH(Dim,Dim),COEFF(Dim,Dim)
        REAL*8 D1(Dim),D2(Dim)
        COEFF=HH
        OPEN(23,FILE='HAM.tmp')
        WRITE(*,*)'MSDFT HAM:'
        DO I=1,DIM
          DO J=1,DIM
            WRITE(*,*)I,J,HH(I,J)
            IF(I==J) THEN
              WRITE(23,*)I,J,'1.00',HH(I,J)
            ELSE
              WRITE(23,*)I,J,'0.00',HH(I,J)
            ENDIF
          END DO
        END DO
        CLOSE(23)
        call tql(dim,dim,COEFF,D1,D2)
        DO I=1,Dim
          WRITE(*,*)'MSDFT ENG:',I,D1(I)
          WRITE(8406,'(A12,I5,F25.12)')'MSDFT ENG:',I,D1(I)
        ENDDO
        OPEN(23,FILE='Energy.tmp')
        WRITE(23,*)D1(1)
        CLOSE(23)
        DO I=1,DIM
          WRITE(*,*)'CI COEFF'
          WRITE(*,*)COEFF(:,I)
          WRITE(8406,'(A27,I5)')'CI Coefficiences for state',I
          WRITE(8406,*)COEFF(:,I)
        ENDDO
        END

