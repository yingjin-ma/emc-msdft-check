      SUBROUTINE CHECKFILE(NUM,iERRO)
      INTEGER NUM,iERRO,I
      CHARACTER*40 CFILE(30)
      iERRO=-1
      CFILE(1)='DETs'
      CFILE(2)='RDMDETs'
      CFILE(3)='INTEGRAL_ORDER_1_IRREPS.grp'
      CFILE(4)='INTEGRAL_ORDER_2_IRREPS.grp'
      CFILE(5)='INTEGRAL_ORDER_4_IRREPS.grp'
      CFILE(6)='INTEGRAL_ORDER_8_IRREPS.grp'
      CFILE(7)='MASORB.orb'
      CFILE(8)='infile.inp'
      CFILE(14)='GAMESSMO.dat'
      CFILE(21)='INT_MO'
      CFILE(22)='FCIDUMP'
      CFILE(23)='ONETWOINT'


      DO I=1,8       
        iERRO=-1
        OPEN(23,FILE=CFILE(I))
        DO WHILE(.TRUE.)
          READ(23,*,END=301)
          iERRO=iERRO+1
        ENDDO
301     CLOSE(23)
        IF(iERRO==-1) THEN
          WRITE(*,*)'ERRO: NO ',CFILE(I) 
          GOTO 399
        ENDIF
      ENDDO

      IF(NUM==0) THEN
       DO I=14,14       
        iERRO=-1
        OPEN(23,FILE=CFILE(I))
        DO WHILE(.TRUE.)
          READ(23,*,END=302)
          iERRO=iERRO+1
        ENDDO
302     CLOSE(23)
        IF(iERRO==-1) THEN
          WRITE(*,*)'ERRO: NO ',CFILE(I)
          GOTO 399
        ENDIF
       ENDDO
      ELSEIF(NUM==1) THEN
       DO I=21,23
        iERRO=-1
        OPEN(23,FILE=CFILE(I))
        DO WHILE(.TRUE.)
          READ(23,*,END=303)
          iERRO=iERRO+1
        ENDDO
303     CLOSE(23)
        IF(iERRO==-1) THEN
          WRITE(*,*)'ERRO: NO ',CFILE(I)
          GOTO 399
        ENDIF
       ENDDO
      ENDIF
399   CONTINUE
      END
