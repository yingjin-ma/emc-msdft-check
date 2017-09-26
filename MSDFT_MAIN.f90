!C*******************************************************************************
!C
!C    MSDFT-CASSCF VERSION 1.0
!C             Jilin Unv. Zexing Qu, 2017.9.4
!C
!C*******************************************************************************
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER,PARAMETER :: DIM=100
      INTEGER iPARA(DIM),iERRO
      CHARACTER*40 cPARA(DIM)
      REAL*8,ALLOCATABLE :: GEOM_0(:,:)
      INTEGER,ALLOCATABLE :: iATOMCHG(:)
      INTEGER IN1,IN2,IN3,IN4,IN5
      INTEGER OUT1

      CALL RD_INCAR(iPARA,cPARA,DIM,iERRO)
      IF(IERRO==-1) GOTO 9999
      IN1=iPARA(7)
      ALLOCATE(GEOM_0(IN1,3))
      ALLOCATE(iATOMCHG(IN1))
      CALL RD_GEOM(IN1,GEOM_0,iATOMCHG,iERRO)
      IF(IERRO==-1) GOTO 9990
      OPEN (8406,FILE=cPARA(1))
      WRITE(8406,*)''
      WRITE(8406,*)'**************************************************'
      WRITE(8406,*)'                                                  '
      WRITE(8406,*)'           MSDFT-SCF VERSION 1.0                  '
      WRITE(8406,*)'             Jilin Univ.  Zexing Qu               '
      WRITE(8406,*)'                                    2017.9.4      '
      WRITE(8406,*)'                                                  '
      WRITE(8406,*)'**************************************************'
      WRITE(8406,*)''
      WRITE(8406,*)'-------- Entering the  MSDFT-SCF Programm  -------'
      WRITE(8406,*)''
      WRITE(8406,*)'------------- STEP 1: Initialization -------------'
      WRITE(8406,*)''
      IN1=iPARA(5)
      IN2=iPARA(13)
      IN3=iPARA(14)
      IN4=iPARA(3)
      IN5=iPARA(4)
      WRITE(8406,*)' Start 1.1 '
      WRITE(8406,*)' Generate CAS configuration'
      CALL CASCONF(IN1,IN2,IN3,IN4,IN5,OUT1,IERRO)
      IF(IERRO==-1) GOTO 9900
      iPARA(21)=OUT1
      IN1=iPARA(22)
      IF(IN1>0 .AND. IN1<=iPARA(21)) THEN
        IN2=iPARA(21)
        IN3=iPARA(1)
        CALL MCCONF(IN1,IN2,IN3,IERRO)
        IF(IERRO==-1) GOTO 9900
        iPARA(21)=IN1
      ENDIF
      CALL find_couple(iPARA(21),iPARA(3)+iPARA(4)+iPARA(5),ierro)
      IF(IERRO==-1) GOTO 9900
      WRITE(8406,*)' End 1.1'
      WRITE(8406,*)''
!######################################################################
!C need to improve
      IN1=iPARA(11)
      IN2=iPARA(6)
      IN3=iPARA(3)
      IN4=iPARA(4)
      WRITE(8406,*)''
      WRITE(8406,*)'Generate the input files for orbital optimization'
      CALL MCFILE(IN1,IN2,IN3,IN4,IERRO)
      WRITE(8406,*)''
      IF(iPARA(91)==0) THEN
!        IF(iPARA(93)==0) THEN
!          CALL Gamess_MO(NAO,NMO,IORDER,NMIN,NUM,CF) 
!        ENDIF
      ENDIF
!######################################################################
      WRITE(8406,*)'------------------- END STEP 1 -------------------'
      WRITE(8406,*)''
      WRITE(8406,*)'------------- STEP 2: CAS(MSDFT)SCF --------------'
      WRITE(8406,*)''
      CALL INTXC_ALLOCATE(iPARA(7))
      CALL BLWCAS(iPARA,DIM,GEOM_0,iATOMCHG,iERRO)
      CALL INTXC_DEALLOCATE()
      WRITE(8406,*)''
      WRITE(8406,*)'------------------- END STEP 2 -------------------'
9900  CONTINUE
      IF(IERRO==-1) THEN
       WRITE(8406,*)'**************** ERROR TERMINATION ***************'
      ELSE
       WRITE(8406,*)'*************** NORMAL TERMINATION ***************'
      ENDIF
      CLOSE(8406)
9990  CONTINUE
      DEALLOCATE(GEOM_0)
      DEALLOCATE(iATOMCHG)
9999  CONTINUE
      END
