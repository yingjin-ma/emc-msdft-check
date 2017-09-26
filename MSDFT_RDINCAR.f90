      SUBROUTINE RD_INCAR(iPARA,FILENAME,DIM,iERRO)
!C#################################################################
!C DIM is the dimenstion of iPARA
!C The initial parameters are saved in iPARA
!C The index of iPARA as follow,
!C iPARA( 1): NAO, total atomic orbitals
!C iPARA( 2): NMO, total molecular orbitals
!C iPARA( 3): FROZEN, frozen orbitals
!C iPARA( 4): CLOSED, closed orbitals
!C iPARA( 5): ACTIVE, active orbitals
!C iPARA( 6): OCC, occupied orbitals
!C iPARA( 7): NATOM, number of atoms
!C iPARA(11): ELE, total electrons
!C iPARA(12): AELE, active electrons
!C iPARA(13): AE, active alpha electrons
!C iPARA(14): BE, active beta  electrons
!C iPARA(15): SPIN, 2S+1
!C iPARA(21): DETS, configurations
!C iPARA(22): CONFs, number of selected configurations
!C iPARA(31): STATES, target states
!C iPARA(32): IWEIGHT, specify the wights (1) or not (0)
!C iPARA(33): iROOT, the number of roots for TQL or Davidson
!C iPARA(41-60): INT, weights of each states (less than 20 states)
!C iPARA(61): IMAX, max stept of MCSCF
!C iPARA(91): RESTART, restart (1) or not (0)
!C iPARA(92): METHOD, CASSCF (0) or MSDFTSCF (1)
!C iPARA(93): Initial Guess for MOs, GAMESS(0),MOLPRO(1),GAUSSIAN(2)
!C iPARA(94): Optimized orbitals (0) or not (1),OPTORB
!C iPARA(95): iMETHOD, (0)HH (1)XHF (2)J+K (3)J
!C iPARA(96): iTOLER, tolerance, default is 15 (10D-15)

!C FILENAME( 1): The name of output file
!C################################################################

      INTEGER DIM
      INTEGER iPARA(DIM)
      CHARACTER*40 FILENAME(DIM)
      CHARACTER*20 CC
      INTEGER DATA,NUM,iLINE,iERRO

!C Setting default values for the initial parameters
      iPARA=0
      IMAX=1
      STATES=1
      iLINE=0
      iERRO=0
      OPEN(23,FILE='INCAR')
      READ(23,*,END=101)FILENAME(1)
      READ(23,*,END=101)
      iLINE=iLINE+1
      DO WHILE(.TRUE.)
        READ(23,*,END=101)CC,DATA
        iLINE=iLINE+1
        CALL iFIND(CC,NUM)
        IF(NUM==-1) THEN
          iLINE=-1
          GOTO 101
        ENDIF
        iPARA(NUM)=DATA
      ENDDO
101   CLOSE(23)
      IF(iLINE==-1) THEN
        iERRO=-1
        WRITE(*,*)'ERROR IN RD_INCAR: READ ERROR'
        GOTO 111
      ELSEIF(iLINE==0) THEN
        CALL iWRITE()
        iERRO=-1
        WRITE(*,*)'ERROR IN RD_INCAR: NO INCAR'
        GOTO 111
      ENDIF

      IF(iPARA(2)==0) iPARA(2)=iPARA(1)
      iPARA(12)=iPARA(11)-2*iPARA(3)-2*iPARA(4)
      iPARA(13)=(iPARA(12)+iPARA(15)-1)/2 
      iPARA(14)=iPARA(12)-iPARA(13)
      iPARA(6)=iPARA(3)+iPARA(4)+iPARA(5)
111   CONTINUE
      END

      SUBROUTINE iWRITE()
      OPEN(24,FILE='INCAR')
      WRITE(24,*)'output_file_name'
      WRITE(24,*)'basis set'
      WRITE(24,*)'METHOD  0    !INT,(0)CASSCF,(1)MSDFT'
      WRITE(24,*)'NATOM   0    !INT,number of atoms'
      WRITE(24,*)'NAO     0    !INT,number of atomic orbitals' 
      WRITE(24,*)'FROZEN  0    !INT,frozen orbitals'
      WRITE(24,*)'CLOSED  0    !INT,closed orbitals'
      WRITE(24,*)'ACTIVE  0    !INT,active orbitals'
      WRITE(24,*)'ELE     0    !INT,total electrons'
      WRITE(24,*)'SPIN    1    !INT,2S+1 eg,singlet:1'
      WRITE(24,*)'STATES  1    !INT,number of target states'
      WRITE(24,*)'IMAX    1    !INT,the max step for orb opt'
      WRITE(24,*)'RESTART 0    !INT,(1)restart or (0)not'
      WRITE(24,*)'OPTORB  1    !INT,(0)orbital opt or (1)not'
      WRITE(24,*)'IWEIGHT 0    !INT,(1)read from WEIGHTS or (0)not'
      WRITE(24,*)'iMETHOD 0    !INT,(0)HH (1)XHF (2)J+K (3)J'
      WRITE(24,*)'CONFs   0    !INT,(0)All DETs (>0) Number of DETs'
      WRITE(24,*)'iROOT   0    !INT,Number of roots for TQL'
      CLOSE(24)
      END

      SUBROUTINE iFIND(CC,NUM)
      CHARACTER*20 CC
      INTEGER NUM
      IF(CC=='METHOD') THEN
        NUM=92
      ELSEIF(CC=='RESTART') THEN
        NUM=91
      ELSEIF(CC=='OPTORB') THEN
        NUM=94
      ELSEIF(CC=='IMAX') THEN
        NUM=61
      ELSEIF(CC=='NAO') THEN
        NUM=1
      ELSEIF(CC=='NMO') THEN
        NUM=2
      ELSEIF(CC=='FROZEN') THEN
        NUM=3
      ELSEIF(CC=='CLOSED') THEN
        NUM=4
      ELSEIF(CC=='ACTIVE') THEN
        NUM=5
      ELSEIF(CC=='ELE') THEN
        NUM=11
      ELSEIF(CC=='SPIN') THEN
        NUM=15
      ELSEIF(CC=='STATES') THEN
        NUM=31
      ELSEIF(CC=='MOGUESS') THEN
        NUM=93
      ELSEIF(CC=='iMETHOD') THEN
        NUM=95
      ELSEIF(CC=='iTOLER') THEN
        NUM=96
      ELSEIF(CC=='CONFs') THEN
        NUM=22
      ELSEIF(CC=='NATOM') THEN
        NUM=7
      ELSEIF(CC=='iROOT') THEN
        NUM=33
      ELSEIF(CC=='IWEIGHT') THEN
        NUM=32
      ELSE 
        NUM=-1
      ENDIF   
      END   
