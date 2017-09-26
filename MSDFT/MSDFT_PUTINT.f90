      SUBROUTINE PUT_INT(iPARA,DIM)
      INTEGER DIM
      INTEGER iPARA(DIM)
      WRITE(8406,*)''
      WRITE(8406,*)'Initial Parameters:'
      WRITE(8406,'(A33,I9)')'NAO, total atomic orbitals'      ,iPARA( 1)
      WRITE(8406,'(A33,I9)')'NMO, total molecular orbitals'   ,iPARA( 2)
      WRITE(8406,'(A33,I9)')'FROZEN, frozen orbitals'         ,iPARA( 3)
      WRITE(8406,'(A33,I9)')'CLOSED, closed orbitals'         ,iPARA( 4)
      WRITE(8406,'(A33,I9)')'ACTIVE, active orbitals'         ,iPARA( 5)
      WRITE(8406,'(A33,I9)')'OCC, occupied orbitals'          ,iPARA( 6)
      WRITE(8406,'(A33,I9)')'ELE, total electrons'            ,iPARA(11)
      WRITE(8406,'(A33,I9)')'AELE, active electrons'          ,iPARA(12)
      WRITE(8406,'(A33,I9)')'AE, active alpha electrons'      ,iPARA(13)
      WRITE(8406,'(A33,I9)')'BE, active beta  electrons'      ,iPARA(14)
      WRITE(8406,'(A33,I9)')'SPIN, 2S+1'                      ,iPARA(15)
      WRITE(8406,'(A33,I9)')'DETS, configurations'            ,iPARA(21)
      WRITE(8406,'(A33,I9)')'STATES, target states'           ,iPARA(31)
      WRITE(8406,'(A33,I9)')'IWEIGHT, read wights(1)or not(0)',iPARA(32)
      WRITE(8406,'(A33,I9)')'IMAX, max stept of MCSCF'        ,iPARA(61)
      WRITE(8406,'(A33,I9)')'RESTART, restart (1) or not (0)' ,iPARA(91)
      WRITE(8406,'(A33,I9)')'METHOD, CASSCF(0) or MSDFTSCF(1)',iPARA(92)
      WRITE(8406,'(A33,I9)')'ORBOPT, optimized(0) or not(1)'  ,iPARA(94)
      WRITE(8406,'(A33,I9)')'MO GUESS,GAMESS(0),MOLPRO(1)'    ,iPARA(93)

      WRITE(8406,*)''
      END
