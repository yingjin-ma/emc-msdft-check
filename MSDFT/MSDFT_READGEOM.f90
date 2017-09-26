      SUBROUTINE RD_GEOM(NATOM,GEOM,ATOMCHG,iERRO)
      IMPLICIT NONE
      INTEGER NATOM,iERRO,I
      INTEGER ATOMCHG(NATOM)
      REAL*8  GEOM(NATOM,3)
      CHARACTER*10 CTYPE(NATOM)
      GEOM=0.0D0
      ATOMCHG=0
      OPEN(23,FILE='GEOM.xyz')
      DO I=1,NATOM
        READ(23,*),CTYPE(I),GEOM(I,1),GEOM(I,2),GEOM(I,3)
      ENDDO
      CLOSE(23)
      DO I=1,NATOM
        IF(CTYPE(I)=='H'  .or. CTYPE(I)=='h')  ATOMCHG(I)=1
        IF(CTYPE(I)=='He' .or. CTYPE(I)=='he') ATOMCHG(I)=2
        IF(CTYPE(I)=='Li' .or. CTYPE(I)=='li') ATOMCHG(I)=3
        IF(CTYPE(I)=='Be' .or. CTYPE(I)=='be') ATOMCHG(I)=4
        IF(CTYPE(I)=='B'  .or. CTYPE(I)=='b')  ATOMCHG(I)=5
        IF(CTYPE(I)=='C'  .or. CTYPE(I)=='c')  ATOMCHG(I)=6
        IF(CTYPE(I)=='N'  .or. CTYPE(I)=='n')  ATOMCHG(I)=7
        IF(CTYPE(I)=='O'  .or. CTYPE(I)=='o')  ATOMCHG(I)=8
        IF(CTYPE(I)=='F'  .or. CTYPE(I)=='f')  ATOMCHG(I)=9
        IF(CTYPE(I)=='Ne' .or. CTYPE(I)=='ne') ATOMCHG(I)=10
        IF(CTYPE(I)=='Na'  .or. CTYPE(I)=='na')  ATOMCHG(I)=11
        IF(CTYPE(I)=='Mg' .or. CTYPE(I)=='mg') ATOMCHG(I)=12
        IF(CTYPE(I)=='Al' .or. CTYPE(I)=='al') ATOMCHG(I)=13
        IF(CTYPE(I)=='Si' .or. CTYPE(I)=='si') ATOMCHG(I)=14
        IF(CTYPE(I)=='P'  .or. CTYPE(I)=='p')  ATOMCHG(I)=15
        IF(CTYPE(I)=='S'  .or. CTYPE(I)=='s')  ATOMCHG(I)=16
        IF(CTYPE(I)=='Cl'  .or. CTYPE(I)=='cl')  ATOMCHG(I)=17
        IF(CTYPE(I)=='Ar'  .or. CTYPE(I)=='ar')  ATOMCHG(I)=18
        IF(CTYPE(I)=='K'  .or. CTYPE(I)=='k')  ATOMCHG(I)=19
        IF(CTYPE(I)=='Ca' .or. CTYPE(I)=='Ca') ATOMCHG(I)=20
      ENDDO      
      END SUBROUTINE RD_GEOM
