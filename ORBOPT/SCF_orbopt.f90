!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Test Program main
!Do orbital optimization
!All the input data in Program main is given by MSDFT (Zexing Qu)
!The SCF_orbopt is used for updata the orbitals (MO) (Yingjin Ma)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!SUBROUTINE CALL_ORBOPT()
!USE GLOBAL_PARA
!
Program main
integer,parameter :: n_orbitals=1
integer,parameter :: n_occ=1
integer n_atoms,icharge,ispin,n_electrons,n_closed,n_active
real*8  oneRDM(n_occ,n_occ),twoRDM(n_occ,n_occ,n_occ,n_occ)
real*8  oneEI(n_orbitals,n_orbitals),twoEI(n_orbitals,n_orbitals,n_orbitals,n_orbitals)
real*8  Fc(n_orbitals,n_orbitals)
real*8  MO(n_orbitals,n_orbitals)
!++++++++++++++++++++++++++++++++++++++++++
! MO will updata according to SCF_orbopt
!++++++++++++++++++++++++++++++++++++++++++
call SCF_orbopt(n_atoms,n_orbitals,icharge,ispin,n_electrons,n_closed,n_active,n_occ, &
                Fc,oneRDM,twoRDM,oneEI,twoEI,MO)
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine SCF_orbopt(n_atoms,n_orbitals,icharge,ispin,n_electrons,n_closed,n_active,n_occ, &
                      Fc,oneRDM,twoRDM,oneEI,twoEI,MO)
integer n_atoms,n_orbitals,icharge,ispin,n_electrons,n_closed,n_active,n_occ
real*8  oneRDM(n_occ,n_occ),twoRDM(n_occ,n_occ,n_occ,n_occ)
real*8  oneEI(n_orbitals,n_orbitals),twoEI(n_orbitals,n_orbitals,n_orbitals,n_orbitals)
real*8  MO(n_orbitals,n_orbitals)
!++++++++++++++++++++++++++++++++++++++++++++++++
! Fc is a N*N matrix including correlation energy 
! which is simiar with oneEI
!++++++++++++++++++++++++++++++++++++++++++++++++
real*8  Fc(n_orbitals,n_orbitals)

End
