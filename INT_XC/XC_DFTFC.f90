subroutine DFT_Fc(ncenters,nconts,nstates,Pa,Pb,Energy_XC,coord,atomchg, &
                  imethod,Coeff,iROOT,iSTATE,WEIGHTS)
implicit none
integer       :: ncenters,nconts,elea,eleb,ele,nstates
integer       :: iROOT,iSTATE,imethod
integer       :: conf_index(iSTATE)
real*8        :: WEIGHTS(iSTATE)
integer       :: imult,icharge,atomchg(ncenters) ! natoms, multiplicity, charge
character     :: functional_in*30  ! functional
real(8)       :: coord(ncenters,3) ! coordinate, atomcharge
real(8)       :: INT_MO(nconts,nconts),Energy_XC(nstates)
! output data zone
real(8)       :: Pa(nconts,nconts,nstates),Pb(nconts,nconts,nstates)
real(8)       :: E_ex,E_corr,Coeff(nstates,nstates)
real(8)       :: MOA(nconts,nconts),MOB(nconts,nconts)   
integer       :: i,j,k,l,DETS(nstates,2,nconts)
real*8        :: TPa(nconts,nconts,iROOT),TPb(nconts,nconts,iROOT)
integer       :: I1,I2
 
imult = 1
icharge = 0
ele=0
do i=1,ncenters
  ele=ele+atomchg(i)
enddo
ele=ele+icharge
elea=(ele+(imult-1))/2
eleb=(ele-(imult-1))/2

!******************************************
! READ DETs and INT_MO from disk
DETS=0
OPEN(23,file='DETs')
READ(23,*)elea,eleb
DO i=1,nstates
  READ(23,*)(DETS(i,1,j),j=1,elea)
  READ(23,*)(DETS(i,2,j),j=1,eleb)
ENDDO
CLOSE(23)
INT_MO=0.0D0
OPEN(23,FILE='INT_MO')
DO I=1,nconts
  READ(23,*)(INT_MO(J,I),J=1,nconts)
ENDDO
CLOSE(23)
!*****************************************
do i=1,nstates
  MOA=INT_MO
  do j=1,elea
    k=DETS(i,1,j)
    MOA(:,j)=INT_MO(:,k)
  enddo
  MOB=INT_MO
  do j=1,eleb
    k=DETS(i,2,j)
    MOB(:,j)=INT_MO(:,k)
  enddo
  call rhoab(nconts,nconts,elea,eleb,MOA,MOB,Pa(:,:,i),Pb(:,:,i))
enddo
functional_in ='PBE_PBE'

print*,ncenters,nconts

if(imethod==1) then
  TPa=0.0D0
  TPb=0.0D0
  do j=1,iROOT
    do i=1,nstates
      TPa(:,:,j)=TPa(:,:,j)+Coeff(i,j)*Coeff(i,j)*Pa(:,:,i)
      TPb(:,:,j)=TPb(:,:,j)+Coeff(i,j)*Coeff(i,j)*Pb(:,:,i)
    enddo
  enddo

  conf_index=0
  IF(ABS(iROOT-iSTATE)>0) THEN
     call CONF_PURE(Coeff(:,1:iROOT),nstates,iROOT,iSTATE,conf_index)
  ENDIF
  Pa=0.0D0
  Pb=0.0D0
  DO I2=1,iSTATE
    IF(ABS(iROOT-iSTATE)>0) THEN
       I1=conf_index(I2)
    ELSE
       I1=I2
    ENDIF
    Pa(:,:,1)=Pa(:,:,1)+TPa(:,:,I1)*WEIGHTS(I2)
    Pb(:,:,1)=Pb(:,:,1)+TPb(:,:,I1)*WEIGHTS(I2)
  ENDDO

  Energy_XC=0.0D0
  call EngineUp(ncenters,imult,icharge,functional_in,&
                    Pa,Pb,nconts,coord,atomchg,&
                    Energy_XC,1)
  print*,'NEW EXC IS:', Energy_XC(1)
else if(imethod==0) then
  call EngineUp(ncenters,imult,icharge,functional_in,&
                    Pa,Pb,nconts,coord,atomchg,&
                    Energy_XC,nstates)
endif
end 


!================================
subroutine EngineUp(ncenters,imult,icharge,functional_in,&
                    Pain,Pbin,ncontsin,coord,atomchg,&
                    Energy_XC,n_states)
! Main program of Engine
!  
!===============================

use GRID_INFO
use MOL_info
implicit none
real(8),parameter   :: ans2bohr = 1.889725989
real(8),parameter   :: PI = 3.14159265359
real(8),parameter   :: covrad(26) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,&
                                     1.36,1.25,1.27,1.39,1.25 ]
!include "parameter.h"
! input data zone
integer       :: ncenters,imult,icharge,n_states
integer       :: atomchg(ncenters)  ! natoms, multiplicity, charge
character     :: functional_in*30  ! functional
real(8)       :: coord(ncenters,3) ! coordinate, atomcharge
! output data zone
real(8)    :: E_ex,E_corr,Energy_XC(n_states)
integer    :: itmax
integer    :: ncontsin
real(8)    :: Pain(ncontsin,ncontsin,n_states)
real(8)    :: Pbin(ncontsin,ncontsin,n_states)
real(8),allocatable ::  Fxca(:,:),Fxcb(:,:)
integer       :: Gtype    ! guesstype, 1: Pseudo Huckel 2:Hcore
                    ! force,                Mullikencharge
integer :: i, j ,k,l
integer :: info

real(8)   :: dist
!###
real(8)   :: rho
!
!
!====

!============
!# 0. Preparing basic data

! prepare data
Natoms = ncenters
Multi = imult
Charge = icharge
Functional = functional_in
nconts = ncontsin
!nContssph = 0


!allocate(atoms(Natoms))

! read coordinate

coreChg = 0
do i = 1,Natoms
  atoms(i)%coor = coord(i,:)
  atoms(i)%charge = atomchg(i) 
  coreChg = atoms(i)%charge + coreChg 
enddo
print *,atoms(1)%base
n_alpha = 0.5*(Multi-Charge+coreChg-1)
n_beta = 0.5*(-Multi-Charge+coreChg+1)

print *,n_alpha,n_beta
! Build linkage matrix 

allocate(linkMat(natoms,natoms))
linkMat = 0
do i = 1,natoms
   do j = i+1,natoms
       dist = sum((atoms(i)%coor - atoms(j)%coor)**2)
       if (dist**0.5 < 1.2*(covrad(atoms(i)%charge) + covrad(atoms(j)%charge))) then 
         linkMat(i,j) = 1
         linkMat(j,i) = 1
       endif
   enddo
enddo



!# 1. Read basis set in and do basis set integral.
!call basshell(info)


allocate(Pa(nconts,nconts))
allocate(Pb(nconts,nconts))
allocate(Fxca(nconts,nconts))
allocate(Fxcb(nconts,nconts))


!# 2. Generate Grid and 
call gridgen(info)
call GTOeval(info)
print*,'2'

!# 3. Calculate E_xc
Energy_XC=0.0D0
do i=1,n_states
  print*,i,n_states
  Pa=Pain(:,:,i)
  Pb=Pbin(:,:,i)
  call  DFT_calc(nconts,E_ex,E_corr,Fxca,Fxcb,info)
  print *,E_ex,E_corr
  Energy_XC(i)=E_ex+E_corr
  Pain(:,:,i)=Fxca
  Pbin(:,:,i)=Fxcb
enddo
!! TEST Density !!
rho = 0
do i = 1,ngrids
   do j = 1,nconts
      do k =1,nconts
         rho = rho + (Pa(j,k) +Pb(j,k)) * grids(i)%val0(j)*grids(i)%val0(k)*grids(i)%weight
      enddo
   enddo
enddo

print *,rho
!!!!!!!!!!!!!!
deallocate(linkMat,Grids)
deallocate(Pa,Pb,Fxca,Fxcb)
!deallocate(atoms)
!call prtMat(Pa+Pb,nconts,"Density") 
end subroutine EngineUp
