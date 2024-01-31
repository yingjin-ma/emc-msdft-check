module coupling_terms

! orbitals in sub-groups
      type::couplings
        integer num             ! integer:: MPSci_num (position)
        integer list            ! integer:: MPSci_list for effective list  (position)
        integer ipos            ! integer:: MPSci_ipos for re-scaled position
        integer n               ! integer:: MPSci counter (only non zero counter)
        integer sign            ! integer:: Phase value
        double precision dv     ! double :: MPS_ci value
        double precision dvp    ! double :: pertubed MPS_ci value
        double precision LV     ! double :: Ci delta/lagrange value
        double precision dsL    ! double :: sign of overlap (<i|mps>)
        double precision dsR    ! double :: sign of overlap (<mps|i>)
        double precision grad   ! double :: gradient for each MPSci vector
        double precision gi     ! double :: Perturbed energy for MPSci vector (except 1)
        double precision,allocatable::A(:,:)  
        double precision,allocatable::AL(:,:),AR(:,:)  
        double precision,allocatable::D(:,:) 
        double precision,allocatable::P(:,:,:,:)  
        double precision,allocatable::DL(:,:),DR(:,:)
        double precision,allocatable::PL(:,:,:,:),PR(:,:,:,:)  
      end type couplings
      type(couplings),allocatable::Ci(:)
     
! dim of coefficient tensor (Davsion vector)
      integer nC,nC_all 
! Hessian (Ci-orb)
      double precision,allocatable::HCR(:,:,:)
      double precision,allocatable::KCR(:,:,:) 
! Hessian (Ci-Ci)
      double precision,allocatable::HCC(:,:)
! Delta Hessian
      double precision,allocatable::deltaH(:,:,:,:)
! Hamiltonian form DMRG
      double precision,allocatable::Hami(:,:),Hvec(:,:),Hval(:)
! Threshold for read in
      double precision,parameter::dthrs=1.0e-12

end module coupling_terms

