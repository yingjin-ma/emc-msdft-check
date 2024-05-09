module matrix

      use coupling_terms

      type::orb_opt1
        double precision,allocatable::R(:,:) ! rotation matrix
        double precision,allocatable::T(:,:) ! 1-integrals in MO basis  
        double precision,allocatable::U(:,:) ! 2-integrals in MO basis  
        double precision,allocatable::A(:,:)
        double precision,allocatable::B(:,:)
        double precision,allocatable::I(:,:) 

        double precision,allocatable::EDeri(:,:)
        double precision,allocatable::Hdiag(:,:)
        double precision,allocatable::H(:,:,:,:)

        double precision,allocatable::D(:,:)
        double precision,allocatable::P(:,:,:,:)

        double precision,allocatable::hD(:,:)
!        double precision,allocatable::JP_KL(:,:,:,:)

      end type orb_opt1
      type(orb_opt1),allocatable::MAT1(:) 
  
      type::orb_opt2
        double precision,allocatable::A(:,:)
        double precision,allocatable::B(:,:)
        double precision,allocatable::A0(:,:)
        double precision,allocatable::B0(:,:)
        double precision,allocatable::T(:,:)
        double precision,allocatable::U(:,:)
        double precision,allocatable::R(:,:)
        double precision,allocatable::deltR(:,:)
        double precision,allocatable::deltT(:,:)

        double precision,allocatable::D(:,:)
        double precision,allocatable::P(:,:,:,:)

        double precision,allocatable::Ederi(:,:)
        double precision,allocatable::H(:,:,:,:)
        double precision,allocatable::G(:,:,:,:)
        double precision,allocatable:: Y(:,:,:,:)
        double precision,allocatable::Y1(:,:)
        double precision,allocatable::Y2(:,:)
        double precision,allocatable::Dh(:,:,:,:)
        double precision,allocatable::Hdiag(:,:)

        double precision,allocatable::Rocc(:,:)
        double precision,allocatable::Aocc(:,:)
        double precision,allocatable::Hocc(:,:,:,:)

! The coupling term
        double precision,allocatable::HCR(:,:,:)
        double precision,allocatable::HCC(:,:)
      end type orb_opt2
      type(orb_opt2)::MAT2

      type::SCForLR_parameters
        ! Number of Davidson elements
        integer nC
        integer nC_list
        character*32,allocatable::position(:)
        integer,allocatable::pos(:,:)
        integer nC_all
        double precision RDM_energy 
        ! orbital Lagrange 
        double precision,allocatable::R(:,:)
        ! gradients contributions and gradients  
        double precision,allocatable::A(:,:)       
        double precision,allocatable::G(:,:)       
        ! density matrix
        double precision,allocatable::d(:,:)
        double precision,allocatable::P(:,:,:,:)
        ! extended "couplings" class from "use coupling_terms"
        ! it including all the Davidson elements related data
        type(couplings),allocatable::Ci(:) 
        ! Hessian (Ci-orb)
        double precision,allocatable::HCR(:,:,:)
        double precision,allocatable::KCR(:,:,:)
        ! Hessian (Ci-Ci)
        double precision,allocatable::HCC(:,:)
        ! Hamiltonian form DMRG
        double precision,allocatable::Hami(:,:),Hvec(:,:),Hval(:)
        ! Lagrange (orbital) corrected 1'- &' 2-RDMs          
        double precision,allocatable::Dorb(:,:) 
        double precision,allocatable::Porb(:,:,:,:) 
        ! Lagrange   (MPSci) corrected 1'- &' 2-RDMs          
        double precision,allocatable::Dmps(:,:) 
        double precision,allocatable::Pmps(:,:,:,:) 
      end type SCForLR_parameters
      type(SCForLR_parameters),allocatable::MPS(:)

      type::MPS_all
        integer,allocatable::ipos(:)
      end type MPS_all
      type(MPS_all),allocatable::mpsall(:)

      type::MPSlist
        character*32 string
        type (MPSlist), pointer :: p
      end type MPSlist
      ! type in matrix
      type (MPSlist), pointer :: mps_head
      type (MPSlist), pointer :: mps_ptr
      type (MPSlist), pointer :: mps_tail

      type::LRlist
        integer ipos
        type (LRlist), pointer :: p
      end type LRlist
      ! type in matrix
      type (LRlist), pointer :: lr_head
      type (LRlist), pointer :: lr_ptr
      type (LRlist), pointer :: lr_ptr1
      type (LRlist), pointer :: lr_ptr2
      type (LRlist), pointer :: lr_tail

      type mps_pos_info
        integer ipos
        integer  pos(3)
      end type  
      type(mps_pos_info),allocatable::mpspos(:)

      type::SCForLR_TDMs
        !gradients contributions 
        double precision,allocatable::A(:,:)       
        ! density matrix
        double precision,allocatable::d(:,:)
        double precision,allocatable::P(:,:,:,:)
      end type SCForLR_TDMs
      type(SCForLR_TDMs),allocatable::MPSTDMs(:,:) 

      type::LR_noredundant
        integer nC       
        integer,allocatable::idx(:)       
        integer,allocatable::pos(:)       
        ! Residual  (RHS)
        double precision,allocatable::Y(:) 
        ! Hessian   (LHS)
        double precision,allocatable::HOO(:,:) 
        double precision,allocatable::HOC(:,:) 
        double precision,allocatable::HCO(:,:) 
        double precision,allocatable::HCC(:,:) 
!       ! LHS & RHS & Lagrange
!        double precision,allocatable::LHS(:,:) 
!        double precision,allocatable::RHS(:) 
!        double precision,allocatable::LAG(:)
        double precision,allocatable::R(:,:) 
        ! Weight Lagrange 
        double precision,allocatable::Dorb(:,:) 
        double precision,allocatable::Dmps(:,:) 
        double precision,allocatable::Porb(:,:,:,:) 
        double precision,allocatable::Pmps(:,:,:,:) 
      end type LR_noredundant
      type(LR_noredundant) L_R
      double precision,allocatable::HLR(:,:)
      double precision,allocatable::MhD(:,:)
      integer,allocatable::iorder(:)
      !Dspmat_CSR
      integer,allocatable::rowoffset(:) 
      integer,allocatable::colindex(:) 
      double precision,allocatable::values(:)
      integer::nonzeros
      !Sspmat_CSR
      integer,allocatable::SHProwoffset(:) 
      integer,allocatable::SHPcolindex(:) 
      real(4),allocatable::SHPvalues(:)
      integer::SHPnonzeros
      !spmat_COO
      integer,allocatable::coo_rows(:)
      integer,allocatable::coo_cols(:)
      double precision,allocatable::coo_values(:) 

end module matrix

