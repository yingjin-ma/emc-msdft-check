module global_control

! May affact some shell commands
      character*64 Linux

      character*64 FCIDUMP
! total molecular orbitals      
      integer norb
! total active orbitals
      integer nact,nact2,nocc
! If frozen core approximation
      logical ifcore,ifcore2
! orbitals in sub-groups
      type::orbitals
        character*3 sym 
        integer nsub
        integer,allocatable::grouptable(:,:)
        integer,allocatable::frozen(:)
        integer,allocatable::closed(:)
        integer,allocatable::   occ(:)
        integer,allocatable::   act(:)
        integer,allocatable::  act2(:)
        integer,allocatable:: extnl(:)
        integer,allocatable:: total(:)
        integer,allocatable::   ele(:)
      end type orbitals
      type(orbitals) :: orb,orbbk

! Redundant variables
      type::redundant_variables
        integer,allocatable::irreps(:)
        integer,allocatable::reflct(:)
        integer,allocatable::offset(:)
        integer,allocatable::orbirp(:)
        double precision,allocatable::sign_mat(:,:)
      end type redundant_variables
      type(redundant_variables)::redu

! DMRG control 
      character*36 dmrg_binary_name
      character*36 dmrg_input_name
      character*36 dmrg_input_deri
      character*36 dmrg_rdm1_name
      character*36 dmrg_rdm2_name
      character*36 dmrg_reorder_name
      character*36 dmrg_output_name
      logical ifmpirun
      character*2 nproc_in_mpirun 
     
      ! If state-averaged or coupled form 
      integer           dmrg_nstates,dmrg_rlxstate,nTDMs 
      double precision  dmrg_weight(99) 
      character*36      dmrg_inputs(99) 
      character*36      dmrg_inputs_TDM(99) 
      character*36      dmrg_outputs(99)

! DMRG-LR control (Coupled Pertubation)
      logical DMRG_SCF  
      logical DMRG_LAG
      logical  CKP_LR

! DMRG-CASSCF control
      integer ncycle
      logical act_omit,restart_maquis 
      logical old_code,debug,LDIIS,LMCSCF
      character*10 method(99)
      double precision,allocatable::threshold(:)
      double precision,allocatable::energies(:)
      double precision,allocatable::Rabs(:)
      integer MAX_iter,Nci_update,NH2_update
      double precision step_damping
      type::thresholds
        double precision E ! energy
        double precision D ! Davidson
        double precision R ! rotation
        double precision G ! gradients
        double precision C ! run MPS update when coupled
        double precision S ! symmetry restrict  
      end type thresholds
      type(thresholds) :: thrs

! DIIS acceleration
      type::DIIS
        double precision,allocatable::Ui(:,:,:)
        double precision,allocatable::Ua(:,:)
      end type DIIS
      type(DIIS),allocatable::Udiis(:) 
      Integer ITERVAL_DIIS
      Logical LDIIS_updated

! sub-symmetry
      type::symmetry
        character*3 point_group 
        integer norb
      end type symmetry
      type(symmetry),allocatable::sym(:) 

! integrals in MO basis
      integer n1e,n2e 
      type::integrals
        double precision dv
        integer i,j,k,l
        type(integrals),pointer::p
      end type integrals
      type(integrals),pointer::head 
      type(integrals),pointer::ptr 
      type(integrals),pointer::tail

! integrals matrix
      double precision HF_ENERGY,SA_RDMenergy
      double precision RDM_ENERGY,Echeck,E0_micro
      double precision FRONRE,FRONRE0,FRONREA 
      double precision,allocatable::T(:,:)
      double precision,allocatable::Fc(:,:)
      double precision,allocatable::U(:,:,:,:) 
      double precision,allocatable::T_origin(:,:)
      double precision,allocatable::U_origin(:,:,:,:) 
      ! save the initial integrals
      double precision,allocatable::T_ini(:,:)
      double precision,allocatable::U_ini(:,:,:,:) 

! 1-index transformed MO-Integrals
      double precision,allocatable::T1dx(:,:)
      double precision,allocatable::J1dx(:,:,:,:)
      double precision,allocatable::K1dx(:,:,:,:)
      double precision,allocatable::J2dx(:,:,:,:)
      double precision,allocatable::K2dx(:,:,:,:)

! 2-index transformed MO-Integrals (second order to T)
      !   full dim
      double precision,allocatable::T2H(:,:)
      double precision,allocatable::U2H(:,:,:,:) 
      ! active dim
      double precision,allocatable::T2nd(:,:)
      double precision,allocatable::U2nd(:,:,:,:)
      logical direct_dmrg

! Active the coupling by integrals/Hessian
      logical CP_integrals 
      integer ith_HESS,ith_INTE
   
! The final transform matrix (this is needed if not from AO)
      double precision,allocatable::UMAT_FINAL(:,:)   
 

end module global_control

