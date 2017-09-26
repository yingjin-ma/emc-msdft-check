module MOL_info
implicit none
type shells
    integer    :: angMoment
    integer    :: nGauss
    real(8),allocatable    :: exponents(:)
    real(8),allocatable    :: contrCoeff(:)
end type shells


type atom
    integer    :: charge
    real(8)    :: coor(3)
    character  :: base*30
    integer    :: nsShell
    integer    :: npShell
    integer    :: ndShell
    integer    :: nfShell
    integer    :: ngShell
    integer    :: nShell
    integer    :: nconts
    type(shells),allocatable :: shell(:)
end type atom
type(atom),allocatable :: atoms(:) 

!   Number of atoms  Number of contractions   Number of  contractions(spherical)
!         V                   V                  V 
integer Natoms,              Nconts ,            ncontssph

!===================
!     Martrix      +
!===================
real(8),allocatable  :: Hcore(:,:),S(:,:)
real(8),allocatable  :: TWOEI(:)
!Basis Normalization
!
real(8),allocatable  :: NorVEC(:)
real(8),allocatable  :: NorVECsph(:)
!                        Density_a   Density_b
!                        V           V
real(8),allocatable  ::  Pa(:,:),    Pb(:,:)

!===================
!    Properties    +
!===================
integer              :: Charge        ! Total charge
integer              :: Multi         ! Multiplicity

integer              :: n_alpha       ! alpha electrons
integer              :: n_beta        ! beta electrons
integer              :: coreChg       ! total core charges

integer,allocatable  :: linkMat(:,:)  ! Linkage matrix

!===================
!   Functional     +
!===================
character Functional*30
end module MOL_info


module GRID_info
!===================
!   Grid           +
!===================
implicit none
type Grid
    real(8)    :: coor(3)
    real(8)    :: weight
    real(8)    :: rho_a,rho_b
    real(8)    :: sigma_aa,sigma_ab,sigma_bb
    real(8),allocatable ::val0(:)
    real(8),allocatable ::val1(:,:)
end type Grid

integer :: ngrids

type(Grid),allocatable :: Grids(:) 
end module GRID_info



