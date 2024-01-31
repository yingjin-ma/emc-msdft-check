module orbopt_version

! stefan: DMRG interface version variables

  implicit none

! ORBOPT version
#define ORBOPT_VERSION "0.2"
!
/* #undef ORBOPT_VERSION_MAJOR */
#define ORBOPT_VERSION_MINOR 2
#define ORBOPT_VERSION_BUILD be144a86e072ce61df994735c7079fa8a0263c3d (U_mult_1index)
#define ORBOPT_GIT_VERSION  "be144a86e072ce61df994735c7079fa8a0263c3d (U_mult_1index)"
!
! ORBOPT version (full string)
#define ORBOPT_VERSION_STRING "ORBOPT - version: 0.2"

character(len=200), public :: orbopt_v = ORBOPT_VERSION_STRING
character(len=200), public :: orbopt_g = ORBOPT_GIT_VERSION

end module orbopt_version
