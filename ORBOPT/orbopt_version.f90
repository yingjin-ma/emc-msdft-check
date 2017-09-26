module orbopt_version

! stefan: DMRG interface version variables

  implicit none

! ORBOPT version
#define ORBOPT_VERSION "0.2"
!
/* #undef ORBOPT_VERSION_MAJOR */
#define ORBOPT_VERSION_MINOR 2
#define ORBOPT_VERSION_BUILD a4f63f6d7524fe573ac33d00cde2a180048ca3c6 (zxqu_msdft)
#define ORBOPT_GIT_VERSION  "a4f63f6d7524fe573ac33d00cde2a180048ca3c6 (zxqu_msdft)"
!
! ORBOPT version (full string)
#define ORBOPT_VERSION_STRING "ORBOPT - version: 0.2"

character(len=200), public :: orbopt_v = ORBOPT_VERSION_STRING
character(len=200), public :: orbopt_g = ORBOPT_GIT_VERSION

end module orbopt_version
