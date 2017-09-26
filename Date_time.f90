module date_time

  implicit none

      integer, public :: systime0
      integer, public :: systime(100)
      real*8 , public :: cputime0
      real*8 , public :: walltime0
      real*8 , public :: walltime(100)
      real*8 , public :: cputime(100)
      real*8 , public :: time_start

end module date_time

real*8 function wtime ( )

!*****************************************************************************80
!
!! WTIME returns the wall clock time.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) WTIME, the wall clock in seconds.
!

  integer(kind = 4) :: clock_max
  integer(kind = 4) :: clock_rate
  integer(kind = 4) :: clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = 8 ) &
        / real ( clock_rate, kind = 8 )

  return
end

 subroutine timing_reporter(funit,string,timing)

 real*8          , intent(in) :: timing
 integer         , intent(in) :: funit
 character(len=*), intent(in) :: string

 write(funit,'(a,a,a,f12.2,a)') ' time needed for ',string,': ',timing,' seconds'

 end subroutine timing_reporter


