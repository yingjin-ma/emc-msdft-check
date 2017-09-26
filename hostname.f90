module hostname

  use iso_c_binding

  implicit none

  public print_host

  interface
    integer( kind = C_INT ) function gethostname( hname, len ) bind( C, name = 'gethostname' )
        use iso_c_binding
        implicit none
        character( kind = C_CHAR ) :: hname( * )
        integer( kind = C_INT ), VALUE :: len
    end function gethostname
  end interface

contains

!-------------------------------------------------------------------------------
  subroutine print_host(luin)

    integer, intent(in)                :: luin
!-------------------------------------------------------------------------------
    integer( kind = C_INT ), parameter :: sl = 200
    character( kind = C_CHAR )         :: hn( sl )
    character( len = sl )              :: fn
    character                          :: c
    integer                            :: res, i, j
!-------------------------------------------------------------------------------

    res = gethostname( hn, sl )
    if ( res == 0 ) then 
        do i = 1, sl
            c = hn( i )
            if ( c == char( 0 ) ) exit
            fn( i: i ) = c
        end do
        do j = i, sl
            fn( j: j ) = ' '
        end do
        write(luin,'(/a,a/)') " Calculation is running on --> ", trim( fn )
    else
        stop 'hostname not available'
    end if
  end subroutine print_host

end module hostname
