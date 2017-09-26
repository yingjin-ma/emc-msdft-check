      Subroutine group_table()

        use global_control

        double precision d2h(8,8)

        d2h(1,:)= (/ 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 /) !Ag
        d2h(2,:)= (/ 2 , 1 , 4 , 3 , 6 , 5 , 8 , 7 /) !B1g
        d2h(3,:)= (/ 3 , 4 , 1 , 2 , 7 , 8 , 5 , 6 /) !B2g
        d2h(4,:)= (/ 4 , 3 , 2 , 1 , 8 , 7 , 6 , 5 /) !B3g
        d2h(5,:)= (/ 5 , 6 , 7 , 8 , 1 , 2 , 3 , 4 /) !Au
        d2h(6,:)= (/ 6 , 5 , 8 , 7 , 2 , 1 , 4 , 3 /) !B1u
        d2h(7,:)= (/ 7 , 8 , 5 , 6 , 3 , 4 , 1 , 2 /) !B2u
        d2h(8,:)= (/ 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1 /) !B3u
 
!        if(orb%sym.eq."d2h")then
          orb%grouptable=d2h
!        end if 

!        write(*,*)orb%grouptable

      end subroutine group_table
