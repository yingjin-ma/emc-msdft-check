      Subroutine Print_parameter(io_num)

        use global_control

        if(method(1).eq."WMKUBAR")then  
          write(io_num,*)" "
          write(io_num,*)"Thresholds (macro-iter.) : "
          write(io_num,"(A38,E8.2)")&
          "                             Energy : ",THRS%E
          write(io_num,"(A38,E8.2)")&
          "                 Rotation(sum|R|/2) : ",THRS%R
          write(io_num,*)"Thresholds (micro-iter.) initial : "
          write(io_num,"(A38,E8.2)")&
          "   Residual between E(val) & E(SCF) : ",THRS%C
!          "   Residual between E(T,c) & E(T,dR): ",THRS%C
          write(io_num,"(A38,E8.2)")&
          "   Residual for max gradient element: ",THRS%G
          write(io_num,*)
          write(io_num,*)&
          "! THRS (micro) will be adjusted automatically !"
          write(io_num,*)
        end if 

        if(dmrg_nstates.gt.1)then
          write(io_num,*)" "
          write(io_num,*)" State-averaged DMRG-SCF with weights: "    
          do i=1,dmrg_nstates
            write(io_num,"(1X,f8.3)",advance='no')dmrg_weight(i)
          end do
          write(io_num,*)" "
        end if 

      End Subroutine Print_parameter

