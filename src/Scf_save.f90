      Subroutine Scf_save(ncycle,dmrg_binary_name)

        integer::ncycle
        character*36::dmrg_binary_name 
        character*7 A1,step

        open(unit=202,file="SCF_SAVE")
          write(202,2021,advance='no')"mkdir Step_"
          if(ncycle.le.9)then
            write(202,"(I1,I1)")0,ncycle
          else
            write(202,"(I2)")ncycle
          end if 
          call system("cp SCF_SAVE SCF_tmp")
          open(unit=203,file='SCF_tmp')
            read(203,*)A1,step
          close(203)
          call system("rm SCF_tmp")
          write(202,2022,advance='no')"mv FCIDUMP " 
          write(202,2030)step
          write(202,2023,advance='no')"mv FCIDUMP_ACTIVE "
          write(202,2030)step
          if(dmrg_binary_name.ne."maquis")then
            write(202,2024,advance='no')"mv spatial_* "
            write(202,2030)step
          else
!            write(202,2025,advance='no')"mv *RDM.*.* " 
!            write(202,2030)step
          end if
!          write(202,2024,advance='no')"mv EU.txt    " 
!          write(202,2030)step
        close(202)

        call system("chmod +x SCF_SAVE")
        call system("./SCF_SAVE")
        call system("mv FCIDUMP_NEW FCIDUMP")
        call system("mv FCIDUMP_NEW_ACTIVE FCIDUMP_ACTIVE")

2021    format(A11)
2022    format(A11)
2023    format(A18)
2024    format(A13)
2025    format(A17)
2030    format(A7)

      end Subroutine Scf_save

