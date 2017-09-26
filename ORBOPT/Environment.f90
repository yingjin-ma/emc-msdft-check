      Subroutine detect_system()

        use global_control 

        call system("cat /etc/issue > linux.version")
 
        linux=""
        open(unit=100,file="linux.version")
          read(100,*)linux
        close(100)

        if(linux(1:6).eq."Ubuntu")then
          linux="ubuntu"
        else
          linux="redhat"
        end if

!        write(6,*)"linux system is ",linux
        call flush(6)

      End Subroutine detect_system
