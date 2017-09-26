      Subroutine Coupling_update_LR(iloop)

        use global_control
        use matrix
    
        integer::iloop 

        call system("rm -rf oneRDM.*.deri*")
        call system("rm -rf twoRDM.*.deri*")
        call RDMderi_update()        

!        write(1,*)"RDM_ENERGY, E0_micro", RDM_ENERGY, E0_micro
!        RDM_ENERGY=E0_micro
!        write(1,*)"RDM_ENERGY <-- E0_micro"

        call PreDMRGLR(iloop)

        stop 

      End Subroutine Coupling_update_LR



      Subroutine Coupling_update(iloop)

        use global_control
        use coupling_terms
        use matrix

        integer::iloop

        !call system("rm -rf info*")
        call system("rm -rf oneRDM.0.deri*")
        call system("rm -rf twoRDM.0.deri*")
        call RDMderi_update() 

        write(1,*)"RDM_ENERGY, E0_micro", RDM_ENERGY, E0_micro
        RDM_ENERGY=E0_micro
        write(1,*)"RDM_ENERGY <-- E0_micro"

        call coupling_prepare(iloop)
        call MPSci_gradient()
        call coupling_CIORB()
        call coupling_CICI()             
!        stop 

      End Subroutine Coupling_update

