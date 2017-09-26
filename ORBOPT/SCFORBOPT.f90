       SUBROUTINE SCF_ORBOPT(OneEI,TwoEI,n_orbital,Energy_nuclei, &
                             Transform_matrix,DFT_Fc)

        use global_control
        use matrix
        use date_time
        use orbopt_header
        logical ifdft 
        integer n_orbital
        real*8 Energy_nuclei
        real*8 DFT_Fc(n_orbital,n_orbital)
        real*8 OneEI(n_orbital,n_orbital)
        real*8 TwoEI(n_orbital,n_orbital,n_orbital,n_orbital)
        real*8 Transform_matrix(n_orbital,n_orbital)
        double precision,allocatable::TM1(:,:)
        double precision,allocatable::GM1(:,:,:,:)

        CALL FCIDUMP_READ()
        T=OneEI
        U=TwoEI
        T_origin=T
        U_origin=U
        T2H=T
        U2H=U
        FRONRE0=Energy_nuclei
        FRONRE=FRONRE0 
        call detect_system()
        call detect_redundancies()
        call HFenergy()
        call print_parameter(6)
        call DMRG_CASSCF()
        T_ini=T_origin
        U_ini=U_origin
 
        !write(6,*)"================================="
        !call print_mat(norb,norb,DFT_Fc,6)
        !stop

        if(DMRG_SCF)then
            i=1
            ifdft=.true.
            call RDMs_READ()
            allocate(TM1(nocc,nocc));           TM1=0.0d0
            allocate(GM1(nocc,nocc,nocc,nocc)); GM1=0.0d0 
            call closed_shell_update(mat2%d,mat2%P,TM1,GM1,1.0d0)  !CCCC
            deallocate(mat2%d,mat2%P)
            allocate(mat2%d(nocc,nocc))
            allocate(mat2%P(nocc,nocc,nocc,nocc))
            mat2%d=TM1   
            mat2%p=GM1 
            deallocate(TM1,GM1) 
            call PreOneSCF(i,DFT_Fc,ifdft)
            call Initial_value()
            call NonLinear_solver(i)
            if(LDIIS)then  
              call DIIS_acceleration(i,ITERVAL_DIIS)
            end if
            if(LDIIS_updated)then
              T=T_ini
              U=U_ini
              mat2%U=Udiis(i)%Ua 
              call print_mat(norb,norb,T,6)
              call print_mat(norb,norb,mat2%U,6)              
            end if
            call CS_dim_fro2cls()  ! CCCC
            call flush(6)
            Transform_matrix(1:norb,1:norb)=mat2%U(1:norb,1:norb)
            deallocate(T,T2H,T1dx,T2nd)
            if(allocated(Fc))then
              deallocate(Fc)
            end if
            deallocate(J1dx,K1dx)
            deallocate(J2dx,K2dx)
            deallocate(U,U2H,U2nd)
            deallocate(T_origin)
            deallocate(U_origin)
            deallocate(mat1)
            deallocate(mat2%T)
            deallocate(mat2%U)
            deallocate(mat2%R)
            deallocate(mat2%A)
            deallocate(mat2%B)
            if(allocated(mat2%Hdiag))then 
              deallocate(mat2%Hdiag)
            end if  
            deallocate(mat2%G)
            deallocate(mat2%H)
            deallocate(mat2%Rocc)
            deallocate(mat2%Aocc)
            deallocate(mat2%Hocc)
            deallocate(mat2%D)
            deallocate(mat2%P)
            deallocate(mat2%Ederi)
            deallocate(mat2%deltR)      
            deallocate(mat2%deltT)      
            deallocate(Umat_FINAL)      
            deallocate(mat2%A0)
            deallocate(mat2%B0)
            deallocate(redu%irreps)
            deallocate(redu%reflct)
            deallocate(redu%offset)
            deallocate(redu%orbirp)
            deallocate(redu%sign_mat)
            deallocate(threshold)
            deallocate(energies)
            deallocate(Rabs)
            deallocate(Udiis)
            deallocate(iorder)
            deallocate(T_ini,U_ini)
        end if

        if(DMRG_LAG)then          
            write(6,*)
            write(6,*)"--------------------------------------- "
            write(6,*)" Entering the linear-response solver "
            write(6,*)"--------------------------------------- "
            write(6,*)
            write(6,"(2X,A,I1,A)")"The ",dmrg_rlxstate,&
            "-th state will be relaxed"
            write(6,*)

            write(2,*)
            write(2,*)"--------------------------------------- "
            write(2,*)" Entering the linear-response solver "
            write(2,*)"--------------------------------------- "
            write(2,*)

            write(1,*)
            write(1,*)"--------------------------------------- "
            write(1,*)" Entering the linear-response solver "
            write(1,*)"--------------------------------------- "
            write(1,*)

            allocate(UMAT_FINAL(norb,norb)); UMAT_FINAL=0.0d0
            if(.NOT.ALLOCATED(mat2%U))then
              allocate(mat2%U(norb,norb))
              mat2%U=0.0d0
            end if
            do i=1,norb
              UMAT_FINAL(i,i)=1.0d0
                  mat2%U(i,i)=1.0d0
            end do 

!            call PreOneSCF(0)
            call FCIDUMP_READ_LR()             
            call run_DMRG_LR()
            call PreDMRGLR() 
            call ConDMRGLR()
            call LagrangeRDMs()

            ! MO Integrals [not relax in LR (only in CP-NR)]
            deallocate(T,U)

            ! Allocated in ConDMRGLR
            deallocate(mat2%D,  mat2%P)
            deallocate(mat2%A)
            deallocate(mat2%Hdiag)
            deallocate(mat2%G, mat2%H)
            ! 
            deallocate(T_origin)
            deallocate(U_origin)
            deallocate(UMAT_FINAL) 
        end if

      end

