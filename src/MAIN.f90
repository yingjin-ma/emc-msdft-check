      Program orbopt

!
        use global_control
        use matrix
        use date_time
        use orbopt_header

!"====================================================================="
!"       Orbital optimization program, version 0.2 (alpha)             "
!  Main references are
!        H.-J. Werner and W. Meyer,      J. Chem. Phys, 73, 2342 (1980)
!        H.-J. Werner and P. J. Knowles, J. Chem. Phys. 82, 5053 (1985)
!        P.G. Szalay, T. Muller, G. Gidofalvi, H. Lischka,
!        and R. Shepard,                    Chem. Rev. 112, 108  (2012)
!"---------------------------------------------------------------------"
!" Code written by                                                     "
!                              Yingjin Ma                              "
!"                                                                     "
!"====================================================================="
! Origional version should around 2009
!                         updated 2013 for combining Block
!                         updated 2015 for combining QCmaquis
!                         updated 2016 for Second-order SCF
!                         updated 2017 for Linear-response
! =====================================================================
! a prototype code for Second-order DMRG-SCF
!                    & DMRG SA-gradients (Hessian matrix part).

! Steps for this prototype code should be:
! 1) combine Maquis                                          (done)
! 2) refine the solver, i.e. write a new Hx=g solver         (done)
! 3) write the remaining part of Hessian                     (done)
! 4) import RDM derivatives                                  (done)
! 5) check if it is work                                     (checked)
! *  improve the integrals part                              (improved)
! 6) closed shell approximation                              (done, should with AH algorithm or MCLR)

! to Stefan : This part not yet done (code done, but have error), which will deduce the cost of transformation during micro-iteratons
! 7) 2-index trans. split into U * matrix, matrix from 1-index transformation. (to do)
! 8) Re-code the RDM-deri part (to do)
! 9) AVX improvement (to do)

! Current solvers for SCF
!   1) Werner-Meyer-Knowles with or without coupling
!   2) (Linear) Newton-Raphson
!   3) Augmented Hessian Newton-Raphson
!   4) Step-restricted Augmented Hessian
!   5) Werner-Meyer's finite difference
!   6) direct inversion
! PS: Only the WMK is well tested and designed
!     Davidson or Full diagionlization is used in 1-4

! Current solvers for LR
!   1) (preconditioner) conjugate gradient

! If it goes well, shift to MOLCAS

!        write(*,*)"starting..."
!       Used for closed/frozen exchanges
        double precision,allocatable::TM1(:,:)
        double precision,allocatable::GM1(:,:,:,:)

        !> 
        n = iargc()
        Cworkdir=""  
        DO i=1,n
          Cworkdir=""  
          CALL getarg(i, Cworkdir)
          WRITE(*,*)"The workdir is ",Cworkdir
          if(n.eq.1)Cworkdir1=CWorkdir
          if(n.eq.2)Cworkdir2=CWorkdir
        END DO
 
        !> program header
        call print_orbopt_header()

        call cpu_time(cputime(1))
        cputime0     = cputime(1)
        walltime(1)  = wtime()
        walltime0    = walltime(1)

        call input()
        !        write(*,*)"after input ..."
        call detect_system()
        call detect_redundancies()
        walltime(2)  = wtime()
        call timing_reporter(3,"redundancies check",walltime(2)-walltime(1))

        call HFenergy()

        walltime(3)  = wtime()
        call timing_reporter(3,"HF-energy check",walltime(3)-walltime(2))

        call print_parameter(6)
        call DMRG_CASSCF()
        !write(*,*)" After DMRG-CASSCF title "
        T_ini=T_origin
        U_ini=U_origin
        ! If allow the DMRG-SCF calculations

        walltime(4)  = wtime()
        call timing_reporter(3,"parameter initialization",walltime(4)-walltime(3))

        ! to Leon : IF DMRG-MCLR, use the convergenced w.f. for MOLCAS;
        !           In order to keep the consistency for MCLR part
        !           (i.e. don't use this DMRG-SCF utility)

        if(DMRG_SCF)then

          do i=1,ncycle

            walltime(10)  = wtime()
            if(i.eq.1)then
              call Run_DMRG(i,0)
            else
              call FCIDUMP_READ(i)              
            end if

            walltime(5)  = wtime()
            write(3,"(A,I3)") &
            " Macro-iter:",i
            write(30,"(A,I3)") &
            " Macro-iter:",i
            call timing_reporter(3,"MPS optimization",walltime(5)-walltime(10))
            call timing_reporter(30,"MPS optimization",walltime(5)-walltime(10))

            !1122      continue
            ! This should be removed later
            !          if(ifcore2.and.i.gt.1)then
            !            ifcore=.true.
            !            nact=nact2
            !            orb%act=orb%act2
            !          end if
                        ! Read in the active RDMs
            !            write(6,*)"before RDMs read in MAIN.f90 "; call flush(6)
            call RDMs_READ()

            !write(6,*)"After RDMs read in", nocc ; call flush(6)
            ! Update them to occupied RDMs
            allocate(TM1(nocc,nocc));           TM1=0.0d0
            allocate(GM1(nocc,nocc,nocc,nocc)); GM1=0.0d0
            !call print_mat(nact,nact,mat2%d,6)
            !call print_gat(nact,nact,nact,nact,mat2%p,6)
            !stop
            call closed_shell_update(mat2%d,mat2%P,TM1,GM1,1.0d0)  !CCCC

            deallocate(mat2%d,mat2%P)


            allocate(mat2%d(nocc,nocc))
            allocate(mat2%P(nocc,nocc,nocc,nocc))
            mat2%d=TM1
            mat2%p=GM1
            deallocate(TM1,GM1)

            walltime(11) = wtime()
            call timing_reporter(3,"closed-shell preparation",walltime(11)-walltime(5))

            write(15,*)"After closed shell update"; call flush(15)
            call print_mat(nocc,nocc,mat2%d,15)
            write(16,*)"After closed shell update"; call flush(16)
            call print_gat(nocc,nocc,nocc,nocc,mat2%p,16)

            !write(6,*)"After RDMs read in"; call flush(6)
            call PreOneSCF(i)

            walltime(12) = wtime()
            call timing_reporter&
            (3,"generation of orb-opt-related operators",walltime(12)-walltime(11))
            !           write(6,*)"After PreOneSCF"; call flush(6)
            !            stop
            call Initial_value()
            !           write(6,*)"After Initial_value "; call flush(6)
            call NonLinear_solver(i)
            !           write(6,*)"After Nonlinear solver "; call flush(6)
            !            call print_mat(norb,norb,mat2%U,6)
            walltime(3) = wtime()
            call timing_reporter(3,"non-linear solver",walltime(3)-walltime(12))

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

            walltime(4) = wtime()
            call timing_reporter(3,"closed-shell update",walltime(4)-walltime(3))

            call transform_INT(.true.)

            walltime(13) = wtime()
            call timing_reporter(3,"4-index transformation",walltime(13)-walltime(4))

            !call system("mv FCIDUMP_NEW_ACTIVE FCIDUMP_ACTIVE")
            call Scf_save(i,dmrg_binary_name) ! Also FCIDUMP and FCIDUMP_AVTIVE
            !            stop
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

            walltime(14) = wtime()
            call timing_reporter(3,"macro-iteration loop",walltime(14)-walltime(10))
            call timing_reporter(30,"macro-iteration loop",walltime(14)-walltime(10))

            if(Rabs(i).lt.thrs%r.and.Echeck.lt.thrs%e)then
              write(*,*)
              write(*,*)"Macro iterations converged!"
              write(*,*)"Orbital optimization has successfully finished"
              write(*,*)"----------------------------------------------"
              write(*,*)"Final ENERGY : ",RDM_ENERGY+FRONRE
              write(*,*)"Final sum|R| : ",Rabs(i)
              exit
            end if
            if(i.eq.ncycle)then
              write(*,*)"Maximum number of macro iterations reached!"
              write(*,*)"Convergence criteria were NOT satisfied!"
              write(*,*)"----------------------------------------"
              write(*,*)"Final ENERGY : ",RDM_ENERGY+FRONRE
              write(*,*)"Final sum|R| : ",Rabs(i)
              exit
            end if
            !if(method(i).eq."cp-hess")stop

          end do
        end if

        if(allocated(Fc))then
          deallocate(Fc)
        end if
        if(allocated(T2H))then
          deallocate(T2H)
        end if
        if(allocated(U2H))then
          deallocate(U2H)
        end if
        if(allocated(J1dx))then
          deallocate(J1dx)
        end if
        if(allocated(K1dx))then
          deallocate(K1dx)
        end if
        if(allocated(J2dx))then
          deallocate(J2dx)
        end if
        if(allocated(K2dx))then
          deallocate(K2dx)
        end if
        if(allocated(U2nd))then
          deallocate(U2nd)
        end if

        deallocate(T_origin)
        deallocate(U_origin)

        deallocate(T_ini,U_ini)
        deallocate(Udiis) ! allocated in DMRG_CASSCF()
!        write(6,*)DMRG_SCF,DMRG_LAG

        walltime(15) = wtime()
        write(6,"(A)")" ===== DMRG-SCF optimization finished ===== "
        write(6,"(A,f12.2)")" ===== time required: ",walltime(15)-walltime0
        call timing_reporter(3,"DMRG-SCF optimization",walltime(15)-walltime0)

        ! to Leon : Entering the linear-response solver

        ! If solving for Lagrange multipliers
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
            else
              mat2%U=0.0d0
            end if
            do i=1,norb
              UMAT_FINAL(i,i)=1.0d0
                  mat2%U(i,i)=1.0d0
            end do

            !            call PreOneSCF(0)
            ! to Leon : read in the integrals in MO basis
            call FCIDUMP_READ_LR()

            ! to Leon : obtain all RDMs and RDM-deris  (need an automatic bypass utility)
            !            call run_DMRG_LR()

            ! to Leon : Calculate all needed elements for constructing LR equations
            call PreDMRGLR()

            ! to Leon : Construct and solve the LR equation
            call ConDMRGLR()

            ! to Leon : Update the RDMs for a specific state in SA space
            call LagrangeRDMs()

            ! Allocated in ConDMRGLR
            deallocate(mat2%D,  mat2%P)
            deallocate(mat2%A)
            deallocate(mat2%Hdiag)
            deallocate(mat2%G, mat2%H)
            !
            deallocate(UMAT_FINAL)

        end if

        ! MO Integrals [not relax in LR (only in CP-NR)]
        deallocate(T,U)
        deallocate(T_origin)
        deallocate(U_origin)

        deallocate(redu%irreps)
        deallocate(redu%sign_mat)

        write(*,*)" "
        write(*,*)"----------- Normal Termination ------------ "
        write(*,*)" "

!       stop
!        write(*,*)"==========================================================="
!        write(*,*)"Next stage:"
!        write(*,*)"   Further usage of abelian group in transforming integrals" (done)
!        write(*,*)"   Stage-average treatment"                                  (done)
!        write(*,*)"   Coupled (one-step) treatment"                             (done)
!        write(*,*)"   Spin-orbit coupling"
!        write(*,*)"==========================================================="
!        write(*,*)""
!        write(*,*)"Thanks to my colleagues in Nanjing/ETHZ"
      end program orbopt

