
      Subroutine Solver_WMK2(icycle)

        use global_control
        use matrix
        use date_time  

        integer::icycle
        logical ciupdate,H2update,Lexit,Lcoupled,LNOcoupled 

        integer,allocatable::porb(:)
        double precision d0,d0_sum,d0_sum0,d2,d2_sum0,d2_sum
        double precision dv,redundant,THdeltaR,ET,ETdR
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:)
        double precision,allocatable::YQ(:)    ! check residual
        double precision,allocatable::G(:),GP(:)  ! The A-A^+, orbital gradient
        double precision,allocatable::H(:,:),HP(:,:),HPA(:,:),HPinv(:,:)
        double precision,allocatable::AP(:,:),BP(:,:)

        double precision,allocatable::T1(:,:),T2(:,:)

! Used for the augmented Hessian part
        double precision lamda               ! damping parameter that control the level shift
        double precision damping             ! in the loop
        double precision D_epsilon,D_lamda   ! epsilon 
        double precision grad_digit,XP1_P    ! in the loop
        double precision Econv,EconvCP0,EconvCP1  ! Check the convergence

! Save the initial A and B
        double precision Aini(norb,norb)
        double precision Bini(norb,norb)
! Gradients(orb) in Matrix form 
        double precision ResidM(norb,norb)
        double precision,allocatable::Gini(:,:,:,:)
!        double precision,allocatable::UBBUP(:,:),AAP(:,:)
        double precision TM1(norb,norb), TM2(norb,norb)
        double precision Ubk(norb,norb),Tbk(norb,norb)
        double precision AA(norb,norb),dRA(norb,norb),AdR(norb,norb)
!Save the valid rotation index
        integer,allocatable::valid(:,:)
        double precision,allocatable::wB(:),wBT(:),ARRA(:),SFT(:)
        logical Lvalid,Lcycle
        logical Lfinish

!       scratch for timing
        character(len=200) :: microlabel

        walltime(25) = wtime()
! Parameter initialization
        lamda=1.0d0 ! "1" means the standard AH 
        ! The larger lamda, the shorter step ??
        damping=step_damping  ! restrict initial

        Lfinish=.false.
        Lcycle=.false.

! initialize of part of Hessian matrix that will be used in this solver
        write(1,*)"entering the WMK solver"
        call flush(1) 

        ndim1=norb

        allocate(X(ndim1**2)); X=0.0d0
        allocate(Y(ndim1**2)); Y=0.0d0
        allocate(G(ndim1**2)); G=0.0d0
        allocate(H(ndim1**2,ndim1**2)); H=0.0d0

!        write(1,*)"after some allocate"
!        call flush(1) 
! Get the initial A, B matrix 
!       call derivaties(0)
!        Aini=mat2%A
!       Bini=mat2%B
! The initial G matrix
        allocate(Gini(norb,norb,norb,norb)); Gini=0.0d0
        Gini=mat2%G

!        write(1,*)"allocated G"
!        call flush(1)
 
        ! This is the orbiatl-rotation gradients 
        redundant=1.0d-9

        if(.false.)then 
          grad_digit=0.0d0
          ij=0
          nij=0
          do i=1,ndim1
            do j=i+1,ndim1 ! only half is needed
            !do j=i,ndim1 ! only half is needed
               ij=ij+1
               ! Residual, only A=UB , same as NR, AH if U=1
               Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
               G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
           
!               write(1,*)"checking",i,j 
!               call flush(1)

               if(dabs(Y(ij)).lt.redundant)then
               else
!                 write(1,*)"in else"
                 IF(LMCSCF)then
                   ! Check for degeneracy irreps
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                     write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                     nij=nij+1
                     if(dabs(G(ij)).gt.grad_digit)then
                       grad_digit=dabs(G(ij))
                     end if
                   end if
                   ! end Check for degeneracy irreps
                 else
                   ! no closed shell 
                   if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
!                     write(1,*)"in else : casscf type"                
!                     call flush(1)
                     ! Check for degeneracy irreps
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then      
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if 
                       end if
                     ! end Check for degeneracy irreps
                     end if
                   ! i is closed shell
                   else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if
                       end if
                     !end if
                   ! j is closed shell 
                   else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if
                       end if  
                     !end if  
                   end if ! close shell
                 end if ! MCSCF or CASSCF
               end if
            end do
          end do

          do i=1,ndim1
            do j=1,i ! only half is needed
            !do j=i,ndim1 ! only half is needed
               ij=ij+1
               ! Residual, only A=UB , same as NR, AH if U=1
               Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
               G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
               if(dabs(Y(ij)).lt.redundant)then
               else
                 IF(LMCSCF)then
                   ! Check for degeneracy irreps
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                     write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                     nij=nij+1
                     if(dabs(G(ij)).gt.grad_digit)then
                       grad_digit=dabs(G(ij))
                     end if
                   end if
                   ! end Check for degeneracy irreps
                 else
                   ! no closed shell 
                   if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       ! Check for degeneracy irreps
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if
                       end if
                       ! end Check for degeneracy irreps
                     end if
                   ! i is closed shell
                   else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if
                       end if
                     !end if
                   ! j is closed shell 
                   else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                         nij=nij+1
                         if(dabs(G(ij)).gt.grad_digit)then
                           grad_digit=dabs(G(ij))
                         end if
                       end if  
                     !end if
                   end if ! check closed shell
                 end if ! MCSCF or CASSCF
               end if             
            end do
          end do

        else
          call valid_rots_count&
              (icycle,norb,mat2%A,orb%nsub,redu%sign_mat,redu%orbirp,1.0d-9,Lmcscf,nij)
          write(6,*)"The effective rotation number : ", nij
        end if     

        nij0=nij
! re-run it to get the valid rotations index
        allocate(valid(nij,2)); valid=0

        if(.false.)then
          ij=0
          kl=0
            !do j=i,ndim1  ! only half is needed
          do i=1,ndim1
            do j=i+1,ndim1  !half
               ij=ij+1
               if(dabs(Y(ij)).lt.redundant)then
               else
                 IF(LMCSCF)then
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then       
                     kl=kl+1
                     valid(kl,1)=i
                     valid(kl,2)=j
                   end if
                 else
                   ! no closed shell 
                   if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then       
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                     end if
                   ! i is closed shell | i-a,i-v
                   else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                     !end if
                   ! j is closed shell | j-a,j-v 
                   else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                     !end if
                   end if
                 end if
               end if
            end do
          end do
          do i=1,ndim1
            do j=1,i  !half
               ij=ij+1
               if(dabs(Y(ij)).lt.redundant)then
               else
                 IF(LMCSCF)then
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then       
                     kl=kl+1
                     valid(kl,1)=i
                     valid(kl,2)=j
                   end if
                 else
                   ! no closed shell 
                   if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                       if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                         kl=kl+1
                         valid(kl,1)=i
                         valid(kl,2)=j
                       end if
                     end if
                   ! i is closed shell | i-a,i-v
                   else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                     !end if
                   ! j is closed shell | j-a,j-v 
                   else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                     !if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                     if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                       kl=kl+1
                       valid(kl,1)=i
                       valid(kl,2)=j
                     end if
                     !end if
                   end if ! 
                 end if
               end if
            end do
          end do
        else
          call valid_rots&
          (norb,mat2%A,orb%nsub,redu%sign_mat,redu%orbirp,1.0d-9,Lmcscf,nij,valid)
        end if

          write(1,*)"valid rotations"
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,1)
          end do
          write(1,*)" "         
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,2)
          end do
          call flush(1) 

          if(nij.eq.0)then
            write(*,*)"No orbital gradients large than ", redundant
            write(*,*)"Convergenced orbitals, no need orb-opt"
            goto 9999
          end if

!        print *," "
         write(2,*)&
         " ===================================================== "
         write(2,*)&
         " Start the micro-iters in the ",icycle,"-th macro-iter"
!        stop  

        write(2,*)"Echeck",Echeck
        write(2,*)"Damping",damping

!        if(icycle.eq.1)then
!          thrs_c=min(Echeck/100d0,thrs%c/(10**(icycle-1))) ! assuming each micro can reduce 1000
!          thrs_c=max(thrs%e/2.0d0,thrs_c)
!        else 
!          thrs_c=thrs%e/2.0d0
!        end if
!        thrs_c=thrs%c
!        write(2,*)"fix the coupling precision as that in input"

!        thrs_g=max(grad_digit/100.0,thrs%g)
        write(2,*)" thrs_coupling  ", thrs%c
        write(2,*)" thrs_energy    ", thrs%e
        write(2,*)" thrs_davidson  ", thrs%d
        write(2,*)" thrs_rotation  ", thrs%r
        call flush(1) 
        call flush(2) 

        d0_sum=0.0d0 
        d2_sum=0.0d0 
        d0_sum0=0.0d0 
        d2_sum0=0.0d0 
        iloop=0
        XP1_P=0.0d0
        iloopmax=max_iter 
        LNOcoupled=.true.
        Econv=RDM_energy*2.0d0 
        EconvCP0=0.0d0
        EconvCP1=Econv
        Nci_update00=0
        Nh2_update00=0
        ith_INTE_BK=ith_INTE 
        do 
          iloop=iloop+1

!          write(6,*)"in loop",iloop," of Solver_WMK, U-mat "
!          call flush(6)
!          call print_mat(norb,norb,mat2%U,6)

          walltime(18) = wtime()

          if(.true.)then ! Standard WMK algorithm  

! decide whether active coulped optimization
            ciupdate=.false.
            H2update=.true.
            !H2update=.false. 
            if(dabs(Econv-E0_micro).lt.thrs%c)then ! Eval and ESCF
              write(2,*)"Econv-E0_micro smaller than thrs_coupling"
              if((iloop-Nci_update00).lt.Nci_update)then
              write(2,*)" checked the interval "
                if(dabs(Econv-E0_micro).lt.thrs%e/10.0d0)then 
                  write(2,*)"active coupling update"
                  ciupdate=.true.
                end if
              end if  
            end if

            if(ciupdate.and.h2update)then
              Lcoupled=.true.
            else
              Lcoupled=.false.
            end if
 
            ! In case AH has problem
            if(CP_integrals.and.Lcycle)then
              Lcoupled=.true.
            end if
            if(.not.CP_integrals.and.Lcycle)then
              Lcycle=.false. 
              exit 
            end if

            write(2,*)iloop,Nci_update00
            if((iloop-Nci_update00).gt.Nci_update)then
              Lcoupled=.true.              
              write(2,*)"Force update MPS-CI after ", &
                        Nci_update,"micro-iters"
              Nci_update00=iloop
            end if

            write(2,*)iloop,Nh2_update00,ith_INTE
            if((iloop-Nh2_update00).gt.Nh2_update.and..not.Lcoupled)then
              Lcoupled=.true.             
              ith_INTE=0     ! control whether update MPS/CI
              write(2,*)"update H(2) after ",Nh2_update,"micro-iters"
            end if

            if(Lcoupled)then
              Nh2_update00=iloop
              XP1_P=0.0d0
              ciupdate=.true.
              H2update=.true. ! Should be updated if use the full Hessian in AH
              LNOcoupled=.false.
              if(ith_INTE.ne.0)then
                Nci_update00=iloop
              end if 
            else
              ciupdate=.false.
              H2update=.false. 
            end if
!            ciupdate=.false.
!            H2update=.false. 

            mat2%deltR=0.0d0
!            print *,"deltaR initial as 0 each loop, nact",nact
            call Operator_update(icycle,ciupdate,H2update,iloop) 
!            print *,"After operator update,nact",nact

            if(ciupdate.or.H2update)then
! to Stefan : Eq.32
!             updated in "Operator_update"
              Gini=mat2%G 
              write(2,*) &
              "updated G-operator is used from new RDMs in loop", iloop 
              ! no need the following line since U=1 in this case
              ! call G_update(Gini) 
            else 
              call G_update(Gini) 
            end if
            if(ith_INTE.ne.0)then
              EconvCP0=EconvCP1
              EconvCP1=E0_micro
            end if 

!            if(Lcoupled)then  ! In case changed
            ith_INTE=ith_INTE_BK
!            end if

!            if(iloop.eq.35)stop 
!           Bini=mat2%B

! =====================================================================
            TM1=0.0d0 
            TM1=mat2%A0+mat2%B0
            call MdXM(norb,mat2%T,TM1,TM2)
            dv=0.0d0
            call trace(norb,TM2,dv)
            dv=dv*0.5d0
            write(2,*)"The Tr[T^+(A+B)]/variation(delta_e2) value : ",dv
            call flush(2) 
            ET=dv+E0_micro+FRONRE0
            ETdR=ET
            write(2,*)"The 2nd energy expression E2(T) at T point : ",ET
            call flush(2) 
            
! =====================================================================
            !print *," === mat2%A === "
            !call print_mat(norb,norb,mat2%A,6)

          else
            !write(*,*)"in loop",iloop," of Solver_WMK "
! Update the B matrix
            call derivatives_B_update(Aini)
            Bini=mat2%B
! Get the ~A , ~B, now they are mat2%A, mat2%B
            call Derivatives_WMK(icycle,iloop,Gini,Bini)
            !write(6,*)"After derivative_wmk" 
            !call flush(6)
          end if
!          call print_mat(norb,norb,mat2%A,6)
!          call print_mat(norb,norb,mat2%B,6)
!          if(iloop.eq.5)stop
         
! Get the initial Hessian for deltaR, same as R if initial
          call Hessian(ndim1,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
! Get UB, BU, (A+A)dR, dR(A+A)
!          write(6,*)"After Hessian " 
          !call print_gat(ndim1,ndim1,ndim1,ndim1,mat2%G,6)
          !call flush(6)

          ij=0
          do i=1,ndim1
            do j=1,ndim1
              ij=ij+1
              kl=0
              do k=1,ndim1
                do l=1,ndim1
                  kl=kl+1
                  H(ij,kl)=mat2%H(i,j,k,l)
                end do
              end do
            end do
          end do

! Detect the valid orbital rotation

! =============================================================================
!  omit the act-act if CAS-type rotation, if "ture", all rotations are allowed
! =============================================================================
! to Stefan : active-active rotations

          if(LMCSCF)then 

            grad_digit=0.0d0
            ij=0
            nij=0
            do i=1,ndim1
              do j=i+1,ndim1 ! only half is needed
              !do j=i,ndim1 ! only half is needed
                 ij=ij+1
                 ! Residual, only A=UB , same as NR, AH if U=1
                 Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   ! Check for degeneracy irreps
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                     write(1,*)icycle,"-th macro-loop, ini-grad",i,j,Y(ij)
                     nij=nij+1
                     if(dabs(G(ij)).gt.grad_digit)then
                       grad_digit=dabs(G(ij))
                     end if
                   end if
                 end if
              end do
            end do
            do i=1,ndim1
              do j=1,i ! only half is needed
              !do j=i,ndim1 ! only half is needed
                 ij=ij+1
                 ! Residual, only A=UB , same as NR, AH if U=1
                 Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 G(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
!                   write(1,*)iloop,"-th loop, grad",i,j,Y(ij)
                     nij=nij+1
                     if(dabs(G(ij)).gt.grad_digit)then
                       grad_digit=dabs(G(ij))
                     end if
                   end if
                 end if
              end do
            end do


            nij0=nij
! re-run it to get the valid rotations index
            if(.not.allocated(valid))then
              allocate(valid(nij,2)); valid=0
            end if
            ij=0
            kl=0
              !do j=i,ndim1  ! only half is needed
            do i=1,ndim1
              do j=i+1,ndim1  !half
                 ij=ij+1
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                     kl=kl+1
                     valid(kl,1)=i
                     valid(kl,2)=j
                   end if
                 end if
              end do
            end do
            do i=1,ndim1
              do j=1,i  !half
                 ij=ij+1
                 if(dabs(Y(ij)).lt.redundant)then
                 else
                   if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                     kl=kl+1
                     valid(kl,1)=i
                     valid(kl,2)=j
                   end if
                 end if
              end do
            end do

            write(1,*)&
            "valid rotations re-check in micro-iters (MCSCF-type)"
            do i=1,nij
              write(1,"(I3)",advance='no')valid(i,1)
            end do
            write(1,*)" "
            do i=1,nij
              write(1,"(I3)",advance='no')valid(i,2)
            end do

          else       !  if  MCSCF
!            nij=nij0
!            allocate(valid(nij,2)); valid=0
 
          end if     !  if CASSCF

!          write(6,*)"The valid index in loop",iloop
!          write(6,*)valid(:,1)
!          write(6,*)valid(:,2)
          !if(iloop.eq.2)stop 
   
          nij=nij0/2
          !nij=nij0
          allocate(porb(nij));   porb=0
          allocate(YP(nij));     YP=0.0d0
          allocate(YQ(nij));     YQ=0.0d0
          allocate(GP(nij));     GP=0.0d0
          allocate(HP(nij,nij)); HP=0.0d0 
        
          ! This is the delta-orbiatl-rotation gradients 
          ! notice G!=Y
          ij=0
          kl=0
          do i=1,ndim1
            do j=1,ndim1 ! only half is needed
            !do j=1,ndim1 ! all
               ij=ij+1
               G(ij)=(mat2%A(i,j)-mat2%A(j,i)) ! gradients
               !Check the valid rotations  
               Lvalid=.false.
               do k=1,nij
                 !write(*,*)"k",k                 
                 if(i.eq.valid(k,1).and.j.eq.valid(k,2))then
                   Lvalid=.true.
                   exit
                 end if
               end do
               if(Lvalid)then  ! can be deleted
                 kl=kl+1
                 porb(kl)=ij
                 !write(*,*)"ij",ij,"kl",kl
               end if
            end do
          end do

          grad_digit=0.0d0 
          do i=1,nij
            GP(i)= G(porb(i))
            YP(i)= G(porb(i))  
            if(dabs(YP(i)).gt.grad_digit)then  ! ??
              grad_digit=dabs(YP(i))
            end if
            do j=1,nij    ! needed Hessian
              HP(i,j)=H(porb(i),porb(j))
            end do
          end do       

          write(1,*)"The valid orb-orb Hessian: "
          call print_mat(nij,nij,HP,1)

           dv=0.0d0 
           do i=1,nij
             dv=dv+dabs(GP(i))/nij
           end do
           write(2,*)" lagest gradient remaining ",grad_digit
           write(2,*)"Average gradient remaining ",dv
           call flush(2) 

          if(.true.)then  ! Here use modified Davidson /AH              

            do iAH=1,1 ! should be several AH steps in Davidson ??

! The new augmented Hessian matrix {except HPA(1,1)}
              allocate(HPA(nij+1,nij+1)); HPA=0.0d0
              HPA(2:nij+1,1)=GP
              HPA(1,2:nij+1)=GP
              HPA(2:nij+1,2:nij+1)=HP
! And the augmented x-vector
              allocate(XP(nij+1))
              allocate(D(nij+1))
              XP=0.0d0
               D=0.0d0

            !write(6,*)XP(1),"initial XP(1) in loop -",iloop
            !call print_mat(nij+1,nij+1,HPA,6)

!============================== Davidson WMK ====================================
              if(.false.)then  ! 
              call Davidson_WMK(nij+1,HPA,XP,D,lamda,valid(1:nij,:),Gini,Bini,10,1.0d-8)
              else
               call Davidson(nij+1,HPA,XP,D,lamda,.true.,nij+1,thrs%D)
              end if
!              call EIGTQL2(nij+1,HPA,D)
!              write(2,*)"From full diagonalization (eig)"
!              call print_mat(1,nij+1,D,2)
!              write(2,*)"From full diagonalization (vec)"
!              call print_mat(nij+1,nij+1,HPA,2)

! ===============================================================================
              write(2,*)D(1),"eig(1), and eig(2:nij) in loop ",iloop
              if(D(1).gt.0.0d0)then
            write(2,*)"================================================"
            write(2,*)"Warning : no negative eigenvalue from Davidson ?"
            write(2,*)"      consider More tight threshold for Davidson"
            write(2,*)"================================================"
              end if
              !call print_mat(1,nij,D(2:nij+1),2)          
              write(2,*)XP(1),"XP(1), and XP(2:nij) "
              call print_mat(1,nij,XP(2:nij+1),2)
              call flush(2) 
              ! check the XP flag, assuming energy should always decreases (delete)
              ! better solution should be more modification in AH/Davidson to shift vec
              ! Already used the shifted in AH
              Lcycle=.false.
              dv=dabs(XP(1))
              if(XP(1).gt.1.0d0)then
!               write(2,*)&
!               "Warning : AH not positive define even with shift"
!               write(2,*)&
!               "       : set rotation matrix as 0, cycle this loop"
!               Lcycle=.true.
!                XP(1)=1.0d0
!               ! random huess for the next 
!               call random_seed()
!               do i=1,nij
!                 call random_number(XP(i+1))
!                 XP(i+1)=XP(i+1)/100000.0d0
!               end do
              end if
              XP1_P=XP(1)
 
              if(1.0d0-dabs(XP(1)).lt.1.0d-15)then
                h2update=.true.
                Lcycle=.true.
                write(2,*)"\lamda close to 1.0d0"
              end if

              call inner_product(nij,GP,XP(2:nij+1),dv)
!              write(*,*)"dv",dv 
!              write(*,*)"epsilon should equ to eigen ",dv/xp(1) 
                D_lamda=1.0d0/xp(1)
              D_epsilon=dv/xp(1)
 
              if(XP(1).lt.0)then
                XP=-1.0d0*XP
!                XP=XP/XP(1)
              end if    
!              stop

              if(dabs(XP(1)).lt.damping)then
!                write(2,*)"Warning: AH(0) element is",XP(1)
!                write(2,*)"      better re-consider the active orbital"  
              end if
!              if(iloop.ge.0)then ! restrict the step length
!                do
!                  call MODIFIED_GRAM_SCHMIDT(nij+1,1,XP)
!                  if(dabs(XP(1)).lt.damping)then
!                    XP(1)=XP(1)*1.25d0
!                    XP(2:nij+1)=XP(2:nij+1)/1.25d0
!                  else
!                   !damping=dabs(XP(1)) ! reset damping
!                    exit
!                  end if
!                end do
!                write(2,*)"damping is",damping
!              end if

 
              do i=1,nij
                X(porb(i))=XP(i+1)
              end do

              write(1,*)

              ij=0; 
              d0=0.0d0; d2=0.0d0; 
              d0_sum=0.0d0; d2_sum=0.0d0
              do i=1,nij
                i1=valid(i,1)
                j1=valid(i,2)
                !write(*,*)"0",i1,j1,mat2%deltR(i1,j1),mat2%deltR(j1,i1) 
                !mat2%deltR(i1,j1)=  mat2%deltR(i1,j1)+ XP(i+1)
                !mat2%deltR(j1,i1)= -mat2%deltR(i1,j1)
! recover the anti-symmetric rotations
                mat2%deltR(i1,j1)=  XP(i+1)/1.0d0
                mat2%deltR(j1,i1)= -XP(i+1)/1.0d0
!                print *,"deltaR = XP related to XP " 
                write(1,*)"deltaR in",iloop,"micro-iter",&
                          i1,j1,mat2%deltR(i1,j1)
                if(dabs(mat2%deltR(i1,j1)).gt.0.10)then
                  write(2,*)"notice : large rotation between",i1,j1,&
                            " : ",mat2%deltR(i1,j1) 
!                  if(iloop.gt.10)then
!                    write(2,*)"reset to 0.01/-0.01 as maxinum" 
!                    if(mat2%deltR(i1,j1).gt.0.10)then
!                      mat2%deltR(i1,j1)= 0.010d0
!                    else
!                      mat2%deltR(i1,j1)=-0.010d0
!                    end if 
!                  end if 
                end if  
                d0=d0+dabs(XP(i+1))
                d2=d2+(XP(i+1))**2
              end do
              d0_sum = d0_sum + d0
              d2_sum = d2_sum + d2

! recover the degenrate rotations
              is0=0
              do isym=1,orb%nsub
                js0=0
                do jsym=1,orb%nsub
                  if(redu%irreps(isym).eq.redu%irreps(jsym).and.isym.lt.jsym)then
                    do i=1,orb%total(isym)
                      do j=1,orb%total(isym)
                        ii=i+is0; jj=j+is0
                        kk=i+js0; ll=j+js0
                ! avoid possible over-shooting by 1/2 (?)
                mat2%deltR(ii,jj)=mat2%deltR(ii,jj)/2.0d0  
                mat2%deltR(kk,ll)=mat2%deltR(ii,jj)*redu%sign_mat(kk,ll)
               write(1,*)ii,jj,mat2%deltR(ii,jj),kk,ll,mat2%deltR(kk,ll)
                      end do
                    end do 
                  end if
                  js0=js0+orb%total(jsym)
                end do
                is0=is0+orb%total(isym)
              end do 
!              stop


! Calculate the residual for orbital rotation

               allocate( wB(nij));   wB=0.0d0
               allocate(wBT(nij));  wBT=0.0d0
              allocate(ARRA(nij)); ARRA=0.0d0
               allocate(SFT(nij));  SFT=0.0d0
              call dB_update(ndim1,nij,orb%nsub,orb%total,orb%occ,&
                      valid(1:nij,:),mat2%B,mat2%G,mat2%deltR,wB,wBT)
!              print *,"~B"
!              call print_mat(1,nij, wB,6)
!              print *,"~B^T"
!              call print_mat(1,nij,wBT,6)
!              print *," mat2%A "
!              call print_mat(1,nij,ARRA,6)

              call ARplusRA(norb,nij,orb%nsub,orb%total,orb%occ,&
                         valid(1:nij,:),mat2%A,mat2%deltR,ARRA)
              write(2,*)"      g/lamda  "
              call print_mat(1,nij, GP/D_lamda,2)
              !call print_mat(1,nij, GP,2)

!              print *," Hessian "
!              call print_mat(nij,nij,HP,6)

              do i=1,nij
                SFT(i)=D_epsilon*XP(i+1)
              end do  
!              print *," Shift "
!              call print_mat(1,nij,SFT,6)

              write(2,*)" YP without g/lamda"
              YP=wB-wBT-0.5d0*ARRA-SFT
              call print_mat(1,nij,YP,2)

              deallocate(WB,WBT,ARRA,SFT)

              deallocate(HPA)
              deallocate(XP)
              deallocate(D)

            end do

! ==========================================================================
            call MdXM(norb,mat2%U,MAT2%B,TM1)
            TM2=TM1+mat2%A
            call MdXM(norb,mat2%deltR,TM2,TM1) 
            dv=0.0d0
            call trace(norb,TM1,dv)
            ETdR=ETdR+0.5d0*dv

            call MXM(norb,mat2%deltR,mat2%deltR,TM1)
            call MXM(norb,TM1,mat2%A,TM2)
            dv=0.0d0
            call trace(norb,TM2,dv)
            ETdR=ETdR+0.5d0*dv
            write(2,*)&
            "The predicted E2(T,dR) for the next micro-iter : ",ETdR
! =========================================================================
            Econv=ETdR-FRONRE0
!            call print_mat(norb,norb,mat2%deltR,6)

            allocate(T1(norb,norb)); T1=0.0d0
            allocate(T2(norb,norb)); T2=0.0d0
            call CAL_T(norb,mat2%deltR,1.0d0,150,T1)
            call T_TO_U(norb,T1,T2) 
            call MXM(norb,mat2%U,T2,T1)
            mat2%U=T1
            mat2%T=mat2%U
            do i=1,norb
              do j=1,norb
                if(dabs(mat2%U(i,j)).lt.1.0d-16)then
                  mat2%U(i,j)=0.0d0
                end if
              end do
              mat2%T(i,i)=mat2%U(i,i)-1.0d0
            end do  
            deallocate(T1,T2) 
  
!            print *,"============ The U matrix : ============"
!            call print_mat(norb,norb,mat2%U,6)
 

!            goto 3333

!3333        continue 

            d0_sum0=0.0d0
            d2_sum0=0.0d0
            do i=1,norb
              do j=1,norb
                d0_sum0=d0_sum0 + dabs(mat2%T(i,j))
                d2_sum0=d2_sum0 + mat2%T(i,j)**2
              end do
            end do 

            write(1,*)"===== cycle  ",icycle,iloop
            write(1,*)"== Current rotation parameters =="            
            write(1,*)"Sum(|deltR|)   ",d0_sum
            write(1,*)"Sum(deltR**2)  ",d2_sum
            write(1,*)"Sum(|T|)       ",d0_sum0
            write(1,*)"Sum(T**2)      ",d2_sum0
            write(1,*)"================================="            

            if(LMCSCF)then
              if(allocated(valid))then
                deallocate(valid)
              end if
            end if

! Check for convergence

          else            ! Here use the directly inversion

            YP=-1.0d0*YP

! And the augmented x-vector
            allocate(XP(nij+1))
            allocate(D(nij))
            XP=0.0d0;XP(1)=1.0d0
             D=0.0d0
            allocate(T1(nij,nij)); T1=0.0d0

            T1=HP
!           call print_mat(nij,nij,HP,6)
            call eigtql2(nij,T1,D)
!            write(*,*)"========= D ============" 
!            call print_mat(1,nij,D,6)

            do i=1,nij
              do k=1,nij
                do l=1,nij
                  if(D(l).gt.1.0d-9)then
                    dv=T1(i,l)*T1(k,l)*YP(k)/D(l)
!                    if(dabs(dv).gt.0.1)then
!                      XP(i+1)=XP(i+1)+ sign(0.1d0,dv) 
!                    else
                      XP(i+1)=XP(i+1)+dv
!                    end if
                  else if(D(l).lt.-1.0d-9)then
                    dv=T1(i,l)*T1(k,l)*YP(k)/D(l)
!                    if(dabs(dv).gt.0.1)then
!                      XP(i+1)=XP(i+1)- sign(0.1d0,dv) 
!                    else
                      XP(i+1)=XP(i+1)-dv
!                    end if
                  end if
                end do
              end do
            end do
            deallocate(T1) 

!            write(6,*)"XP - 0  ",XP
            ! Set XP(1) as 1 as the flag in Augmented Hessian algorithm
!            call MODIFIED_GRAM_SCHMIDT(nij+1,1,XP)
!            write(6,*)iloop,"iloop,XP - 1  ",XP
!            stop

            if(iloop.ge.0)then
              do
                call MODIFIED_GRAM_SCHMIDT(nij+1,1,XP)
                if(dabs(XP(1)).lt.damping)then
                  XP(1)=XP(1)*1.25d0
                  XP(2:nij+1)=XP(2:nij+1)/1.25d0
                else
                 !damping=dabs(XP(1)) ! reset damping
                  exit
                end if
              end do
            end if
!            write(6,*)iloop,"iloop,XP - 1  ",XP

            do i=1,nij
              X(porb(i))=XP(i+1)
            end do

            ij=0;
            d0=0.0d0; d2=0.0d0;
            d0_sum=0.0d0; d2_sum=0.0d0
            do i=1,ndim1
              do j=1,ndim1
                ij=ij+1
                mat2%deltR(i,j)=X(ij)
                d0=d0+dabs(X(ij))
                d2=d2+(X(ij))**2
              end do
            end do
            d0_sum = d0_sum + d0/2.0d0
            d2_sum = d2_sum + d2/2.0d0

!           call B_update(ndim1,orb%nsub,orb%total,orb%occ,mat2%B,mat2%G,mat2%deltR,waveB)
!           call residual_WMK(ndim1,0.0d0,mat2%U,waveB,mat2%A,mat2%deltR,ResidM)

!           call print_mat(norb,norb,ResidM,6)
!           stop
            ! save previous U,T matrix          
            !Ubk=mat2%U 
            !Tbk=mat2%T 

            allocate(T1(norb,norb)); T1=0.0d0
            allocate(T2(norb,norb)); T2=0.0d0
            call CAL_T(norb,mat2%deltR,d0,50,T1)
            call T_TO_U(norb,T1,T2)
            call MXM(norb,mat2%U,T2,T1)

            mat2%U=T1
            mat2%T=mat2%U
            do i=1,norb
              do j=1,norb
                if(dabs(mat2%U(i,j)).lt.1.0d-16)then
                  mat2%U(i,j)=0.0d0
                end if
              end do
              mat2%T(i,i)=mat2%U(i,i)-1.0d0
            end do
            deallocate(T1,T2)

            deallocate(D)
            deallocate(XP)

          end if

          write(1,*)"GP"
          call print_mat(nij,1,GP,1)
!          write(1,*)"YP"
!          call print_mat(nij,1,YP,1)
!          write(6,*)"before deallocate "
!          call flush(6)

          deallocate(porb)
          deallocate(HP)
          deallocate(GP,YP,YQ)
!          deallocate(valid)

!          write(6,*)"after deallocate "
!          call flush(6)

          d0_sum0=0.0d0
          d2_sum0=0.0d0
          do i=1,norb
            do j=1,norb
              d0_sum0=d0_sum0 + dabs(mat2%T(i,j))
              d2_sum0=d2_sum0 + mat2%T(i,j)**2
            end do
          end do

!          write(6,*)"before keep data "
!          call flush(6)

          write(1,*)"===== cycle  ",icycle,iloop
          write(1,*)"== Current rotation parameters =="
          write(1,*)"Sum(|deltR|)   ",d0_sum
          write(1,*)"Sum(deltR**2)  ",d2_sum
          write(1,*)"Sum(|T|)       ",d0_sum0
          write(1,*)"Sum(T**2)      ",d2_sum0
          write(1,*)"================================="
       
          write(2,*)"Energy converged(SCF)   : ", Econv+FRONRE0
          write(2,*)"Energy micro-DMRG/CI    : ", E0_micro+FRONRE0
          write(2,*)"Energy change : ", Econv-E0_micro
          if(CP_integrals.and.ith_INTE.ne.0)then 
            write(2,*)"Energy change between two MPS/CI: ",&
             EconvCP1-EconvCP0
          end if
          walltime(24) = wtime()
          microlabel = ""; write(microlabel,'(a,i4)')"microiteration #",iloop
          call timing_reporter(3,trim(microlabel),walltime(24)-walltime(18))
          call timing_reporter(30,trim(microlabel),walltime(24)-walltime(18))

          if(iloop.eq.iloopmax)then
            write(2,*)"reach the maximum loop"
            exit 
          end if
  
!          if(iloop.gt.10)then
            if(CP_integrals.and.ith_INTE.eq.0)then            
              if(dabs(Econv-E0_micro).lt.thrs%e)then
                write(2,*)&
            "  : reach the desired energy accuracy in micro-loops"
                exit
              end if
            end if
            if(CP_integrals.and.ith_INTE.ne.0)then            
              if(dabs(EconvCP0-EconvCP1).lt.thrs%e)then
                write(2,*)&
            "CP: reach the desired energy accuracy in micro-loops"
                exit
              end if
            end if
!          end if
!          if(dabs(D_epsilon).lt.0)then
!            write(*,*)"tiny rotation (0)"
!            exit 
!          end if

        end do

!        print *," ====== if never coupled case =======" 
        if(LNOcoupled)then
          UMAT_FINAL=mat2%U
        end if

        d0_sum0=0.0d0
        d2_sum0=0.0d0
        do i=1,norb
          do j=1,norb
            if(i.ne.j)then
              d0_sum0=d2_sum0+dabs(UMAT_FINAL(i,j))
              d2_sum0=d2_sum0+UMAT_FINAL(i,j)**2
            end if
          end do       
        end do 

        write(*,2434)Econv+fronre0, E0_micro+fronre0, dsqrt(d2_sum0)
2434    format(4X,f14.8,4X,f14.8,4X,f8.4)

        write(1,*)" Finished the ",icycle,"all micro-iteration"
        write(2,*)" Finished the ",icycle,"all micro-iteration"
!        stop

        Rabs(icycle)=d0_sum0
        write(2,*)icycle,d0_sum,d2_sum,d0_sum0,d2_sum0 
!        call print_mat(ndim1,ndim1,mat2%T,6)
!        stop

! Recover the A,B,G matrix if changed 
        mat2%A=Aini
        mat2%B=Bini
        mat2%G=Gini

9999    continue

        deallocate(Gini)
        if(allocated(valid))then
          deallocate(valid)
        end if
        deallocate(X,Y,G,H)

        walltime(17) = wtime()
        microlabel = ""
        write(microlabel,'(a,i4,a)')"microiteration"//&
        " section to complete (total # of micro-iterations:",iloop,")"
        call timing_reporter(3,trim(microlabel),walltime(17)-walltime(25))
        call timing_reporter(30,trim(microlabel),walltime(17)-walltime(25))

      End Subroutine Solver_WMK2


