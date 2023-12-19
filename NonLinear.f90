      Subroutine NonLinear_solver(icycle)

! Solving the Nonlinear equatuon

        use global_control
        use matrix        

        integer::icycle

        if(method(icycle).eq."microit")then 
          call derivaties(0)
          call iteration(icycle,.false.) 
        else if(method(icycle).eq."direct2")then
          call derivaties2(0,0)
          call iteration(icycle,.false.)  
        else if(method(icycle).eq."direct3")then
          call derivaties3(0,0)
          call iteration(icycle,.false.)          
        else if(method(icycle).eq."orbhess")then
          call solver_hessian(norb,icycle)
        else if(method(icycle).eq."occhess")then
          call solver_hessian_occupy(nact,icycle)
          call derivaties(0)
          call iteration(icycle,.true.)          
        else if(method(icycle).eq."cp-hess")then
! to Stefan : Coupled Newton-Raphson algorithm
!             Only used for checking Hessian (e.g. positive define)
!             (commented at P15 of the draft)
          write(1,*)" -->  Start the coupled (hessian) terms" 
          call coupling(icycle)
!          if(ith_hess.eq.2)then
            ! Use the origonal WM solver
            !call derivaties(0)
            !call iteration(icycle,.false.) 
            T=T_origin 
            U=U_origin 
            mat2%U=Umat_final
!          end if 
          write(1,*)" --> Finish the coupled (hessian) terms" 
!          stop
        else if(method(icycle).eq."ItdeltR")then
          call Solver_Augmented_Hessian_deltaR(icycle)
        else if(method(icycle).eq."WMKUBAR")then
!          write(6,*)"Before Solver WMK "; call flush(6)
! to Stefan : Werner-Meyer-Knowles algorithm
          call Solver_WMK(icycle)
!          write(6,*)" After Solver WMK "; call flush(6)
          T=T_origin 
          U=U_origin 
          mat2%U=Umat_final
!          write(6,*)"The Umat_final"
!          call print_mat(norb,norb,Umat_final,6)
        else if(method(icycle).eq."augment")then
!          write(6,*)"Before Solver AH "; call flush(6)
!          write(6,*)"This is a simple/debugging AH, no redundant limit"
          call Solver_Augmented_Hessian(icycle) 
!          write(6,*)" After Solver AH "; call flush(6)
        else
          call derivaties4(0,0)
          call iteration(icycle,.false.)  
        end if


      end subroutine Nonlinear_solver


      Subroutine DIIS_acceleration(icycle,m)
! Used for the macro-iters-diis

        use global_control
        use matrix

        integer::icycle,M
  
        double precision TM1(norb,norb)
        double precision TM2(norb,norb)
        double precision TM3(norb,norb)

        double precision,allocatable::dP(:,:)
        double precision,allocatable::BM(:,:),eigBM(:)
        double precision,allocatable::b(:),c(:)

        TM1=0.0d0
        TM2=0.0d0
        TM3=0.0d0 

        allocate(Udiis(icycle)%Ua(norb,norb))  
        Udiis(icycle)%Ua=0.0d0

        if(icycle.le.1)then
          Udiis(icycle)%Ua=Umat_final
        else
          TM1=Udiis(icycle-1)%Ua
          call MXM(norb,TM1,mat2%U,Udiis(icycle)%Ua)
        end if

        call print_mat(norb,norb,Udiis(icycle)%Ua,6)  

        if(icycle.lt.m+2)then
          LDIIS_UPDATED=.false.
        else
          LDIIS_UPDATED=.true.

          write(6,*)"Entering the diis",m 
          allocate(dP(norb**2,m))
          dP=0.0d0
          ! get previous M error vectors
          do i=1,m 
            kl=0            
            do k=1,norb
              do l=1,norb
                kl=kl+1
                dp(kl,i) = Udiis(icycle-m+i)%Ua(k,l) &
                        -Udiis(icycle-m+i-1)%Ua(k,l)
              end do
            end do            
          end do  
          write(6,*)"After the error vectors" 
          ! Get the C_i vectors
          allocate(BM(m+1,m+1));   BM=0.0d0
          allocate(eigBM(m+1)); eigBM=0.0d0
          allocate(b(m+1));         b=0.0d0
          allocate(c(m));           c=0.0d0
          do i=1,m
            do j=1,m
              do kl=1,norb**2              
                BM(i,j)=BM(i,j)+dp(kl,i)*dp(kl,j) 
              end do
            end do
          end do
          write(6,*)"After the BM matrix" 
          BM(m+1,1:m)=-1.0d0 
          BM(1:m,m+1)=-1.0d0

          ! if never ill-condition, the inversion can be used
          if(.true.)then
            !call print_mat(m+1,m+1,BM,6)
            call inv(m+1,BM) 
            !write(6,*)"============"
            !call print_mat(m+1,m+1,BM,6)
            c=-1.0d0*BM(1:m,m+1)
            write(2,*)"===   DIIS  ====="
            write(2,*)c
          else
          ! Should be the same way as the maxtrix function

          !call print_mat(m+1,m+1,BM,6)
!          call eigtql2(m+1,BM,eigBM)
!          call print_mat(m+1,m+1,BM,6)
!          call print_mat(1,m+1,eigBM,6)

          end if
          write(*,*)"DIIS not fully implemented yet" 
          stop 

          ! Overwrite the icycle-th U matrix          
          do i=1,m
            Udiis(icycle)%Ua=    Udiis(icycle)%Ua + &
                            c(i)*Udiis(icycle-m+i-1)%Ua
          end do
          call GRAM_SCHMIDT(norb,Udiis(icycle)%Ua)
          !stop

          deallocate(eigBM)          
          deallocate(BM)          
          deallocate(b)          
          deallocate(c)          
          deallocate(dp)                 
        end if

   
      End Subroutine DIIS_acceleration



