! Special Davidson only used in Werner-Meyer-Knowles algorithm
      Subroutine Davidson_WMK(ndim,Hdav,v,eig,lamda,&
                              valid,Gini,Bini,nMAX,residual)

        use global_control
        use matrix

        integer::ndim,nMAX
        double precision::Hdav(ndim,ndim)
        double precision::v(ndim)
        double precision::eig(ndim)
        double precision::residual

! Used in Augmented Hessian solver
        double precision::lamda  ! If AH 
        double precision epsilon ! level-shift parameter
        double precision g(ndim-1)

! Used in Werner-Meyer-Knowles
!        integer::porb(ndim)
        integer::valid(ndim-1,2)
        double precision Y(ndim),YM(ndim,nMAX)
        double precision SM(norb,norb)
        double precision deltaR(norb,norb)
!        double precision AdR(norb,norb),dRA(norb,norb)
!        double precision UBBU(norb,norb),AA(norb,norb)        
        double precision Gini(norb,norb,norb,norb)
        double precision Bini(norb,norb)
        double precision ResidM(norb,norb) 

        logical ciupdate,Lexit 

! Used in normal davidson
        double precision dv,dt
        double precision Hb(ndim,nMAX),b(ndim,nMAX)
        double precision Q(ndim,nMAX),X(ndim),D(ndim)
        double precision,allocatable::A(:,:),eigA(:)

! used in residual calculation
        double precision,allocatable::wB(:),wBT(:),ARRA(:),SFT(:)


          Hb=0.0d0
           b=0.0d0
           Q=0.0d0
           X=0.0d0
         eig=0.0d0

          write(6,*)"entering the Davidson_WMK"
          write(6,*)valid(:,1)
          write(6,*)valid(:,2)
          call print_mat(ndim,ndim,Hdav,6) 
          call flush(6)

        ! Get gradient 
          g=Hdav(2:ndim,1)

        ! Initialize eigen vector
          v   =0.0d0
          v(1)=1.0d0/lamda
          !call random_seed()
          do i=2,ndim
            !call random_number(v(i))
            !v(i)=v(i)/50.0d0
            v(i)=-g(i-1)/Hdav(i,i)
          end do

! Put the converged vector to check
!          v(1)= 0.79349327351601939
!          v(2)= 0.11774709E-01
!          v(3)= -.33680119E-01 
!          v(4)= 0.60334024     
!          v(5)= -.71245954E-01

          b(:,1)=v


        ! Modify the Hessian
        !  Hdav(2:ndim,2:ndim)=Hdav(2:ndim,2:ndim)/lamda !all

        ! Start the davidson
          do idvs=1,nMAX

            if(idvs.gt.1)then  ! New vector
              b(:,idvs)=X
            end if
            call MODIFIED_GRAM_SCHMIDT(ndim,idvs,b(:,1:idvs))

            ! only for the idvs.eq.1 case
            if(idvs.eq.1)v=b(:,1)
          ! Get level-shift parameter  ???
            ! idvs-th orthonormalized vector
             call inner_product(ndim-1,g,b(2:ndim,idvs),dv)
            ! idvs-th corrected vector
            !call inner_product(ndim-1,g,v(2:ndim),dv)

            lamda=1.0d0/b(1,idvs) 
            write(6,*)"curretn  lamda ",   lamda
            epsilon=lamda*dv   
            write(6,*)"curretn epsilon", epsilon

!            if(idvs.eq.1)then
              !update mat2%deltR 
              !mat2%deltR=0.0d0

! Calculate the residual for orbital rotation
              nij=ndim-1 

              ! idvs-th orthonormalized vector
              call vec2deltaR(norb,nij,Valid,b(2:ndim,idvs),mat2%deltR)
              ! idvs-th       corrected vector
              !call vec2deltaR(norb,nij,Valid,v(2:ndim),mat2%deltR)
              write(6,*)"current deltaR" 
              call print_mat(norb,norb,mat2%deltR,6) 
              !stop 

               allocate( wB(nij));   wB=0.0d0
               allocate(wBT(nij));  wBT=0.0d0
              allocate(ARRA(nij)); ARRA=0.0d0
               allocate(SFT(nij));  SFT=0.0d0
              call dB_update(norb,nij,orb%nsub,orb%total,orb%occ,&
                      valid(1:nij,:),mat2%B,mat2%G,mat2%deltR,wB,wBT)
!              print *,"~B"
!              call print_mat(1,nij, wB,6)
!              print *,"~B^T"
!              call print_mat(1,nij,wBT,6)

              call ARplusRA(norb,nij,orb%nsub,orb%total,orb%occ,&
                         valid(1:nij,:),mat2%A,mat2%deltR,ARRA)

!              print *," ARRA "
!              call print_mat(1,nij,ARRA,6)

              print *,"g/lamda"
              call print_mat(1,nij, g/lamda,6)

!              print *," Hessian "
!              call print_mat(nij,nij,HP,6)

              do i=1,nij
                ! idvs-th orthonormalized vector
                SFT(i)=epsilon*b(i+1,idvs)
                ! idvs-th       corrected vector
                !SFT(i)=epsilon*v(i+1)
              end do 
!              print *," Shift "
!              call print_mat(1,nij,SFT,6)

              ! calculate the residual
              YM(:,idvs)=0 
              YM(2:ndim,idvs)=wB-wBT-0.5d0*ARRA-SFT
              print *," YM without g/lamda "
              call print_mat(1,nij,YM(2:ndim,idvs),6)
              YM(2:ndim,idvs)=YM(2:ndim,idvs)+g/lamda
              print *," ======= YM in davidson loop ======== ",idvs
              write(6,*)"current residual ",YM(:,idvs)
!              stop
!              call print_mat(1,nij,YM(2:ndim,idvs),6)
!              write(6,*)"11 Av-residual YP in micro-AH",iAH,"cycle",dv

              deallocate(WB,WBT,ARRA,SFT)

!              stop

              !ciupdate=.false.
              !idavidson=0 
              !call Operator_update(idavidson,ciupdate,Lexit)
!              stop

!              if(idvs.eq.1)then
!              else
                Hdav(2:ndim,1)=YM(2:ndim,idvs)
                Hdav(1,2:ndim)=YM(2:ndim,idvs)
!              end if

!            end if

            !call print_mat(ndim,ndim,Hdav,6) 
            !call flush(6)
            !stop
         
            call MXV(ndim,ndim,Hdav,b(:,idvs),Hb(:,idvs))

            !write(6,*)" b",  b(:,idvs)
            !write(6,*)"Hb", Hb(:,idvs)
            !call flush(6)
!           if(idvs.eq.2)stop

            allocate(eigA(idvs));  eigA=0.0d0
            allocate(A(idvs,idvs));   A=0.0d0
            do i=1,idvs
              do j=1,idvs
                call inner_product(ndim,b(:,i),Hb(:,j),A(i,j))
              end do
            end do
!            call print_mat(idvs,idvs,A,6)

            call eigtql2(idvs,A,eigA)
!            write(*,*)"Vec-A(1) and eigA(1)",eigA            
!            call print_mat(idvs,idvs,A,6)

! Calculated eigenvalue and eigenvector (keep positive eigE and eigV ??)
!            eig(1)=eigA(ipos)
            ipos=1 
 
            v=0.0d0
            Y=0.0d0
            do i=1,idvs
              v=v+A(i,ipos)*b(:,i)      ! updated DeltaR
              Y=Y+A(i,ipos)*YM(:,idvs)  !       Residual
              !Y=YM(:,idvs)  !       Residual
            end do
            !write(6,*)"Current vector  ", v
            write(6,*)"Current residual", Y

            call LENGTH(Y,ndim,dv)
            write(6,*)"length of residual  :  ", dv
            !call flush(6)

            X=0.0d0  ! Get the new expansion vector
            do i=2,ndim
              X(i)=-1.0d0*Y(i)/(Hdav(i,i)-epsilon)
            end do

            if(dv.lt.residual.or.idvs.eq.nMAX)then
              exit
              !convergenced
            end if

            deallocate(eigA)
            deallocate(A)

          end do

      End Subroutine Davidson_WMK

! ==============================================================

      Subroutine Davidson(ndim,H,v,eig,lamda,AH,nMAX,residual)

        integer::ndim,nMAX
        double precision::H(ndim,ndim)
        double precision::v(ndim)  
        double precision::eig(ndim)  
        double precision::residual

! Use in Augmented Hessian solver
        double precision::lamda  ! If AH 
        double precision epsilon ! level-shift parameter
        double precision g(ndim-1) 
        logical::AH ! If augmented Hessian algorithm

! Use in normal davidson
        double precision dv,dt
        double precision Hb(ndim,nMAX),b(ndim,nMAX)
        double precision Q(ndim,nMAX),X(ndim),D(ndim)
        double precision,allocatable::A(:,:),eigA(:)

          Hb=0.0d0 
           b=0.0d0
           Q=0.0d0
           X=0.0d0
         eig=0.0d0

!        write(6,*)
!        write(6,*) "Start the Davidson part with residual",residual
!        call flush(6) 

        if(AH)then

          ! Initialize eigen vector
          v   =0.0d0
          v(1)=1.0d0/lamda 
          call random_seed()
          do i=2,ndim 
            call random_number(v(i))
            v(i)=v(i)/50.0d0 
          end do
          b(:,1)=v

          
          g=H(2:ndim,1)

          !  write(*,*)"  lamda is ", lamda
          ! Modify the Hessian to control step-length
          !  H(2:ndim,2:ndim)=H(2:ndim,2:ndim)!/lamda  

          write(2,*)" gradients vector :"
          
          call print_mat(1,ndim-1,g,2)


          
          do idvs=1,nMAX

            ! Get level-shift parameter
            call inner_product(ndim-1,g,v(2:ndim),dv)
            epsilon=lamda*dv/v(1)  ! 
            !            epsilon=0.0 
            !            write(2,*)"epsilon is ", epsilon
            !            do i=1,ndim    ! 1/2  First  half
            !              H(i,i)=H(i,i) 
            !            end do
        

            if(idvs.gt.1)then
              b(:,idvs)=X
            end if

            call MODIFIED_GRAM_SCHMIDT(ndim,idvs,b(:,1:idvs))

            !            lamda=1.0d0/v(1)  ! Don't know how to reduce the Hessian ?
            !            H(2:ndim,2:ndim)=H(2:ndim,2:ndim)/lamda !all

            call MXV(ndim,ndim,H,b(:,idvs),Hb(:,idvs))

            !            write(6,*)" b",  b(:,idvs)
            !            write(6,*)"Hb", Hb(:,idvs)
            !            call flush(6)

            allocate(eigA(idvs));  eigA=0.0d0
            allocate(A(idvs,idvs));   A=0.0d0
            do i=1,idvs
              do j=1,idvs
                call inner_product(ndim,b(:,i),Hb(:,j),A(i,j))
              end do
            end do

            

            call eigtql2(idvs,A,eigA)

            write(2,*)"eig-A"
            call print_mat(1,1,eigA,2)

            ipos=1  ! choose the positive eigvalue and eigvector
            if(dabs(eigA(1)).gt.0.1d0)then
              !              do i=1,idvs
              !                if(dabs(eigA(i)).lt.0.1d0)then
              !                  ipos=i
              !                  write(1,*)"Invalid eigvalue detected, shift to the ",&
              !                  ipos,"eigvector in Davidson(AH case)"
              !                  exit
              !                end if              
              !              end do
            end if

            !            write(6,*)"eigA",eigA(1)
            !            write(1,*)"Vec-A"
            !            call print_mat(idvs,min(idvs,5),A,1)

            ! Q is the residual
           
            Q(:,idvs)=0.0d0
            do i=1,idvs
             Q(:,idvs)=Q(:,idvs)+A(i,ipos)*Hb(:,i)-A(i,ipos)*eigA(ipos)*b(:,i)
            end do
            

            X=0.0d0
            do i=1,ndim
              !X(i)=1.0d0/(eigA(ipos)-H(i,i))*Q(i,idvs) ! Standard Davidson
              if(i.ne.1)then
                X(i)=-1.0d0*Q(i,idvs)/H(i,i)  ! More efficient for AH 
              end if
            end do

            call LENGTH(Q(:,idvs),ndim,dv)
            !            write(6,*)"Finish",idvs,"-th iteration with |r|",dv
            !            call flush(6)          
 
            ! Calculated eigenvalue and eigenvector 
            eig(ipos)=eigA(ipos)
            v=0.0d0
            do i=1,idvs
              v=v+A(i,ipos)*b(:,i)
            end do
            !            write(*,*)"v  :  ",v
            !            write(*,*)"   "
                      !if(dv.lt.residual)then
            if(dv.lt.residual.or.idvs.eq.nMAX)then
              dv=epsilon-eig(ipos)
                write(2,*)"eigenvalue(diag) and eigenvalue(shift) "   
                write(2,*) eigA(ipos),epsilon
                write(2,*)"eigenvalue "   
                call print_mat(1,idvs,eigA,2)
                !                write(2,*)"eigenvec   "   
                !                call print_mat(idvs,idvs,A,2)
              if(dabs(dv).lt.1.0d-5)then !! ???? 
                write(2,*)"resonable solution found"
              else
                write(2,*)"Warning : wrong in Davidson ?? "
                !                write(2,*)"          cycle this AH step"
                !                V=0.0d0
                !                v(1)=100.0d0
              end if
              exit
              !convergenced
            end if

            deallocate(eigA)
            deallocate(A)

            

          end do

        else

          
          v   =0.0d0
          call random_seed()
          do i=1,ndim
            call random_number(v(i))
            v(i)=v(i)*0.01 
          end do
          V(1)=1.0d0
          b(:,1)=v

          
          do idvs=1,nMAX

            if(idvs.gt.1)then
              b(:,idvs)=X
            end if

            call MODIFIED_GRAM_SCHMIDT(ndim,idvs,b(:,1:idvs))
            call MXV(ndim,ndim,H,b(:,idvs),Hb(:,idvs))

            !            write(6,*)" b",  b(:,idvs)
            !            write(6,*)"Hb", Hb(:,idvs)
            !            call flush(6)

            allocate(eigA(idvs));  eigA=0.0d0
            allocate(A(idvs,idvs));   A=0.0d0
            do i=1,idvs
              do j=1,idvs
                call inner_product(ndim,b(:,i),Hb(:,j),A(i,j))
              end do
            end do
            call eigtql2(idvs,A,eigA)

            !            write(6,*)"eigA",eigA(1)
            !            call print_mat(idvs,idvs,A,6)

            ! Q is the residual
            do i=1,idvs              
             Q(:,idvs)=Q(:,idvs)+A(i,1)*Hb(:,i)-A(i,1)*eigA(1)*b(:,i)
            end do
          !            write(*,*)"Q  :  ",Q(:,idvs)

            X=0.0d0
            do i=1,ndim
              X(i)=1.0d0/(eigA(1)-H(i,i))*Q(i,idvs)
            end do

            call LENGTH(Q(:,idvs),ndim,dv)
          !            write(6,*)"Finish",idvs,"-th iteration with |r|",dv
          !            call flush(6)           
 
          ! Calculated eigenvalue and eigenvector 
            eig(1)=eigA(1)
            v=0.0d0
            do i=1,idvs 
              v=v+A(i,1)*b(:,i)
            end do
            if(dv.lt.residual.or.idvs.eq.nMAX)then
              exit
              !convergenced
            end if

            deallocate(eigA)
            deallocate(A)

          end do

        end if

      end Subroutine Davidson

! 

      Subroutine vec2deltaR(ndimA,ndimE,Valid,v,deltaR)

        integer::ndimA,ndimE
        integer::Valid(ndimE,2)
        double precision::v(ndimE)
        double precision::deltaR(ndimA,ndimA)
 
        deltaR=0.0d0

        do i=1,ndimE
          i1=valid(i,1)
          j1=valid(i,2)
          deltaR(i1,j1)=   v(i)
          deltaR(j1,i1)= -deltaR(i1,j1)
        end do

      End Subroutine vec2deltaR


