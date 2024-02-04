      Subroutine Solver_Augmented_Hessian_deltaR(icycle)

        use global_control
        use matrix

        integer::icycle

        integer,allocatable::porb(:)
        double precision d0,d0_sum,d0_sum0,d2,d2_sum0,d2_sum
        double precision dv,redundant,XP1last,phase,phaselast,THdeltaR
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:)
        double precision,allocatable::H(:,:),HP(:,:),HPA(:,:)

        double precision,allocatable::T1(:,:),T2(:,:)

        ! Used for the augmented Hessian part
        double precision lamda        ! damping parameter that control the level shift
        double precision levelshift   ! in the loop

        ! Save the initial A and B
        double precision Aini(norb,norb)
        double precision Bini(norb,norb)
        double precision,allocatable::Gini(:,:,:,:)
        !Save the valid rotation index
        integer,allocatable::valid(:,:)
        logical Lvalid
        logical Lfinish

        double precision,allocatable::G(:)


        ! Parameter initialization
        !lamda=1.55d0 ! "1" means the standard AH    : For H2
        !levelshift=0.99999999d0  ! restrict initial : For H2
        lamda=1.0d0 ! "1" means the standard AH 
        ! The larger lamda, the shorter step
        levelshift=0.999999  ! restrict initial

        Lfinish=.false.
          ! initialize of part of Hessian matrix that will be used in this solver

        ndim1=norb

        allocate(X(ndim1**2)); X=0.0d0
        allocate(Y(ndim1**2)); Y=0.0d0
        allocate(G(ndim1**2)); G=0.0d0
        allocate(H(ndim1**2,ndim1**2)); H=0.0d0

        ! Get the initial A, B matrix 
        !        call derivaties(0)
        Aini=mat2%A
        !        Bini=mat2%B
        ! The initial G matrix
        allocate(Gini(norb,norb,norb,norb)); Gini=0.0d0
        Gini=mat2%G

        ! This is the orbiatl-rotation gradients 
        redundant=1.0d-9

          grad_digit=0.0d0
          ij=0
          nij=0
          do i=1,ndim1
            do j=i+1,ndim1 ! only half is needed
                !do j=i,ndim1 ! only half is needed
               ij=ij+1
               
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

          nij0=nij

! re-run it to get the valid rotations index
          allocate(valid(nij,2)); valid=0
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

          write(1,*)"valid rotations"
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,1)
          end do
          write(1,*)" "
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,2)
          end do
          call flush(1)

          nij=nij0/2

!        ij=0
!        nij=0
!        do i=1,ndim1
!          do j=1,ndim1
!             ij=ij+1
!             Y(ij)=1.0d0*(mat2%A(i,j)-mat2%A(j,i))
!             if(dabs(Y(ij)).lt.redundant)then
!             else
!               nij=nij+1                
!             end if
!          end do
!        end do
!        nij0=nij
! re-run it to get the valid rotations index
!        allocate(valid(nij,2)); valid=0
!        ij=0
!        kl=0
!        do i=1,ndim1
!          do j=1,ndim1
!             ij=ij+1
!             if(dabs(Y(ij)).lt.redundant)then
!             else
!               kl=kl+1
!               valid(kl,1)=i
!               valid(kl,2)=j
!             end if
!          end do
!        end do
!
!        nij=nij0/2

        write(*,*)
        d0_sum=0.0d0 
        d2_sum=0.0d0 
        d0_sum0=0.0d0 
        d2_sum0=0.0d0 
        iloop=0
        nsmooth=0
        XP1last=0.0d0
        phaselast=1.0d0
        do 
          iloop=iloop+1
          !write(*,*)"in loop",iloop
          ! Update the B matrix
          call derivatives_B_update(Aini)
          Bini=mat2%B
          ! Get the ~A , same as A if initial U.eq.1
          call Derivaties_deltaR(icycle,iloop,Gini,Bini)
          !write(*,*)"After derivative deltaR" 
          !call flush(6)
          ! Get the initial Hessian for deltaR, same as R if initial
          call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)

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
  
          !write(*,*)"Before get the delta-orbiatl-rot" 
          !call flush(6) 

          allocate(porb(nij));   porb=0
          allocate(YP(nij));     YP=0.0d0
          allocate(HP(nij,nij)); HP=0.0d0
 
          ! This is the delta-orbiatl-rotation gradients 
          !          redundant=1.0d-9
          ij=0
          kl=0
          do i=1,ndim1
            do j=1,ndim1
               ij=ij+1
               Y(ij)= 1.0d0*(mat2%A(i,j)-mat2%A(j,i))
               ! Check the valid rotations  
               Lvalid=.false.
               do k=1,nij
                 if(i.eq.valid(k,1).and.j.eq.valid(k,2))then
                   Lvalid=.true.
                   exit
                 end if
               end do
               if(Lvalid)then
                 kl=kl+1
                 porb(kl)=ij
                 !write(*,*)"ij",ij,"kl",kl
               end if
            end do
          end do
   
          do i=1,nij
            YP(i)=Y(porb(i))
            do j=1,nij
              HP(i,j)=H(porb(i),porb(j))
            end do
          end do

          if(.true.)then  ! Here use modified Davidson
            !YP=-1.0d0*YP
            allocate(HPA(nij+1,nij+1))
            HPA=0.0d0

            ! The new augmented Hessian matrix {except HPA(1,1)}
            HPA(2:nij+1,1)=YP
            HPA(1,2:nij+1)=YP
            HPA(2:nij+1,2:nij+1)=HP
            ! And the augmented x-vector
            allocate(XP(nij+1))
            allocate(D(nij+1))
            XP=0.0d0
             D=0.0d0
            !write(*,*)" Before Davidson "

            call Davidson(nij+1,HPA,XP,D,lamda,.true.,100,1.0d-12)
      
            write(1,*)XP(1),"initial XP(1) in loop -",iloop
            !            call print_mat(1,nij+1,XP,6)
            !            if(iloop.eq.1)then

            if(XP(1).lt.0)then
              XP=-1.0d0*XP
            end if    
            if(XP1last.gt.dabs(XP(1)))then
              phase=-1.0d0
            else 
              phase= 1.0d0
            end if     
            if(phaselast/phase.lt.0.0d0)then 
            !              levelshift=0.9+levelshift/10.0d0
            end if
            XP1last=dabs(XP(1)) ! XP last step to ensure a valid direction
            phaselast=phase

            !            end if
            if(dabs(XP(1)).gt.0.99999999d0)then
              Lfinish=.true.
            !              exit
            end if
            if(dabs(XP(1)).gt.0.999999d0)then  ! Currently work good for me 
              nsmooth=nsmooth+1
            !              levelshift=0.99999 
            end if 
            if(iloop.ne.1)then
              do 
                call MODIFIED_GRAM_SCHMIDT(nij+1,1,XP)
                if(XP(1).lt.0)then
                  XP=-1.0d0*XP
                end if
                if(dabs(XP(1)).lt.levelshift)then
                  XP(1)=XP(1)*1.25d0
                  XP(2:nij+1)=XP(2:nij+1)/1.25d0
                else
                 !levelshift=dabs(XP(1)) ! reset levelshift
                  exit 
                end if 
              end do
            end if
            if(phase.lt.0)then
              XP=XP*phase/2.0d0
            end if

            !            write(*,*)iloop," - level_shift ",levelshift
            !            write(*,*)"YP in loop", iloop
            !            call print_mat(1,nij,YP,6)
            !            write(*,*)"XP in loop",iloop
            !            call print_mat(1,nij+1,XP,6)

            do i=1,nij
              X(porb(i))=XP(i+1)
            end do


              ij=0;
              d0=0.0d0; d2=0.0d0;
              d0_sum=0.0d0; d2_sum=0.0d0
              do i=1,nij
                i1=valid(i,1)
                j1=valid(i,2)
                ! recover the anti-symmetric rotations
                mat2%deltR(i1,j1)=  XP(i+1)/1.0d0
                mat2%deltR(j1,i1)= -XP(i+1)/1.0d0
                write(1,*)"rotation R in",iloop,"micro-iter",&
                          i1,j1,mat2%deltR(i1,j1)
                if(dabs(mat2%deltR(i1,j1)).gt.0.10)then
                  write(2,*)"notice : large rotation between",i1,j1,&
                            " : ",mat2%deltR(i1,j1)
                end if
                d0=d0+dabs(XP(i+1))
                d2=d2+(XP(i+1))**2
              end do
              d0_sum = d0_sum + d0
              d2_sum = d2_sum + d2

            
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


              !            ij=0; 
              !            d0=0.0d0; d2=0.0d0; 
              !            d0_sum=0.0d0; d2_sum=0.0d0
              !            do i=1,ndim1
              !              do j=1,ndim1
              !                 ij=ij+1
              !                 mat2%deltR(i,j)=X(ij)
              !                 d0=d0+dabs(X(ij))
              !                 d2=d2+(X(ij))**2
              !              end do
              !            end do
              !            d0_sum = d0_sum + d0/2.0d0
              !            d2_sum = d2_sum + d2/2.0d0
 
            call print_mat(norb,norb,mat2%deltR,6)

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
 
            deallocate(porb)
            deallocate(HPA)
            deallocate(HP)
            deallocate(YP)
            deallocate(XP)
            deallocate(D)

            d0_sum0=0.0d0
            d2_sum0=0.0d0
            do i=1,norb
              do j=1,norb
                d0_sum0=d0_sum0 + dabs(mat2%T(i,j))
                d2_sum0=d2_sum0 + mat2%T(i,j)**2
              end do
            end do 
           
            write(1,*)"== Current rotation parameters =="            
            write(1,*)"Sum(|deltR|)   ",d0_sum
            write(1,*)"Sum(deltR**2)  ",d2_sum
            write(1,*)"Sum(|T|)       ",d0_sum0
            write(1,*)"Sum(T**2)      ",d2_sum0
            write(1,*)"================================="            

              
            if(icycle.eq.1)then
              minloop=20
              THdeltaR=1.0d-8
            else if(icycle.gt.1.and.icycle.le.3)then
              minloop=10
              THdeltaR=1.0d-8
            else
              minloop=5
              THdeltaR=1.0d-8
            end if
            if(iloop.gt.minloop)then
              if(iloop.gt.MAX_iter)then
                write(1,*)"micro-loop reach MAX_iter",MAX_iter 
                exit 
              end if 
              if(d2_sum.lt.THdeltaR)Then
                write(1,*)"Small sum of deltaR detected"
                exit          
              end if 
              if(d2_sum0.gt.nact*1.0d0)then
                write(1,*)"some elements of T**2 is too large"
                exit 
              end if 
              if(Lfinish)then
            write(1,*)"Tiny orbital rotations, should close to converge"
                exit          
              end if 
            end if
            if(nsmooth.gt.3)then
              write(1,*)"Smooth enough MEP"
              exit
            end if
          else            ! Here should be exact Werner-Mayer
          end if
        
        end do 
        Rabs(icycle)=d0_sum0
        write(*,2435,advance='no')"    ",Rabs(icycle)  
        write(*,2436)"   ",method(icycle)  
2435    format(A4,f8.4)
2436    format(A3,A7)

        ! Recover the A,B,G matrix if changed 
        mat2%A=Aini
        mat2%B=Bini
        mat2%G=Gini

        deallocate(Gini)
        deallocate(valid)
        deallocate(X,Y,G,H)

!        stop 

      End Subroutine Solver_Augmented_Hessian_deltaR

! ========================================================================

      Subroutine Solver_Augmented_Hessian(icycle)

        use global_control
        use matrix

        integer::icycle
 
        integer,allocatable::porb(:)
        double precision d0,dv,redundant
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:) 
        double precision,allocatable::H(:,:),HP(:,:),HPA(:,:)

        double precision,allocatable::T1(:,:),T2(:,:),T3(:,:) 

        ! Used for the augmented Hessian part
        double precision lamda   ! damping parameter that control the level shift

        integer,allocatable::valid(:,:)
        double precision,allocatable::G(:)
        logical Lvalid

        ! Parameter initialization
        lamda=1.0d0 ! "1" means the standard AH 
        ! The larger lamda, the shorter step

        ! initialize of part of Hessian matrix that will be used in this solver

        ndim1=norb

        allocate(X(ndim1**2)); X=0.0d0
        allocate(Y(ndim1**2)); Y=0.0d0
        allocate(G(ndim1**2)); G=0.0d0
        allocate(H(ndim1**2,ndim1**2)); H=0.0d0

 
        iHessian=1
        if(iHessian.eq.0)then 
          call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
        else
          if(.NOT.ALLOCATED(mat2%Dh))then
            allocate(mat2%Dh(norb,norb,norb,norb)) ! use the same form as G
          else
            mat2%Dh=0.0d0 
          end if
          if(.NOT.ALLOCATED(mat2%Y))then
            allocate(mat2%Y(norb,norb,norb,norb)) ! use the same form as G
          else
            mat2%Y=0.0d0 
          end if
          ! Same way as in "Molecular electronic-structure theory" 
          write(2,*)"Dh_gen start"
          call flush(2)
          write(2,*)"orb%nsub" ,  orb%nsub
          call flush(2)
          write(2,*)"orb%act"  ,  orb%act
          call flush(2)
          write(2,*)"orb%total",  orb%total
          call flush(2)
          write(2,*)"nact"     ,     nact
          call flush(2)
          write(2,*)"norb"     ,     norb
          call flush(2)
          call Dh_gen(orb%nsub,orb%act,orb%total,&
                     nact,norb,T,mat2%D,mat2%Dh,orb%grouptable)
          write(2,*)"Dh_gen finish"
          call flush(2)
          call Y_gen(orb%nsub,orb%act,orb%total,&
                     nact,norb,U,mat2%P,mat2%Y,orb%grouptable)
          write(2,*)"Y_gen finish"
          call flush(2)
          call Hessian3(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)
        end if

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

        ! This is the orbiatl gradients 
        redundant=1.0d-9

        !        ij=0
        !        nij=0
        !        do i=1,ndim1
        !          do j=1,ndim1
        !             ij=ij+1
        !             Y(ij)=-1.0d0*(mat2%A(i,j)-mat2%A(j,i))
        !             if(dabs(Y(ij)).lt.redundant)then
        !             else
        !               nij=nij+1
        !             end if
        !             write(6,"(f9.5)",advance='no')Y(ij)
        !          end do
        !          write(6,"(f9.5)")Y(ij)
        !        end do

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


          nij0=nij
          ! re-run it to get the valid rotations index
          allocate(valid(nij,2)); valid=0
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

          write(1,*)"valid rotations"
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,1)
          end do
          write(1,*)" "
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,2)
          end do

        nij=nij0/2

        allocate(porb(nij));   porb=0
        allocate(YP(nij));     YP=0.0d0
        allocate(HP(nij,nij)); HP=0.0d0

        !        ij=0
        !        kl=0
        !        do i=1,ndim1
        !          do j=1,ndim1
        !             ij=ij+1
        !             if(dabs(Y(ij)).lt.redundant)then
        !             else 
        !               kl=kl+1
        !               porb(kl)=ij             
        !             end if
        !          end do
        !        end do

        !        do i=1,nij
        !          YP(i)=Y(porb(i))
        !          do j=1,nij
        !            HP(i,j)=H(porb(i),porb(j))
        !          end do
        !        end do

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
            YP(i)= G(porb(i))
            if(dabs(YP(i)).gt.grad_digit)then  ! ??
              grad_digit=dabs(YP(i))
            end if
            do j=1,nij    ! needed Hessian
              HP(i,j)=H(porb(i),porb(j))
            end do
          end do

          !write(1,*)"The valid orb-orb Hessian: "
          !call print_mat(nij,nij,HP,1)         


        !        if(.false.)then  ! Here use modified Davidson
        if(.true.)then  ! Here use modified Davidson
          !YP=-1.0d0*YP
          allocate(HPA(nij+1,nij+1))          
          HPA=0.0d0
          
          ! The new augmented Hessian matrix {except HPA(1,1)}
          HPA(2:nij+1,1)=YP
          HPA(1,2:nij+1)=YP
          HPA(2:nij+1,2:nij+1)=HP
          ! And the augmented x-vector
          allocate(XP(nij+1))
          allocate(D(nij+1))
          XP=0.0d0
           D=0.0d0
          call Davidson(nij+1,HPA,XP,D,lamda,.true.,100,1.0d-13) 

          if(XP(1).lt.0)then
            XP=-1.0d0*XP
          end if
 
          write(1,*)"AH XP"
          call print_mat(nij,1,XP,1) 
          !          stop

          do i=1,nij
            X(porb(i))=XP(i+1)
          end do 



              ij=0;
              d0=0.0d0; d2=0.0d0;
              d0_sum=0.0d0; d2_sum=0.0d0
              do i=1,nij
                i1=valid(i,1)
                j1=valid(i,2)
                ! recover the anti-symmetric rotations
                mat2%R(i1,j1)=  XP(i+1)/1.0d0
                mat2%R(j1,i1)= -XP(i+1)/1.0d0
                write(1,*)"rotation R in",iloop,"micro-iter",&
                          i1,j1,mat2%R(i1,j1)
                if(dabs(mat2%R(i1,j1)).gt.0.10)then
                  write(2,*)"notice : large rotation between",i1,j1,&
                            " : ",mat2%R(i1,j1)
                  ! reduce the redundant affect 
                  !mat2%R=mat2%R/2.0d0
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
                mat2%R(ii,jj)=mat2%R(ii,jj)/2.0d0
                mat2%R(kk,ll)=mat2%R(ii,jj)*redu%sign_mat(kk,ll)
               write(1,*)ii,jj,mat2%R(ii,jj),kk,ll,mat2%R(kk,ll)
                      end do
                    end do
                  end if
                  js0=js0+orb%total(jsym)
                end do
                is0=is0+orb%total(isym)
              end do



          !          ij=0; d0=0.0d0
          !          do i=1,ndim1
          !            do j=1,ndim1
          !               ij=ij+1
          !                mat2%R(i,j)=X(ij)
          !               d0=d0+dabs(X(ij))
          !            end do
          !          end do
          Rabs(icycle)=d0

          deallocate(HPA)
          deallocate(XP)
          deallocate(D)

        else            ! Here use Full(alternative) inverse

          YP=-1.0d0*YP
          allocate(HPA(nij+1,nij+1))
          HPA=0.0d0

          ! The new augmented Hessian matrix {except HPA(1,1)}
          HPA(2:nij+1,1)=YP
          HPA(1,2:nij+1)=YP
          HPA(2:nij+1,2:nij+1)=HP
          ! And the augmented x-vector
          allocate(XP(nij+1))
          allocate(D(nij+1))
          XP=0.0d0
           D=0.0d0
          allocate(T1(nij+1,nij+1)); T1=0.0d0
        
          T1=HPA
          !          call print_mat(nij,nij,HP,6)

          call eigtql2(nij+1,T1,D) 

          !          call Davidson(nij,T2,XP,D,1.0d0,.false.,8,1.0d-12) 
   
          XP=0.0d0
          do i=1,nij
            do k=1,nij
              do l=1,nij
                if(D(l).gt.1.0d-9)then
                  XP(i+1)=XP(i+1)+T1(i+1,l)*T1(k+1,l+1)*YP(k)/D(l+1)
                else if(D(l).lt.-1.0d-9)then
                  XP(i+1)=XP(i+1)-T1(i+1,l)*T1(k+1,l+1)*YP(k)/D(l+1)
                end if
              end do
            end do
          end do

          !          call LENGTH(XP,nij+1,dv)

          do i=1,nij
            X(porb(i))=-XP(i+1)
            !X(porb(i))=-XP(i+1)/dv
          end do 

          write(1,*)"AH XP ",XP
          !         stop
    
          ij=0; d0=0.0d0
          do i=1,ndim1
            do j=1,ndim1
               ij=ij+1
               mat2%R(i,j)=X(ij)
               d0=d0+dabs(X(ij)) 
            end do
          end do
          Rabs(icycle)=d0

          deallocate(D)
          deallocate(XP)
          deallocate(T1)

        end if
  
        if(.true.)then
          call CAL_T(norb,mat2%R,d3,150,mat2%T)
          call T_TO_U(norb,mat2%T,mat2%U)
        else  ! debugging, delete 
          allocate(T1(norb,norb));T1=0.0d0
          allocate(T2(norb,norb));T2=0.0d0
          call CAL_T(norb,mat2%R,d3,150,T1)
          call MXM(norb,mat2%U,T1,T2)
          mat2%T=mat2%T+T2
          call T_TO_U(norb,mat2%T,mat2%U)
          call print_mat(norb,norb,mat2%U,6)
          deallocate(T1) 
          deallocate(T2) 
        end if
        call print_mat(norb,norb,mat2%R,6)
        write(*,2435,advance='no')"    ",Rabs(icycle)  
        write(*,2436)"   ",method(icycle)  
2435    format(A4,f8.4)
2436    format(A3,A7)

        deallocate(X)
        deallocate(Y,YP)
        deallocate(H,HP)

      end Subroutine Solver_Augmented_Hessian
