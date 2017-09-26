      Subroutine Solver_Hessian_Occupy(ndim1,icycle)

        use global_control
        use matrix

        integer::ndim1,icycle
 
        integer,allocatable::porb(:)
        double precision d0,redundant
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:) 
        double precision,allocatable::H(:,:),HP(:,:)

        double precision,allocatable::T1(:,:),T2(:,:) 

! initialize of part of Hessian matrix that will be used in this solver
        allocate(X(ndim1**2)); X=0.0d0
        allocate(Y(ndim1**2)); Y=0.0d0
        allocate(H(ndim1**2,ndim1**2)); H=0.0d0

        !need to be improved using subspace when constructing part of Hessian
        call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag) 

        i0=0 
        ix=0
        do i=1,orb%nsub           
          j0=0 
          jx=0
          do j=1,orb%nsub

            do ii=1,orb%occ(i)                          
              do jj=1,orb%occ(j)           
                mat2%Aocc(ix+ii,jx+jj)=mat2%A(i0+ii,j0+jj)
              end do 
            end do
            j0=j0+orb%total(j)
            jx=jx+orb%occ(j)
          end do
          i0=i0+orb%total(i)
          ix=ix+orb%occ(i)
        end do
         
!        call print_mat(norb,norb,mat2%A)
!        call print_mat(ndim1,ndim1,mat2%Aocc)

        i0=0 
        ix=0
        do i=1,orb%nsub           
         j0=0 
         jx=0
         do j=1,orb%nsub
          k0=0 
          kx=0
          do k=1,orb%nsub
           l0=0 
           lx=0
           do l=1,orb%nsub
            
            do ii=1,orb%occ(i)                          
             do jj=1,orb%occ(j)
              do kk=1,orb%occ(k)
               do ll=1,orb%occ(l)
                 mat2%Hocc(ix+ii,jx+jj,kx+kk,lx+ll)=&
                    mat2%H(i0+ii,j0+jj,k0+kk,l0+ll)
               end do
              end do
             end do
            end do  
            
            l0=l0+orb%total(l)
            lx=lx+orb%occ(l)
           end do
           k0=k0+orb%total(k)
           kx=kx+orb%occ(k)
          end do
          j0=j0+orb%total(j)
          jx=jx+orb%occ(j)
         end do
         i0=i0+orb%total(i)
         ix=ix+orb%occ(i)
        end do

!        call print_gat(nact,nact,nact,nact,mat2%Hocc,7)
!        stop

        ij=0
        do i=1,ndim1
          do j=1,ndim1
            ij=ij+1
            kl=0
            do k=1,ndim1
              do l=1,ndim1
                kl=kl+1
                H(ij,kl)=mat2%Hocc(i,j,k,l) 
              end do
            end do
          end do
        end do

        ! This is the mcscf Fcok matrix 
        redundant=1.0e-9
        ij=0
        nij=0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             Y(ij)=-1.0d0*(mat2%Aocc(i,j)-mat2%Aocc(j,i))
             if(dabs(Y(ij)).lt.redundant)then
             else
               nij=nij+1
             end if
!             H(ij,ij)=100+ij ! when testing
          end do
        end do

        allocate(porb(nij));   porb=0
        allocate(YP(nij));     YP=0.0d0
        allocate(HP(nij,nij)); HP=0.0d0

        ij=0
        kl=0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             if(dabs(Y(ij)).lt.redundant)then
             else 
               kl=kl+1
               porb(kl)=ij             
             end if
!             H(ij,ij)=100+ij ! when testing
          end do
        end do

        do i=1,nij
          YP(i)=Y(porb(i))
          do j=1,nij
            HP(i,j)=H(porb(i),porb(j))
          end do
        end do

        do i=1,nij
          do j=1,nij
!            write(321,*)i,j,HP(i,j) 
          end do 
        end do
        call flush(6)

        allocate(D(nij));       D=0.0d0
        allocate(XP(nij));     XP=0.0d0
        allocate(T1(nij,nij)); T1=0.0d0
        allocate(T2(nij,nij)); T2=0.0d0

        T1=HP
!       call eigtql2(ndim1**2,T1,D) 
        call eigtql2(nij,T1,D) 

        do i=1,nij
!          write(322,*)"hessian eigval",i,D(i) 
        end do

        XP=0.0d0
        do i=1,nij
          do k=1,nij
            do l=1,nij
              if(D(l).gt.1.0e-9)then
                XP(i)=XP(i)+T1(i,l)*T1(k,l)*YP(k)/D(l)         
              else if(D(l).lt.-1.0e-9)then
                XP(i)=XP(i)-T1(i,l)*T1(k,l)*YP(k)/D(l)         
              end if
            end do
          end do
!          write(323,*)"Value of X",i,XP(i) 
          X(porb(i))=XP(i) 
        end do
  
        ij=0; d0=0.0d0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             mat2%Rocc(i,j)=X(ij)
             d0=d0+dabs(X(ij)) 
          end do
        end do
        Rabs(icycle)=d0
!        call print_mat(nact,nact,mat2%Rocc)
!        stop

        call CAL_T(norb,mat2%R,d3,150,mat2%T)
        call T_TO_U(norb,mat2%T,mat2%U)

!        call print_mat(norb,norb,mat2%U)

        deallocate(T1,T2)
        deallocate(X,XP)
        deallocate(Y,YP)
        deallocate(D)
        deallocate(H,HP)

!        stop

      end Subroutine Solver_Hessian_Occupy


      Subroutine Solver_Hessian(ndim1,icycle)

        use global_control
        use matrix

        integer::ndim1,icycle
 
        integer,allocatable::porb(:)
        double precision d0,redundant
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:) 
        double precision,allocatable::H(:,:),HP(:,:)

        double precision,allocatable::T1(:,:),T2(:,:) 

! initialize of part of Hessian matrix that will be used in this solver
        allocate(X(ndim1**2)); X=0.0d0
        allocate(Y(ndim1**2)); Y=0.0d0
        allocate(H(ndim1**2,ndim1**2)); H=0.0d0

        call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)

        do i=1,norb
          do j=1,norb
            do k=1,norb
              do l=1,norb
!                write(31,*)i,j,k,l,mat2%H(i,j,k,l) 
              end do
            end do
          end do
        end do
 
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

        ! This is the mcscf Fock matrix 
        redundant=1.0e-9
        ij=0
        nij=0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             Y(ij)=-1.0d0*(mat2%A(i,j)-mat2%A(j,i))
             if(dabs(Y(ij)).lt.redundant)then
             else
               nij=nij+1
             end if
!             H(ij,ij)=100+ij ! when testing
          end do
        end do

        allocate(porb(nij));   porb=0
        allocate(YP(nij));     YP=0.0d0
        allocate(HP(nij,nij)); HP=0.0d0

        ij=0
        kl=0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             if(dabs(Y(ij)).lt.redundant)then
             else 
               kl=kl+1
               porb(kl)=ij             
             end if
!             H(ij,ij)=100+ij ! when testing
          end do
        end do

        do i=1,nij
          YP(i)=Y(porb(i))
          do j=1,nij
            HP(i,j)=H(porb(i),porb(j))
          end do
        end do

        do i=1,nij
          do j=1,nij
!            write(321,*)i,j,HP(i,j) 
          end do 
        end do
        call flush(6)

!        call inv(ndim1**2,H)
      
!        do i=1,ndim1**2
!          do j=1,ndim1**2
!            write(321,*)"hessian -1 ",i,j,H(i,j) 
!            !write(*,*)"hessian diagonal",i,j,H(i,j) 
!          end do 
!        end do
!        call flush(6)
!        stop  

!        X=1.0d0*Y
!        call MXMG(ndim1**2,1,H,X)

!        allocate(T1(ndim1**2,ndim1**2)); T1=0.0d0
!        allocate(T2(ndim1**2,ndim1**2)); T2=0.0d0
!        allocate(D(ndim1**2)); D=0.0d0
        allocate(D(nij));       D=0.0d0
        allocate(XP(nij));     XP=0.0d0
        allocate(T1(nij,nij)); T1=0.0d0
        allocate(T2(nij,nij)); T2=0.0d0

        T1=HP
!       call eigtql2(ndim1**2,T1,D) 
        call eigtql2(nij,T1,D) 

!        write(*,*)"hessian eigval" 
!        do i=1,ndim1**2
        do i=1,nij
!          if(dabs(D(i)).lt.1.0e-12)then
!             D(i)=0.0d0
!          end if
!          write(322,*)"hessian eigval",i,D(i) 
        end do
!        call MdXM(nij,T1,HP,T2)
!        call MXM(nij,T2,T1,HP)

!        do i=1,nij
!          do j=1,nij
!            write(322,*)"hessian -2 ",i,j,HP(i,j) 
            !write(*,*)"hessian diagonal",i,j,H(i,j) 
!          end do 
!        end do

!        call flush(6)      
!        stop 

        XP=0.0d0
        do i=1,nij
          do k=1,nij
            do l=1,nij
              if(D(l).gt.1.0e-9)then
                XP(i)=XP(i)+T1(i,l)*T1(k,l)*YP(k)/D(l)         
              else if(D(l).lt.-1.0e-9)then
                XP(i)=XP(i)-T1(i,l)*T1(k,l)*YP(k)/D(l)         
              end if
            end do
          end do
!          write(323,*)"Value of X",i,XP(i) 
          X(porb(i))=XP(i) 
        end do
  
        ij=0; d0=0.0d0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             mat2%R(i,j)=X(ij)
             d0=d0+dabs(X(ij)) 
          end do
        end do
        Rabs(icycle)=d0

!        call print_mat(norb,norb,mat2%R)
!        stop
!       Now, deal with the remaining rotations

        call CAL_T(norb,mat2%R,d3,150,mat2%T)
        call T_TO_U(norb,mat2%T,mat2%U)

        write(*,2435,advance='no')"    ",Rabs(icycle)  
        write(*,2436)"   ",method(icycle)  
2435    format(A4,f8.4)
2436    format(A3,A7)

        deallocate(T1,T2)
        deallocate(X,XP)
        deallocate(Y,YP)
        deallocate(D)
        deallocate(H,HP)

!        stop

      end Subroutine Solver_Hessian
