! 2-index integrals transformation -- work in micro-iter.
!      h_ij=(U^+hU)_ij
!   (ij|kl)=-(ij|kl)+(U^+ J^kl U)_ij + (U^+ J^ij U)_kl + (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl
!  only for integrals in active space
      Subroutine U2nd_integrals&
                 (nsub,act,tot,nact,norb,U,UJ,UK,Uact,group)

        INTEGER::norb
        DOUBLE PRECISION:: U(norb,norb,norb,norb)
        DOUBLE PRECISION::UJ(norb,norb,norb,norb)
        DOUBLE PRECISION::UK(norb,norb,norb,norb)

        integer::act(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
! For J part
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)
! For K part
        double precision,allocatable::TG3(:,:,:,:),TG4(:,:,:,:)
        double precision,allocatable::TG5(:,:,:,:),TG6(:,:,:,:)
! Actually, only 1 is enough (for checking I use 6)

        double precision::Uact(nact,nact,nact,nact)

! For 2-intergals, two index means accurate to second order in T       
!     1) The -(ij|kl) part
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                   Uact(ia+1:ia+act(i),ja+1:ja+act(j),&
                        ka+1:ka+act(k),la+1:la+act(l))&
                     =U(i0+1:i0+act(i),j0+1:j0+act(j),&
                        k0+1:k0+act(k),l0+1:l0+act(l))
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO
! Change the sign 
        Uact=-1.0d0*Uact         

!     2) (U^+ J^kl U)_ij  and  (U^+ J^ij U)_kl part 
!       basing on <i|J^kl|j> = (ij|kl)
        ALLOCATE(TG1(nact,nact,nact,nact));TG1=0.0d0
        ALLOCATE(TG2(nact,nact,nact,nact));TG2=0.0d0
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                   ! IJKL 
                   TG1(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UJ(i0+1:i0+act(i),j0+1:j0+act(j),&
                       k0+1:k0+act(k),l0+1:l0+act(l))
                   ! KLIJ
                   TG2(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UJ(k0+1:k0+act(k),l0+1:l0+act(l),&
                       i0+1:i0+act(i),j0+1:j0+act(j))
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

!     3) (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl part  
!       basing on <i|K^kl|j> = (ik|lj) 
        ALLOCATE(TG3(nact,nact,nact,nact));TG3=0.0d0
        ALLOCATE(TG4(nact,nact,nact,nact));TG4=0.0d0
        ALLOCATE(TG5(nact,nact,nact,nact));TG5=0.0d0
        ALLOCATE(TG6(nact,nact,nact,nact));TG6=0.0d0

        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                   ! JIKL 
                   TG3(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UK(j0+1:j0+act(j),i0+1:i0+act(i),&
                       k0+1:k0+act(k),l0+1:l0+act(l))
                   ! IJKL
                   TG4(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UK(i0+1:i0+act(i),j0+1:j0+act(j),&
                       k0+1:k0+act(k),l0+1:l0+act(l))
                   ! JILK
                   TG5(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UK(j0+1:j0+act(j),i0+1:i0+act(i),&
                       l0+1:l0+act(l),k0+1:k0+act(k))
                   ! IJLK
                   TG6(ia+1:ia+act(i),ja+1:ja+act(j),&
                       ka+1:ka+act(k),la+1:la+act(l))&
                   =UK(i0+1:i0+act(i),j0+1:j0+act(j),&
                       l0+1:l0+act(l),k0+1:k0+act(k))
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

        Uact=Uact+TG1+TG2+TG3+TG4+TG5+TG6        

        deallocate(TG1)
        deallocate(TG2)
        deallocate(TG3)
        deallocate(TG4)
        deallocate(TG5)
        deallocate(TG6)


      End Subroutine U2nd_integrals 

! transform the integrals base on book "molecular electronic theory"
! formula 12.5.24 & 12.5.25

      Subroutine transform_1index_kappa&
                 (nsub,act,tot,nact,norb,Trans,T,U,T1,U1,group)

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb) ! Trans is equal to kappa
        DOUBLE PRECISION::T(norb,norb)
        DOUBLE PRECISION::U(norb,norb,norb,norb)

        integer::act(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)

        type::transform
          double precision,allocatable::U(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision::T1(norb,norb)
        DOUBLE PRECISION::U1(norb,norb,norb,norb)

        integer o

        write(1,*)"entering transform_1index_kappa"
        call flush(1)

! Distribute the transform matrix base on group
!        allocate(UM(nsub))
!        i0=0
!        do i=1,nsub
!          allocate(UM(i)%U(tot(i),tot(i)))
!          UM(i)%U=0.0d0
!          UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
!!          call print_mat(tot(i),tot(i),UM(i)%U)
!          i0=i0+tot(i)
!        end do       

        T1=0.0d0
        U1=0.0d0 

        write(1,*)"After U=trans"
        call flush(1)

! 1 index transformation
          if(.true.)then ! just for work

            do i=1,norb
              do j=1,norb
                do k=1,norb 
                  T1(i,j)=T1(i,j)+TRANS(i,k)*T(k,j)
                  T1(i,j)=T1(i,j)+TRANS(j,k)*T(i,k)
                end do
              end do
            end do

            do i=1,norb
              do j=1,norb
                do k=1,norb
                  do l=1,norb
                    do o=1,norb
                      U1(i,j,k,l)=U1(i,j,k,l)+TRANS(i,o)*U(o,j,k,l)
                      U1(i,j,k,l)=U1(i,j,k,l)+TRANS(j,o)*U(i,o,k,l)
                      U1(i,j,k,l)=U1(i,j,k,l)+TRANS(k,o)*U(i,j,o,l)
                      U1(i,j,k,l)=U1(i,j,k,l)+TRANS(l,o)*U(i,j,k,o)
                    end do
                  end do
                end do
              end do
            end do
 
          else  ! Subgroup and symmetry 

          end if


      End Subroutine transform_1index_kappa

! 1-index integrals transformation for operators
!     h_rj=(hU)_rj
!  J^kl_rj=(J^kl U)_rj
!  K^kl_rj=(h^kl U)_rj
!  r loops all orbitals, others only internal
!(delete later, since Operator_update do the same thing) // 

      Subroutine transform_1index&
                 (nsub,act,tot,nact,norb,TRANS,T,U,Tact,Uact,group)

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb)
        DOUBLE PRECISION::T(norb,norb)
        DOUBLE PRECISION::U(norb,norb,norb,norb)

        integer::act(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        type::transform
          double precision,allocatable::U(:,:)
          double precision,allocatable::T(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
! For J and K part
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)

! I use the same matrix as these in 2-index for new T and new U, but different dimension
        double precision::Tact(norb,norb)
        double precision::Uact(norb,norb,norb,norb)

          allocate(TM(norb,norb));   TM=0.0d0

! Distribute the transform matrix base on group
        allocate(UM(nsub))
        i0=0
        do i=1,nsub
          allocate(UM(i)%U(tot(i),tot(i)))
          allocate(UM(i)%T(tot(i),tot(i)))
          UM(i)%U=0.0d0
          UM(i)%T=0.0d0
          UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
          UM(i)%T=UM(i)%U
!         T=U-1 in sub-irreps
          do j=1,tot(i)
            UM(i)%T(j,j)=UM(i)%U(j,j)-1.0d0
          end do
!          call print_mat(tot(i),tot(i),UM(i)%U)
          i0=i0+tot(i)
        end do

! For 1-intergals, half index transformation
!       h_rj=(hU)_rj
        i0=0
        do i=1,nsub
          allocate(TM1(tot(i),tot(i)));TM1=0.0d0
          allocate(TM2(tot(i),tot(i)));TM2=0.0d0
          TM2=T(i0+1:i0+tot(i),i0+1:i0+tot(i))
!          call MdXM(tot(i),UM(i)%U,TM1,TM2)
          call  MXM(tot(i),TM2,UM(i)%U,TM1)
          TM(i0+1:i0+tot(i),i0+1:i0+tot(i))=TM1
          deallocate(TM1)
          deallocate(TM2)
          i0=i0+tot(i)
        end do
        Tact=TM
        deallocate(TM)
!        call print_mat(norb,norb,Tact,6)

! For 2-intergals, one index means only affect the non-linear equation later
!     2) (U^+ J^kl U)_rj  and  (T^+ K^kl T)_rj part 
!       basing on <i|J^kl|j> = (ij|kl)
!       basing on <i|K^kl|j> = (ik|lj) 
        ALLOCATE(TG1(norb,norb,norb,norb));TG1=U
        ALLOCATE(TG2(norb,norb,norb,norb));TG2=U
        ! --> (U^+ J^kl U)_rj
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do ii=1,tot(i)
                    do jj=1,act(j) 

!                      allocate(TM1(tot(k),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(k),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(k),tot(l)));TM3=0.0d0

                      TM2=U(I0+ii,J0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))

!                      write(*,*)"i,j,k,l",i,j,k,l 
!                    call MXMG(tot(k),tot(l),tot(k),UM(k)%U,TM1,TM2,"TN")
                    call MXMG(tot(k),tot(l),tot(l),TM2,UM(l)%U,TM3,"NN")

                    TG1(I0+ii,J0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))&
                               =TM3(1:tot(k),1:tot(l))

!                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO

!        goto 1000

        ! --> (T^+ K^kl T)_rj   ! modified  T -> U , diff to the reference
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do ii=1,tot(i)
                    do jj=1,act(j)

                      !allocate(TM1(tot(i),tot(k)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(k)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(k)));TM3=0.0d0

                      TM2=TG1(i0+ii,k0+1:k0+tot(k),l0+1:l0+tot(l),j0+jj)

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    !call MXMG(tot(i),tot(k),tot(i),UM(i)%T,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(k),tot(k),TM2,UM(k)%U,TM3,"NN")

                    TG2(i0+ii,k0+1:k0+tot(k),l0+1:l0+tot(l),j0+jj)&
                               =TM3(1:tot(i),1:tot(k))

                      !deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
      
! transformed integrals to Uact 
        Uact=TG2 

!1000    continue  !(only transform the J part)
        Uact=TG1
!        call print_gat(norb,norb,norb,norb,Uact,6)

        deallocate(UM)
        deallocate(TG1)
        deallocate(TG2)

      End subroutine transform_1index

! 2-index integrals transformation
!      h_ij=(U^+hU)_ij
!   (ij|kl)=-(ij|kl)+(U^+ J^kl U)_ij + (U^+ J^ij U)_kl + (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl
!  only for integrals in active space
      Subroutine transform_2index&
                 (nsub,act,tot,nact,norb,TRANS,T,U,Tact,Uact,group)

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb)
        DOUBLE PRECISION::T(norb,norb)
        DOUBLE PRECISION::U(norb,norb,norb,norb)

        integer::act(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        type::transform
          double precision,allocatable::U(:,:)
          double precision,allocatable::T(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
! For J part
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)
! For K part
        double precision,allocatable::TG3(:,:,:,:),TG4(:,:,:,:)
        double precision,allocatable::TG5(:,:,:,:),TG6(:,:,:,:)
! Actually, only 1 is enough (for checking I use 6)

        double precision::Tact(nact,nact) 
        double precision::Uact(nact,nact,nact,nact) 

          allocate(TM(norb,norb));   TM=0.0d0       
! Active Integrals for T and U(V) 
!        allocate(Tact(nact,nact)); Tact=0.0d0     
!        allocate(Uact(nact,nact,nact,nact));  Uact=0.0d0   
!        call print_mat(norb,norb,TRANS,6)

! Distribute the transform matrix base on group
        allocate(UM(nsub))
        i0=0
        do i=1,nsub
          allocate(UM(i)%U(tot(i),tot(i)))
          allocate(UM(i)%T(tot(i),tot(i)))
          UM(i)%U=0.0d0
          UM(i)%T=0.0d0
          UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
          UM(i)%T=UM(i)%U 
!         T=U-1 in sub-irreps
          do j=1,tot(i)
            UM(i)%T(j,j)=UM(i)%U(j,j)-1.0d0
          end do
!          call print_mat(tot(i),tot(i),UM(i)%U)
          i0=i0+tot(i)
        end do

! For 1-intergals, two index means fully transformed
        i0=0
        do i=1,nsub
          allocate(TM1(tot(i),tot(i)));TM1=0.0d0
          allocate(TM2(tot(i),tot(i)));TM2=0.0d0
          TM1=T(i0+1:i0+tot(i),i0+1:i0+tot(i))
          call MdXM(tot(i),UM(i)%U,TM1,TM2)
          call  MXM(tot(i),TM2,UM(i)%U,TM1)
          TM(i0+1:i0+tot(i),i0+1:i0+tot(i))=TM1
          deallocate(TM1)
          deallocate(TM2)
          i0=i0+tot(i)
        end do
! The transformed 1-body integrals (only active to Tact)
        i0=0  ! offset for   full 1-integrals
        ia=0  ! offset for active 1-integrals
        do i=1,nsub
          Tact(ia+1:ia+act(i),ia+1:ia+act(i))&
           =TM(i0+1:i0+act(i),i0+1:i0+act(i))
          i0=i0+tot(i)
          ia=ia+act(i)
        end do 
        deallocate(TM) 
!       call print_mat(nact,nact,Tact,6)
!       call flush(6) 
        
! For 2-intergals, two index means accurate to second order in T       
!     1) The -(ij|kl) part
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0 
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                   Uact(ia+1:ia+act(i),ja+1:ja+act(j),&
                        ka+1:ka+act(k),la+1:la+act(l))&
                     =U(i0+1:i0+act(i),j0+1:j0+act(j),&
                        k0+1:k0+act(k),l0+1:l0+act(l))
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO              
! Change the sign 
        Uact=-1.0d0*Uact         

!        call print_gat(nact,nact,nact,nact,Uact,11) 
!        write(*,*)1,2,3,4,Uact(1,2,3,4)
!        write(*,*)2,1,3,4,Uact(2,1,3,4)
!        write(*,*)1,2,4,3,Uact(1,2,4,3)
!        write(*,*)2,1,4,3,Uact(2,1,4,3)
!        stop 
      
!     2) (U^+ J^kl U)_ij  and  (U^+ J^ij U)_kl part 
!       basing on <i|J^kl|j> = (ij|kl)
        ALLOCATE(TG1(nact,nact,nact,nact));TG1=0.0d0
        ALLOCATE(TG2(nact,nact,nact,nact));TG2=0.0d0

        ! --> (U^+ J^kl U)_ij          
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do kk=1,act(k)
                    do ll=1,act(l)

                      allocate(TM1(tot(i),tot(j)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(j)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(j)));TM3=0.0d0

                      TM1=U(i0+1:i0+tot(i),j0+1:j0+tot(j),k0+kk,l0+ll)

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(j),tot(i),UM(i)%U,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(j),tot(j),TM2,UM(j)%U,TM3,"NN")

                    TG2(ia+1:ia+act(i),ja+1:ja+act(j),ka+kk,la+ll)&
                               =TM3(1:act(i),1:act(j))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO        

        ! --> (U^+ J^ij U)_kl
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do ii=1,act(i)
                    do jj=1,act(j)

                      allocate(TM1(tot(k),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(k),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(k),tot(l)));TM3=0.0d0

                      TM1=U(I0+ii,J0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(k),tot(l),tot(k),UM(k)%U,TM1,TM2,"TN")
                    call MXMG(tot(k),tot(l),tot(l),TM2,UM(l)%U,TM3,"NN")

                    TG1(Ia+ii,Ja+jj,ka+1:ka+act(k),la+1:la+act(l))&
                               =TM3(1:act(k),1:act(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

!        write(6,*)"Finish 2-index transformation on J"
!        call flush(6)

!     3) (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl part  
!       basing on <i|K^kl|j> = (ik|lj) 
        ALLOCATE(TG3(nact,nact,nact,nact));TG3=0.0d0
        ALLOCATE(TG4(nact,nact,nact,nact));TG4=0.0d0
        ALLOCATE(TG5(nact,nact,nact,nact));TG5=0.0d0
        ALLOCATE(TG6(nact,nact,nact,nact));TG6=0.0d0

!        goto 2222 

        ! --> (T^+ K^ik T)_jl
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do ii=1,act(i)
                    do kk=1,act(k)

!                      allocate(TM1(tot(l),tot(j)));TM1=0.0d0
!                      allocate(TM2(tot(l),tot(j)));TM2=0.0d0
!                      allocate(TM3(tot(l),tot(j)));TM3=0.0d0
                      allocate(TM1(tot(j),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(j),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(j),tot(l)));TM3=0.0d0

!                      TM1=U(i0+ii,k0+kk,l0+1:l0+tot(l),j0+1:j0+tot(j))
                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+kk,l0+1:l0+tot(l)) ! (ijkl)
!                      TM1=U(i0+ii,j0+1:j0+tot(j),l0+1:l0+tot(l),k0+kk) ! Ron? 

!                      write(*,*)"i,j,k,l",i,j,k,l 

!                    call MXMG(tot(l),tot(j),tot(l),UM(l)%T,TM1,TM2,"TN")
!                    call MXMG(tot(l),tot(j),tot(j),TM2,UM(j)%T,TM3,"NN")
                    call MXMG(tot(j),tot(l),tot(j),UM(j)%T,TM1,TM2,"TN")
                    call MXMG(tot(j),tot(l),tot(l),TM2,UM(l)%T,TM3,"NN")

!                    TG3(ia+ii,ka+kk,la+1:la+act(l),ja+1:ja+act(j))&
                    TG3(ia+ii,ja+1:ja+act(j),ka+kk,la+1:la+act(l))&
!                               =TM3(1:act(l),1:act(j))
                               =TM3(1:act(j),1:act(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

        ! --> (T^+ K^jk T)_il
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do jj=1,act(j)
                    do kk=1,act(k)

                      allocate(TM1(tot(i),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(l)));TM3=0.0d0

!                      TM1=U(i0+1:i0+tot(i),k0+kk,l0+1:l0+tot(l),j0+jj)
                      TM1=U(i0+1:i0+tot(i),j0+jj,k0+kk,l0+1:l0+tot(l)) ! (ijkl)
!                      TM1=U(j0+jj,i0+1:i0+tot(i),l0+1:l0+tot(l),k0+kk) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(l),tot(i),UM(i)%T,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(l),tot(l),TM2,UM(l)%T,TM3,"NN")

!                    TG4(ia+1:ia+act(i),ka+kk,la+1:la+act(l),ja+jj)&
                    TG4(ia+1:ia+act(i),ja+jj,ka+kk,la+1:la+act(l))&
                               =TM3(1:act(i),1:act(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

        ! --> (T^+ K^il T)_jk
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do ii=1,act(i)
                    do ll=1,act(l)

!                      allocate(TM1(tot(k),tot(j)));TM1=0.0d0
!                      allocate(TM2(tot(k),tot(j)));TM2=0.0d0
!                      allocate(TM3(tot(k),tot(j)));TM3=0.0d0
                      allocate(TM1(tot(j),tot(k)));TM1=0.0d0
                      allocate(TM2(tot(j),tot(k)));TM2=0.0d0
                      allocate(TM3(tot(j),tot(k)));TM3=0.0d0

!                      TM1=U(i0+ii,k0+1:k0+tot(k),l0+ll,j0+1:j0+tot(j))
                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+1:k0+tot(k),l0+ll) ! (ijkl)
!                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+1:k0+tot(k),l0+ll) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

!                    call MXMG(tot(k),tot(j),tot(k),UM(k)%T,TM1,TM2,"TN")
!                    call MXMG(tot(k),tot(j),tot(j),TM2,UM(j)%T,TM3,"NN")
                    call MXMG(tot(j),tot(k),tot(j),UM(j)%T,TM1,TM2,"TN")
                    call MXMG(tot(j),tot(k),tot(k),TM2,UM(k)%T,TM3,"NN")

!                    TG5(ia+ii,ka+1:ka+act(k),la+ll,ja+1:ja+act(j))&
                    TG5(ia+ii,ja+1:ja+act(j),ka+1:ka+act(k),la+ll)&
!                               =TM3(1:act(k),1:act(j))
                               =TM3(1:act(j),1:act(k))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

        ! --> (T^+ K^jl T)_ik
        i0=0
        ia=0
        DO i=1,nsub
          j0=0
          ja=0
          DO j=1,nsub
            k0=0
            ka=0
            DO k=1,nsub
              l0=0
              la=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for active part
                  do jj=1,act(j)
                    do ll=1,act(l)

                      allocate(TM1(tot(i),tot(k)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(k)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(k)));TM3=0.0d0

!                      TM1=U(i0+1:i0+tot(i),k0+1:k0+tot(k),l0+ll,j0+jj) !(iklj)
                      TM1=U(i0+1:i0+tot(i),j0+jj,k0+1:k0+tot(k),l0+ll) ! (ijkl)
!                      TM1=U(j0+jj,i0+1:i0+tot(i),k0+1:k0+tot(k),l0+ll) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(k),tot(i),UM(i)%T,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(k),tot(k),TM2,UM(k)%T,TM3,"NN")

!                    TG6(ia+1:ia+act(i),ka+1:ka+act(k),la+ll,ja+jj)&
                    TG6(ia+1:ia+act(i),ja+jj,ka+1:ka+act(k),la+ll)&
                               =TM3(1:act(i),1:act(k))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
                la=la+act(l)
              END DO
              k0=k0+tot(k)
              ka=ka+act(k)
            END DO
            j0=j0+tot(j)
            ja=ja+act(j)
          END DO
          i0=i0+tot(i)
          ia=ia+act(i)
        END DO

!        write(6,*)"Finish 2-index transformation on K"
!        call flush(6)
!        call print_GAT(nact,nact,nact,nact,Uact,6)
!        write(6,*) "=================||||=====================" 

!2222    continue    

        Uact=Uact+TG1+TG2+TG3+TG4+TG5+TG6
!        call print_GAT(nact,nact,nact,nact,Uact,6)

        deallocate(UM)
        deallocate(TG1)
        deallocate(TG2)
        deallocate(TG3)
        deallocate(TG4)
        deallocate(TG5)
        deallocate(TG6)       


      end subroutine transform_2index

! 2-index integrals transformation
!      h_ij=(U^+hU)_ij
!   (ij|kl)=-(ij|kl)+(U^+ J^kl U)_ij + (U^+ J^ij U)_kl + (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl
!  only for integrals in active space
      Subroutine transform_2index_ALL&
                 (nsub,tot,norb,TRANS,T,U,Ttot,Utot,group)

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb)
        DOUBLE PRECISION::T(norb,norb)
        DOUBLE PRECISION::U(norb,norb,norb,norb)

        integer::tot(nsub)
        integer::group(8,8)

        type::transform
          double precision,allocatable::U(:,:)
          double precision,allocatable::T(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:),TM(:,:)
! For J part
        double precision,allocatable::TG1(:,:,:,:),TG2(:,:,:,:)
! For K part
        double precision,allocatable::TG3(:,:,:,:),TG4(:,:,:,:)
        double precision,allocatable::TG5(:,:,:,:),TG6(:,:,:,:)
! Actually, only 1 is enough (for checking I use 6)

        double precision::Ttot(norb,norb) 
        double precision::Utot(norb,norb,norb,norb) 

        Ttot=0.0d0 
        Utot=0.0d0 

          allocate(TM(norb,norb));   TM=0.0d0       
! Active Integrals for T and U(V) 
!        allocate(Ttot(norb,norb)); Ttot=0.0d0     
!        allocate(Utot(norb,norb,norb,norb));  Utot=0.0d0   
!        call print_mat(norb,norb,TRANS,6)

! Distribute the transform matrix base on group
        allocate(UM(nsub))
        i0=0
        do i=1,nsub
          allocate(UM(i)%U(tot(i),tot(i)))
          allocate(UM(i)%T(tot(i),tot(i)))
          UM(i)%U=0.0d0
          UM(i)%T=0.0d0
          UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
          UM(i)%T=UM(i)%U 
!         T=U-1 in sub-irreps
          do j=1,tot(i)
            UM(i)%T(j,j)=UM(i)%U(j,j)-1.0d0
          end do
!          call print_mat(tot(i),tot(i),UM(i)%U)
          i0=i0+tot(i)
        end do

! For 1-intergals, two index means fully transformed
        i0=0
        do i=1,nsub
          allocate(TM1(tot(i),tot(i)));TM1=0.0d0
          allocate(TM2(tot(i),tot(i)));TM2=0.0d0
          TM1=T(i0+1:i0+tot(i),i0+1:i0+tot(i))
          call MdXM(tot(i),UM(i)%U,TM1,TM2)
          call  MXM(tot(i),TM2,UM(i)%U,TM1)
          TM(i0+1:i0+tot(i),i0+1:i0+tot(i))=TM1
          deallocate(TM1)
          deallocate(TM2)
          i0=i0+tot(i)
        end do
! The transformed 1-body integrals
        i0=0  ! offset for   full 1-integrals
        do i=1,nsub
          Ttot(i0+1:i0+tot(i),i0+1:i0+tot(i))&
           =TM(i0+1:i0+tot(i),i0+1:i0+tot(i))
          i0=i0+tot(i)
        end do 
        deallocate(TM) 
!        call print_mat(norb,norb,Ttot,6)
!        call flush(6) 
        
! For 2-intergals, two index means accurate to second order in T       
!     1) The -(ij|kl) part
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0 
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                   Utot(i0+1:i0+tot(i),j0+1:j0+tot(j),&
                        k0+1:k0+tot(k),l0+1:l0+tot(l))&
                     =U(i0+1:i0+tot(i),j0+1:j0+tot(j),&
                        k0+1:k0+tot(k),l0+1:l0+tot(l))
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO              
! Change the sign 
        Utot=-1.0d0*Utot         

!        call print_gat(norb,norb,norb,norb,Utot,11) 
!        write(*,*)1,2,3,4,Utot(1,2,3,4)
!        write(*,*)2,1,3,4,Utot(2,1,3,4)
!        write(*,*)1,2,4,3,Utot(1,2,4,3)
!        write(*,*)2,1,4,3,Utot(2,1,4,3)
!        stop 
      
!     2) (U^+ J^kl U)_ij  and  (U^+ J^ij U)_kl part 
!       basing on <i|J^kl|j> = (ij|kl)
        ALLOCATE(TG1(norb,norb,norb,norb));TG1=0.0d0
        ALLOCATE(TG2(norb,norb,norb,norb));TG2=0.0d0

        ! --> (U^+ J^kl U)_ij          
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do kk=1,tot(k)
                    do ll=1,tot(l)

                      allocate(TM1(tot(i),tot(j)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(j)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(j)));TM3=0.0d0

                      TM1=U(i0+1:i0+tot(i),j0+1:j0+tot(j),k0+kk,l0+ll)

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(j),tot(i),UM(i)%U,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(j),tot(j),TM2,UM(j)%U,TM3,"NN")

                    TG2(i0+1:i0+tot(i),j0+1:j0+tot(j),k0+kk,l0+ll)&
                               =TM3(1:tot(i),1:tot(j))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO        
!        call print_gat(norb,norb,norb,norb,TG2,12) 

        ! --> (U^+ J^ij U)_kl
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do ii=1,tot(i)
                    do jj=1,tot(j)

                      allocate(TM1(tot(k),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(k),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(k),tot(l)));TM3=0.0d0

                      TM1=U(i0+ii,j0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(k),tot(l),tot(k),UM(k)%U,TM1,TM2,"TN")
                    call MXMG(tot(k),tot(l),tot(l),TM2,UM(l)%U,TM3,"NN")

                    TG1(i0+ii,j0+jj,k0+1:k0+tot(k),l0+1:l0+tot(l))&
                               =TM3(1:tot(k),1:tot(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
!        call print_gat(norb,norb,norb,norb,TG1,11) 

!        write(6,*)"Finish 2-index transformation on J"
!        call flush(6)

!     3) (1+\tau_ij)(1+\tau_kl)(T^+ K^ik T)_jl part  
!       basing on <i|K^kl|j> = (ik|lj) 
        ALLOCATE(TG3(norb,norb,norb,norb));TG3=0.0d0
        ALLOCATE(TG4(norb,norb,norb,norb));TG4=0.0d0
        ALLOCATE(TG5(norb,norb,norb,norb));TG5=0.0d0
        ALLOCATE(TG6(norb,norb,norb,norb));TG6=0.0d0

!        goto 2222 

        ! --> (T^+ K^ik T)_jl
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do ii=1,tot(i)
                    do kk=1,tot(k)

!                      allocate(TM1(tot(l),tot(j)));TM1=0.0d0
!                      allocate(TM2(tot(l),tot(j)));TM2=0.0d0
!                      allocate(TM3(tot(l),tot(j)));TM3=0.0d0
                      allocate(TM1(tot(j),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(j),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(j),tot(l)));TM3=0.0d0

!                      TM1=U(i0+ii,k0+kk,l0+1:l0+tot(l),j0+1:j0+tot(j))
                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+kk,l0+1:l0+tot(l)) ! (ijkl)
!                      TM1=U(i0+ii,j0+1:j0+tot(j),l0+1:l0+tot(l),k0+kk) ! Ron? 

!                      write(*,*)"i,j,k,l",i,j,k,l 

!                    call MXMG(tot(l),tot(j),tot(l),UM(l)%T,TM1,TM2,"TN")
!                    call MXMG(tot(l),tot(j),tot(j),TM2,UM(j)%T,TM3,"NN")
                    call MXMG(tot(j),tot(l),tot(j),UM(j)%T,TM1,TM2,"TN")
                    call MXMG(tot(j),tot(l),tot(l),TM2,UM(l)%T,TM3,"NN")

!                    TG3(i0+ii,k0+kk,l0+1:l0+tot(l),j0+1:j0+tot(j))&
                    TG3(i0+ii,j0+1:j0+tot(j),k0+kk,l0+1:l0+tot(l))&
!                               =TM3(1:tot(l),1:tot(j))
                               =TM3(1:tot(j),1:tot(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
!        call print_gat(norb,norb,norb,norb,TG3,13) 

        ! --> (T^+ K^jk T)_il
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do jj=1,tot(j)
                    do kk=1,tot(k)

                      allocate(TM1(tot(i),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(l)));TM3=0.0d0

!                      TM1=U(i0+1:i0+tot(i),k0+kk,l0+1:l0+tot(l),j0+jj)
                      TM1=U(i0+1:i0+tot(i),j0+jj,k0+kk,l0+1:l0+tot(l)) ! (ijkl)
!                      TM1=U(j0+jj,i0+1:i0+tot(i),l0+1:l0+tot(l),k0+kk) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(l),tot(i),UM(i)%T,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(l),tot(l),TM2,UM(l)%T,TM3,"NN")

!                    TG4(i0+1:i0+tot(i),k0+kk,l0+1:l0+tot(l),j0+jj)&
                    TG4(i0+1:i0+tot(i),j0+jj,k0+kk,l0+1:l0+tot(l))&
                               =TM3(1:tot(i),1:tot(l))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
!        call print_gat(norb,norb,norb,norb,TG2,14) 

        ! --> (T^+ K^il T)_jk
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do ii=1,tot(i)
                    do ll=1,tot(l)

!                      allocate(TM1(tot(k),tot(j)));TM1=0.0d0
!                      allocate(TM2(tot(k),tot(j)));TM2=0.0d0
!                      allocate(TM3(tot(k),tot(j)));TM3=0.0d0
                      allocate(TM1(tot(j),tot(k)));TM1=0.0d0
                      allocate(TM2(tot(j),tot(k)));TM2=0.0d0
                      allocate(TM3(tot(j),tot(k)));TM3=0.0d0

!                      TM1=U(i0+ii,k0+1:k0+tot(k),l0+ll,j0+1:j0+tot(j))
                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+1:k0+tot(k),l0+ll) ! (ijkl)
!                      TM1=U(i0+ii,j0+1:j0+tot(j),k0+1:k0+tot(k),l0+ll) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

!                    call MXMG(tot(k),tot(j),tot(k),UM(k)%T,TM1,TM2,"TN")
!                    call MXMG(tot(k),tot(j),tot(j),TM2,UM(j)%T,TM3,"NN")
                    call MXMG(tot(j),tot(k),tot(j),UM(j)%T,TM1,TM2,"TN")
                    call MXMG(tot(j),tot(k),tot(k),TM2,UM(k)%T,TM3,"NN")

!                    TG5(i0+ii,k0+1:k0+tot(k),l0+ll,j0+1:j0+tot(j))&
                    TG5(i0+ii,j0+1:j0+tot(j),k0+1:k0+tot(k),l0+ll)&
!                               =TM3(1:tot(k),1:tot(j))
                               =TM3(1:tot(j),1:tot(k))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
!        call print_gat(norb,norb,norb,norb,TG5,15) 

        ! --> (T^+ K^jl T)_ik
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  ! only for totive part
                  do jj=1,tot(j)
                    do ll=1,tot(l)

                      allocate(TM1(tot(i),tot(k)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(k)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(k)));TM3=0.0d0

!                      TM1=U(i0+1:i0+tot(i),k0+1:k0+tot(k),l0+ll,j0+jj) !(iklj)
                      TM1=U(i0+1:i0+tot(i),j0+jj,k0+1:k0+tot(k),l0+ll) ! (ijkl)
!                      TM1=U(j0+jj,i0+1:i0+tot(i),k0+1:k0+tot(k),l0+ll) ! Ron

!                      write(*,*)"i,j,k,l",i,j,k,l 

                    call MXMG(tot(i),tot(k),tot(i),UM(i)%T,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(k),tot(k),TM2,UM(k)%T,TM3,"NN")

!                    TG6(i0+1:i0+tot(i),k0+1:k0+tot(k),l0+ll,j0+jj)&
                    TG6(i0+1:i0+tot(i),j0+jj,k0+1:k0+tot(k),l0+ll)&
                               =TM3(1:tot(i),1:tot(k))

                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end do
                  end do
                end if
                l0=l0+tot(l)
              END DO
              k0=k0+tot(k)
            END DO
            j0=j0+tot(j)
          END DO
          i0=i0+tot(i)
        END DO
!        call print_gat(norb,norb,norb,norb,TG6,16) 

!        write(6,*)"Finish 2-index transformation on K"
!        call flush(6)
!        call print_GAT(norb,norb,norb,norb,Utot,6)
!        write(6,*) "=================||||=====================" 

!2222    continue    

        Utot=Utot+TG1+TG2+TG3+TG4+TG5+TG6
!        call print_GAT(norb,norb,norb,norb,Utot,20)

        deallocate(UM)
        deallocate(TG1)
        deallocate(TG2)
        deallocate(TG3)
        deallocate(TG4)
        deallocate(TG5)
        deallocate(TG6)       

      end subroutine transform_2index_ALL


