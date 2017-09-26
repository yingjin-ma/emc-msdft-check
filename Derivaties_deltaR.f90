!=======================================================================
! REF H.-J. Werner, Adv. Chem. Phys. LXIX, 1 (1987).
! REF H.-J. Werner and P. J. Knowles, J. Chem. Phys. 82, 5053 (1985).
! ======================================================================

      Subroutine B_delta(norb,nsub,tot,occ,GM,RT,deltaB)

        double precision:: Bini(norb,norb)
        double precision::   RT(norb,norb)
        double precision::   GM(norb,norb,norb,norb)

        double precision::deltaB(norb,norb)
        double precision     TM1(norb,norb)

        integer :: nsub
        integer :: tot(nsub)
        integer :: occ(nsub)

        TM1=0.0d0
        deltaB=0.0d0 

!  Generate the delta B-matrix
        noffset=0
        do n=1,nsub; do n1=1,tot(n)

          moffset=0
          do m=1,nsub; do m1=1,tot(m)

              joffset=0
              do j=1,nsub; do j1=1,tot(j)
                koffset=0
                do k=1,nsub; do k1=1,tot(k)
!  Induce for integrals
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset
!  waveB part
                  TM1(ni,mi)=TM1(ni,mi)&
                             +GM(ni,mi,ji,ki)*RT(ki,ji)
                end do
                koffset=koffset+tot(k)
                end do
              end do
              joffset=joffset+tot(j)
              end do
!            end if

          end do
          moffset=moffset+tot(m)
          end do
        end do
        noffset=noffset+tot(n)
        end do

! return the waveB matrix 
        deltaB=TM1


      End Subroutine B_delta

      Subroutine ARplusRA(norb,ndim,nsub,tot,occ,valid,AM,RT,ARRA)

        integer :: norb,nsub
        integer :: tot(nsub)
        integer :: occ(nsub)
        integer :: valid(ndim,2)

        logical    lvalid

        double precision::AM(norb,norb)
        double precision::RT(norb,norb)
        double precision::ARRA(ndim)

        double precision  AA(norb,norb)
        double precision dv1,dv2

!        write(*,*)nsub,tot,occ
!        write(*,*)"valid",valid(:,1)
!        write(*,*)"valid",valid(:,2) 

        do i=1,norb
          do j=1,norb
            AA(i,j)=AM(i,j)+AM(j,i)
          end do 
        end do

        iarra=0
        do i=1,norb
          do j=1,norb

            Lvalid=.false.
            do jk=1,ndim
              if(i.eq.valid(jk,1).and.j.eq.valid(jk,2))then
                Lvalid=.true.
                exit
              end if
            end do
      
            if(Lvalid)then     
              iarra=iarra+1
              dv1=0.0d0 
              dv2=0.0d0 
              do ij=1,norb
                dv1=dv1+AA(i,ij)*RT(ij,j)
                dv2=dv2+RT(i,ij)*AA(ij,j)
              end do
              ARRA(iarra)=dv1+dv2
            end if

          end do
        end do   


      End Subroutine

      Subroutine dB_update(norb,ndim,nsub,tot,occ,valid,Bini,GM,RT,wB,wBT)

        double precision:: Bini(norb,norb)
        double precision::   RT(norb,norb)
        double precision::   GM(norb,norb,norb,norb)

        double precision:: wB(ndim)
        double precision::wBT(ndim)
        double precision TM1(norb,norb)

        integer :: norb,nsub
        integer :: tot(nsub)
        integer :: occ(nsub) 
        integer :: valid(ndim,2)

        logical lvalid 

!        write(*,*)nsub,tot,occ
!        write(*,*)"valid",valid(:,1) 
!        write(*,*)"valid",valid(:,2)

         wB=0.0d0 
        wBT=0.0d0 

! waveB = B 
!        TM1=Bini
!        call print_mat(norb,norb,Bini,6) 
! waveB = waveB + G*T(R)

!  Generate the B-matrix
        noffset=0
        nn1=0 
        nm=0
        do n=1,nsub; do n1=1,tot(n)
          nn1=nn1+1

          moffset=0
          mm1=0
          do m=1,nsub; do m1=1,tot(m)
            mm1=mm1+1         

            Lvalid=.false.
            do jk=1,ndim
              if(nn1.eq.valid(jk,1).and.mm1.eq.valid(jk,2))then
                Lvalid=.true.
                exit
              end if
            end do

            if(Lvalid)then
              ni=n1+noffset
              mi=m1+moffset
              nm=nm+1
!               wB(nm)=Bini(ni,mi)  ! 
!              wBT(nm)=Bini(mi,ni)  ! 
              !write(*,*)"valid??? ",nn1,mm1,ni,mi,Bini(mi,ni) 

              joffset=0
              jj1=0
              do j=1,nsub; do j1=1,tot(j)
                jj1=jj1+1
 
                koffset=0
                kk1=0
                do k=1,nsub; do k1=1,tot(k)
                  kk1=kk1+1
 
                  ! should be all teh R, thus no need redundant check
!                  Lvalid=.false.
!                  do jk1=1,ndim
!                    if(jj1.eq.valid(jk1,1).and.kk1.eq.valid(jk1,2))then
!                    Lvalid=.true.
!                    exit
!                    end if
!                  end do

!                  if(Lvalid)then
!  Induce for integrals
                    ji=j1+joffset
                    ki=k1+koffset
!  wB part
                     wB(nm)= wB(nm)+GM(ni,mi,ji,ki)*RT(ki,ji)!*2.0d0
                    wBT(nm)=wBT(nm)+GM(mi,ni,ji,ki)*RT(ki,ji)!*2.0d0
!                  end if

                end do
                koffset=koffset+tot(k)
                end do
              end do
              joffset=joffset+tot(j)
              end do
            end if

          end do
          moffset=moffset+tot(m)
          end do
        end do
        noffset=noffset+tot(n)
        end do

! return the waveB matrix 

!        call print_mat(norb,norb,waveB,6) 

      End Subroutine dB_update



      Subroutine Y_update(deltaR,Gini,Bini,Y)

        use global_control
        use matrix

        double precision::Y(norb**2)
        double precision::Bini(norb,norb)
        double precision::deltaR(norb,norb)
        double precision::Gini(norb,norb,norb,norb)

        double precision YM(norb,norb)

        double precision AA(norb,norb)
        double precision UdR(norb,norb)
        double precision UB(norb,norb),BU(norb,norb)
        double precision AdR(norb,norb),dRA(norb,norb)

        double precision Bwave(norb,norb)

! update UdR
        UdR=0.0d0
        call MXM(norb,mat2%U,deltaR,UdR)

!        call print_gat(norb,norb,norb,norb,Gini,6)
        !call print_mat(norb,norb,Bini,6)

! ~B = B + G * UdR
        Bwave=Bini
!       call print_mat(norb,norb,Bwave,6)
!       stop
!  Generate the B-matrix
        noffset=0
        do n=1,orb%nsub; do n1=1,orb%total(n)
          moffset=0
          do m=1,orb%nsub; do m1=1,orb%total(m)
            joffset=0
            do j=1,orb%nsub; do j1=1,orb%occ(j)
              koffset=0
              do k=1,orb%nsub; do k1=1,orb%total(k)
!  Induce for integrals
                ni=n1+noffset
                mi=m1+moffset
                ji=j1+joffset
                ki=k1+koffset
!  B part
                Bwave(ni,mi)=Bwave(ni,mi)+Gini(ni,mi,ji,ki)*UdR(ki,ji)
              end do
              koffset=koffset+orb%total(k)
              end do
            end do
            joffset=joffset+orb%total(j)
            end do
          end do
          moffset=moffset+orb%total(m)
          end do
        end do
        noffset=noffset+orb%total(n)
        end do

        call MdXM(norb,mat2%U,Bwave,UB) 
        call MdXM(norb,Bwave,mat2%U,BU)

        AA=0.0d0
        do i=1,ndim1
          do j=1,ndim1
            AA(i,j)=mat2%A(i,j)+mat2%A(j,i)
          end do
        end do
        
        call MXM(norb,AA,deltaR,AdR) 
        call MXM(norb,deltaR,AA,dRA)

        YM=UB-BU-0.5d0*(AdR+dRA)

!        write(*,*)"The Ym in matrix form"
!        call print_mat(norb,norb,YM,6)
        ij=0
        do i=1,norb
          do j=1,norb
            ij=ij+1
            Y(ij)=YM(i,j)
          end do
        end do 
 
      End Subroutine Y_update

! Only update the G-operator

      Subroutine G_update(Gini)

        use global_control
        use matrix

        double precision::Gini(norb,norb,norb,norb)
        double precision,allocatable::Gd(:,:,:,:)

! ~G^ij = U^+ G^ij U
        allocate(Gd(norb,norb,norb,norb)) ;  Gd=0.0d0
        call UdGU(orb%nsub,orb%act,orb%total,&
                  nact,norb,MAT2%U,Gini,Gd,orb%grouptable)
        mat2%G=Gd
        deallocate(Gd)

      End Subroutine G_update

      Subroutine Derivatives_WMK(icycle,iloop,Gini,Bini)

        use global_control
        use matrix
 
        integer::icycle,iloop
        double precision::Bini(norb,norb)
        double precision::Gini(norb,norb,norb,norb)
 
        double precision UdR(norb,norb)
        double precision,allocatable::Gd(:,:,:,:)

! ~A = U^+ B
        call MdXM(norb,mat2%U,Bini,mat2%A)

! update UdR
        UdR=0.0d0
        call MXM(norb,mat2%U,mat2%deltR,UdR)        

! ~B = B + G * UdR
        mat2%B=Bini
!  Generate the B-matrix
        noffset=0
        do n=1,orb%nsub; do n1=1,orb%total(n)
          moffset=0
          do m=1,orb%nsub; do m1=1,orb%total(m)
            joffset=0
            do j=1,orb%nsub; do j1=1,orb%occ(j)
              koffset=0
              do k=1,orb%nsub; do k1=1,orb%total(k)
!  Induce for integrals
                ni=n1+noffset
                mi=m1+moffset
                ji=j1+joffset
                ki=k1+koffset
!  B part   ! 0.5d0 due to the G is actually 2*G
                mat2%B(ni,mi)=mat2%B(ni,mi)&
                           +Gini(ni,mi,ji,ki)*UdR(ki,ji)*0.5d0  
              end do
              koffset=koffset+orb%total(k)
              end do
            end do
            joffset=joffset+orb%total(j)
            end do
          end do
          moffset=moffset+orb%total(m)
          end do
        end do
        noffset=noffset+orb%total(n)
        end do        

! ~G^ij = U^+ G^ij U
        allocate(Gd(norb,norb,norb,norb)) ;  Gd=0.0d0
        call UdGU(orb%nsub,orb%act,orb%total,&
                  nact,norb,MAT2%U,Gini,Gd,orb%grouptable)
        mat2%G=Gd
        deallocate(Gd)

      End subroutine Derivatives_WMK

      Subroutine Derivaties_deltaR(icycle,iloop,Gini,Bini)

        use global_control
        use matrix
 
        integer::icycle,iloop
        double precision::Bini(norb,norb)
        double precision::Gini(norb,norb,norb,norb)

        double precision,allocatable::Gd(:,:,:,:)

! ~A = U^+ B
        call MdXM(norb,mat2%U,Bini,mat2%A)
! ~G^ij = U^+ G^ij U
        allocate(Gd(norb,norb,norb,norb)) ;  Gd=0.0d0
        call UdGU(orb%nsub,orb%act,orb%total,&
                  nact,norb,MAT2%U,Gini,Gd,orb%grouptable)
        mat2%G=Gd
        deallocate(Gd) 

!        write(*,*)"pass the UdGU " 
!        call flush(6) 

      end subroutine Derivaties_deltaR

!=======================================================================

! G-matrix transformation for deltaR
!      ~G^ij= U+ G^ij U
!  only for integrals in active space
      Subroutine UdGU(nsub,act,tot,nact,norb,TRANS,G,Gact,group)

        INTEGER::norb
        DOUBLE PRECISION::TRANS(norb,norb)
        DOUBLE PRECISION::G(norb,norb,norb,norb)

        integer::act(nsub)
        integer::tot(nsub)
        integer::group(8,8)

        type::transform
          double precision,allocatable::U(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TG1(:,:,:,:)

        double precision::Gact(norb,norb,norb,norb)
        Gact=0.0d0

! Distribute the transform matrix base on group
        allocate(UM(nsub))
        i0=0
        do i=1,nsub
          allocate(UM(i)%U(tot(i),tot(i)))
          UM(i)%U=0.0d0
          UM(i)%U=TRANS(i0+1:i0+tot(i),i0+1:i0+tot(i))
          i0=i0+tot(i)
        end do
!        write(*,*)"pass the U matrix"
!        call flush(6)

        ALLOCATE(TG1(norb,norb,norb,norb));TG1=0.0d0
        ! base on the integrals (index) 
        ! --> (U^+ G^ij U) , ij only active
        i0=0
        DO i=1,nsub
          j0=0
          DO j=1,nsub
            k0=0
            DO k=1,nsub
              l0=0
              DO l=1,nsub
                if(group(i,group(j,group(k,l))).eq.1)then
                  do jj=1,act(j)
                    do kk=1,act(k)

                      allocate(TM1(tot(i),tot(l)));TM1=0.0d0
                      allocate(TM2(tot(i),tot(l)));TM2=0.0d0
                      allocate(TM3(tot(i),tot(l)));TM3=0.0d0

                      TM1=G(i0+1:i0+tot(i),j0+jj,k0+kk,l0+1:l0+tot(l))

                    call MXMG(tot(i),tot(l),tot(i),UM(i)%U,TM1,TM2,"TN")
                    call MXMG(tot(i),tot(l),tot(l),TM2,UM(l)%U,TM3,"NN")

                    TG1(i0+1:i0+tot(i),J0+jj,k0+kk,l0+1:l0+tot(l))&
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

!        write(*,*)"pass the trans matrix"
!        call flush(6)

        Gact=TG1

        deallocate(UM)
        deallocate(TG1)

      End Subroutine UdGU

