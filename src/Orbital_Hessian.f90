      Subroutine Hessian3(norb,dh,F,Y,Hess,Hdiag)

        integer::norb,nact
        double precision::    F(norb,norb)
        double precision::   dh(norb,norb,norb,norb)
        double precision::    Y(norb,norb,norb,norb)
        double precision:: Hess(norb,norb,norb,norb)
        double precision::Hdiag(norb,norb)

        double precision dij,drj,dis,drs,adelt
        double precision,allocatable::aplus(:,:)
 
        integer r,s

        allocate(aplus(norb,norb))

        dij=0.0d0
        drj=0.0d0
        dis=0.0d0
        drs=0.0d0
        Hdiag=0.0d0
        Hess =0.0d0
        aplus=0.0d0

        do i=1,norb
          do j=1,norb
            aplus(i,j)=F(i,j)+F(j,i)
          end do
        end do

        write(2,*)" ==================== Warning ===================== "
        write(2,*)"  The Hessian not fully-tested for irreps cases, "
        write(2,*)"      especially for closed-shell approximation. "
        write(2,*)" ==================== Warning ===================== "
        call flush(2)

        ! Hessian_ri,sj  ::  r,s all type, i,j active       
        do r=1,norb
          do i=1,norb    ! only for c1 if nact
            do s=1,norb
              do j=1,norb   ! only for c1 if nact

                adelt=0.0d0
                call delt(i,j,dij)
                call delt(r,j,drj)
                call delt(i,s,dis)
                call delt(r,s,drs)

                adelt=aplus(r,s)*dij-aplus(i,s)*drj-aplus(r,j)*dis+aplus(i,j)*drs   ! all
                hess(r,i,s,j)=hess(r,i,s,j)+dh(r,s,i,j)-dh(i,s,r,j)-dh(r,j,i,s)+dh(i,j,r,s)
                hess(r,i,s,j)=hess(r,i,s,j)+Y(r,i,s,j)-Y(i,r,s,j)-Y(r,i,j,s)+Y(i,r,j,s)
                hess(r,i,s,j)=hess(r,i,s,j)-0.5d0*adelt
              end do
            end do
            Hdiag(r,i)=hess(r,i,r,i)
          end do
        end do

      End Subroutine Hessian3

      Subroutine HessianHF(norb,Dh,F,Y,Hess,Hdiag)

! Use the way of generating HF Hessian, should be same in Molcas's MCLR
!   The version for tesing, not consider effcience

        integer::norb
        double precision::    F(norb,norb)
        double precision::   Dh(norb,norb,norb,norb)
        double precision::    Y(norb,norb,norb,norb)
        double precision:: Hess(norb,norb,norb,norb)
        double precision::Hdiag(norb,norb)

        double precision djl,djk,dil,dik,adelt
        double precision,allocatable::aplus(:,:)
! should be a_dag ...

        allocate(aplus(norb,norb))

        djl=0.0d0
        djk=0.0d0
        dil=0.0d0
        dik=0.0d0
        Hdiag=0.0d0
        Hess =0.0d0
        aplus=0.0d0

        do i=1,norb
          do j=1,norb
            aplus(i,j)=F(i,j)+F(j,i)
            !write(*,*)i,j,aplus(i,j)
          end do
        end do         

!
        do i=1,norb
          do j=1,norb
            do k=1,norb
              do l=1,norb
                adelt=0.0d0
                call delt(j,l,djl)
                call delt(j,k,djk)
                call delt(i,l,dil)
                call delt(i,k,dik)

                adelt=aplus(i,k)*djl-aplus(i,l)*djk-aplus(j,k)*dil&
                     +aplus(j,l)*dik

                hess(i,j,k,l)= Dh(i,k,j,l)-Dh(j,k,i,l) &
                              -Dh(i,l,j,k)+Dh(j,l,i,k)

                hess(i,j,k,l)=hess(i,j,k,l)+Y(i,j,k,l)-Y(j,i,k,l)&
                                           -Y(i,j,l,k)+Y(j,i,l,k)
                hess(i,j,k,l)=hess(i,j,k,l)-0.5d0*adelt

!                write(2323,2324)i,j,k,l,hess(i,j,k,l)

              end do
            end do
            Hdiag(i,j)=hess(i,j,i,j)
          end do
        end do

2324    format(i2,i2,i2,i2,f12.8)
!        write(*,*)"preonescf",iorb
!        stop

        deallocate(aplus)


      End Subroutine HessianHF 

      subroutine Hessian2()

        use global_control
        use matrix 

        double precision djl,djk,dil,dik,adelt
        double precision,allocatable::aplus(:,:)
        
        type::mat_hessian
          double precision,allocatable::aplus(:,:) 
        end type mat_hessian
        type(mat_hessian),allocatable::MT(:) 
              
        allocate(MT(orb%nsub)) 

        do i=1,orb%nsub
          allocate(MT(i)%aplus(orb%total(i),orb%total(i)))
          MT(i)%aplus=0.0d0
          do k=1,orb%total(i)
            do l=1,orb%total(i)
              MT(i)%aplus(k,l)=mat1(i)%A(k,l)+mat2%A(l,k)
            end do
          end do
        end do

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
            if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
              do ii=1,orb%total(i)
              do jj=1,orb%total(j)
                do kk=1,orb%total(k)
                do ll=1,orb%total(l)

                  !i1=
                  ! not finish, since the optimized way will be so tedious 
                  ! and should be on-the-fly when solving, so currently 
                  ! back to the old form
 
                  adelt=0.0d0

                  call delt(j,l,djl)
                  call delt(j,k,djk)
                  call delt(i,l,dil)
                  call delt(i,k,dik)

                  adelt=aplus(i,k)*djl-aplus(i,l)*djk-aplus(j,k)*dil&
                       +aplus(j,l)*dik

!                  hmat(i,j,k,l)=gmat(i,j,l,k)-gmat(i,j,k,l)&
!                               -gmat(j,i,l,k)+gmat(j,i,k,l)-0.5d0*adelt

                  !write(2323,2324)i,j,k,l,hmat(i,j,k,l)
                end do
                end do
              end do
              end do
            end if
            l0=l0+orb%total(l)
            lx=lx+orb%occ(l)
          end do
          k0=k0+orb%total(k)
          kx=kx+orb%occ(k)
          end do
!          dmat(i,j)=hmat(i,j,i,j)
          j0=j0+orb%total(j)
          jx=jx+orb%occ(j)
        end do
        i0=i0+orb%total(i)
        ix=ix+orb%occ(i)
        end do

2324    format(i2,i2,i2,i2,f12.8)
        deallocate(MT) 

      end subroutine Hessian2


! to Stefan : Eq.52
!             It is work for all type Hessian elements
!             When implementing, people may separate it into i-a,a-v,i-v type rotations part (e.g. in Molcas)

      subroutine hessian(iorb,amat,gmat,hmat,dmat)

        integer::iorb
        double precision::amat(iorb,iorb),dmat(iorb,iorb)
        double precision::gmat(iorb,iorb,iorb,iorb)
        double precision::hmat(iorb,iorb,iorb,iorb)
        double precision djl,djk,dil,dik,adelt
        double precision,allocatable::aplus(:,:)
        ! should be a_dag ...

        allocate(aplus(iorb,iorb))
 
        djl=0.0d0
        djk=0.0d0
        dil=0.0d0
        dik=0.0d0
        dmat=0.0d0
        hmat=0.0d0 
        aplus=0.0d0

        do i=1,iorb
          do j=1,iorb
            aplus(i,j)=amat(i,j)+amat(j,i)
            !write(*,*)i,j,aplus(i,j)
          end do
        end do
        !stop

         !do i=1,iorb
         !  do j=1,iorb
         !    do k=1,iorb
         !      do l=1,iorb
         !        write(700,*)i,j,k,l,gmat(i,j,k,l)
         !      end do
         !    end do 
         !  end do
         !end do
         !stop 
        !
        do i=1,iorb
          do j=1,iorb
            do k=1,iorb
              do l=1,iorb
                adelt=0.0d0
                call delt(j,l,djl)
                call delt(j,k,djk)
                call delt(i,l,dil)
                call delt(i,k,dik)

                adelt=aplus(i,k)*djl-aplus(i,l)*djk-aplus(j,k)*dil&
                     +aplus(j,l)*dik

                hmat(i,j,k,l)=gmat(i,j,l,k)-gmat(i,j,k,l)&
                             -gmat(j,i,l,k)+gmat(j,i,k,l)-0.5d0*adelt

               !write(2323,2324)i,j,k,l,hmat(i,j,k,l)

              end do
            end do
            dmat(i,j)=hmat(i,j,i,j)
          end do
        end do

2324    format(i2,i2,i2,i2,f12.8)
!        write(*,*)"preonescf",iorb
!        stop

        deallocate(aplus)

      end subroutine hessian

! -------------------------------------

       subroutine delt(i,j,d)

         integer::i,j
         double precision d

         if(i.eq.j)then
           d=1.0d0
         else
           d=0.0d0
         end if

       end subroutine delt

! ------------------------------------------

      subroutine hessian_test(iorb,amat,gmat,hmat,dmat)

        integer::iorb
        double precision::amat(iorb,iorb),dmat(iorb,iorb)
        double precision::am2(iorb,iorb),dm2(iorb,iorb)
        double precision::gmat(iorb,iorb,iorb,iorb)
        double precision::hmat(iorb,iorb,iorb,iorb)
        double precision::gm2(iorb,iorb,iorb,iorb)
        double precision djl,djk,dil,dik,adelt
        double precision,allocatable::aplus(:,:)
! should be a_dag ...

!       in the test
        integer od(iorb)

        od=0

        open(unit=117,file="orb.reorder")
          do i=1,iorb
            read(117,*)od(i)
          end do
        close(117)
        do i=1,iorb
          write(*,*)i,od(i)
        end do
        !stop

        allocate(aplus(iorb,iorb))

        djl=0.0d0
        djk=0.0d0
        dil=0.0d0
        dik=0.0d0
        dmat=0.0d0
        hmat=0.0d0
        aplus=0.0d0

        do i=1,iorb
          do j=1,iorb
            am2(od(i),od(j))=amat(i,j)
          end do
        end do 
        amat=am2

        do i=1,iorb
          do j=1,iorb
            do k=1,iorb
              do l=1,iorb
                gm2(od(i),od(j),od(k),od(l))=gmat(i,j,k,l)
              end do
            end do
          end do
        end do 
        gmat=gm2
 
        do i=1,iorb
          do j=1,iorb
            aplus(i,j)=amat(i,j)+amat(j,i)
            !write(*,*)i,j,aplus(i,j)
          end do
        end do
        !stop

!         do i=1,iorb
!           do j=1,iorb
!             do k=1,iorb
!               do l=1,iorb
!                 write(700,2326)i,j,k,l,dabs(gmat(i,j,k,l))
!               end do
!             end do 
!           end do
!         end do
!         stop 
!2326    format(i2,i2,i2,i2,f8.3)
                

        do i=1,iorb
          do j=1,iorb
            do k=1,iorb
              do l=1,iorb
                adelt=0.0d0
                call delt(j,l,djl)
                call delt(j,k,djk)
                call delt(i,l,dil)
                call delt(i,k,dik)

                adelt=aplus(i,k)*djl-aplus(i,l)*djk-aplus(j,k)*dil&
                     +aplus(j,l)*dik

                hmat(i,j,k,l)=gmat(i,j,l,k)-gmat(i,j,k,l)&
                             -gmat(j,i,l,k)+gmat(j,i,k,l)-0.5d0*adelt

!               write(2323,2324)i,j,k,l,hmat(i,j,k,l)

              end do
            end do
            dmat(i,j)=hmat(i,j,i,j)
          end do
        end do

2324    format(i2,i2,i2,i2,f12.8)
!        write(*,*)iorb
!        stop

        deallocate(aplus)

      end subroutine hessian_test













