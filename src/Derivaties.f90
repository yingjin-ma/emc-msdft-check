!REF H.-J. Werner and W. Meyer, J. Chem. Phys. 73, 2342 (1980).

      Subroutine Derivaties(icycle)

        use global_control
        use matrix

        double precision UB(norb,norb)
        double precision BU(norb,norb)
        !write(*,*)mat2%B 
        !write(*,*)mat2%T 
        
        mat2%B=0.0d0
!  Generate the B-matrix
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%total(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%total(k)
!      induces for rdm
                  nr=n1+noffset1
                  mr=m1+moffset1
                  jr=j1+joffset1
                  kr=k1+koffset1
!      induce  for integrals
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset
!              B   part
                  mat2%B(ni,mi)=mat2%B(ni,mi)&
                             +mat2%G(ni,mi,ji,ki)*mat2%T(ki,ji)  ! i->r Mar.16 2015 (not yet)
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
             end do
             joffset=joffset+orb%total(j)
             joffset1=joffset1+orb%occ(j)
             end do
             !write(*,*)n1+noffset,m1+moffset,mat2%B(n1+noffset,m1+moffset)
           end do
           moffset=moffset+orb%total(m)
           moffset1=moffset1+orb%occ(m)
           end do
         end do
         noffset=noffset+orb%total(n)
         noffset1=noffset1+orb%occ(n)
         end do

         mat2%B=mat2%B+mat2%A

!         do i=1,norb
!           do j=1,norb
!             write(3333,*)i,j,mat2%B(i,j)
!           end do
!         end do
!         call print_mat(norb,norb,mat2%B)
!         stop

         !write(*,*)"T U ==================___++++_+_)+"

!          do i=1,norb
!            do j=1,norb
!              write(704,*)i,j,mat2%T(i,j),mat2%B(i,j)
!            end do
!          end do

         UB=0.0d0
         BU=0.0d0 
         mat2%Ederi=0.0d0

         call MdXM(norb,mat2%U,mat2%B,UB)
         call MdXM(norb,mat2%B,mat2%U,BU)

         !write(*,*)"UBBU1"
         !stop

         mat2%Ederi=UB-BU

! mat1%B is decleared in PreOneSCF
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
             !write(*,*)i1+ioffset,j1+joffset,mat2%A(i1+ioffset,j1+joffset),mat2%B(i1+ioffset,j1+joffset) 
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
               mat1(i)%Hdiag(i1,j1)=mat2%Hdiag(i1+ioffset,j1+ioffset)
               mat1(i)%Ederi(i1,j1)=mat2%Ederi(i1+ioffset,j1+ioffset)
!               write(707,*)i,i1,j1,mat1(i)%Ederi(i1,j1)
             end do
           end do
           ioffset=ioffset+orb%total(i)
!           call print_mat(orb%total(i),orb%total(i),mat1(i)%Hdiag)
         end do

      end subroutine Derivaties

! ======================================================================

      Subroutine derivatives_B_update(Bini)

        use global_control
        use matrix
        
        double precision::Bini(norb,norb)

! B = A 
        mat2%B=Bini
! B = A + G*T

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
                mat2%B(ni,mi)=mat2%B(ni,mi)&
                           +mat2%G(ni,mi,ji,ki)*mat2%T(ki,ji)  
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

      End Subroutine derivatives_B_update

!=======================================================================
! REF H.-J. Werner, Adv. Chem. Phys. LXIX, 1 (1987).
! REF H.-J. Werner and P. J. Knowles, J. Chem. Phys. 82, 5053 (1985).

      Subroutine Derivaties2(icycle,iloop)

        use global_control
        use matrix
 
        integer::icycle,iloop
        double precision UB(norb,norb)
        double precision BU(norb,norb)
        double precision UdR(norb,norb)
        double precision Bd(norb,norb)
        double precision AA(norb,norb)
        double precision AdR(norb,norb)
        double precision dRA(norb,norb)

        Bd=0.0d0
        UdR=0.0d0

        !write(*,*)mat2%U
        !write(*,*)"============"
        !write(*,*)mat2%deltR
        !do i=1,norb
        !  do j=1,norb
        !    write(*,*)"deltR",i,j,mat2%deltR(i,j)
        !  end do
        !end do
        !stop
        !write(*,*)"============"

        call MXM(norb,mat2%U,mat2%deltR,UdR)

        !do i=1,norb
        !  do j=1,norb
        !    write(*,*)i,j,UdR(i,j)
        !  end do
        !end do
        !stop
        !write(*,*)mat2%B 
        !write(*,*)mat2%T 
        
        mat2%B=0.0d0
!  Generate the B-matrix
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%total(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%total(k)
!      induces for rdm
                  nr=n1+noffset1
                  mr=m1+moffset1
                  jr=j1+joffset1
                  kr=k1+koffset1
!      induce  for integrals
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset
!              B   part
                  mat2%B(ni,mi)=mat2%B(ni,mi)&
                             +mat2%G(ni,mi,ji,ki)*mat2%T(ki,ji)  ! i -> r Mar. 16, 2015 (not yet)
                  Bd(ni,mi)=Bd(ni,mi)+mat2%G(ni,mi,ji,ki)*mat2%T(ki,ji)&
                                     +mat2%G(ni,mi,ji,ki)*UdR(ki,ji) 
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
             end do
             joffset=joffset+orb%total(j)
             joffset1=joffset1+orb%occ(j)
             end do
             !write(*,*)n1+noffset,m1+moffset,mat2%B(n1+noffset,m1+moffset)
           end do
           moffset=moffset+orb%total(m)
           moffset1=moffset1+orb%occ(m)
           end do
         end do
         noffset=noffset+orb%total(n)
         noffset1=noffset1+orb%occ(n)
         end do

         mat2%B=mat2%B+mat2%A
         Bd=Bd+mat2%A 
! The ~A=UdB
         call MdXM(norb,mat2%U,mat2%B,mat2%A)
! The ~B=B+G*UdR
         mat2%B=Bd 

         if(icycle.eq.0.or.mod(iloop,20).eq.0)then
           call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
           do i=1,norb
             do j=1,norb
               write(700,*)"Hdiag",i,j,mat2%Hdiag(i,j)   
             end do
           end do
         end if
          
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
               mat1(i)%A(i1,j1)=mat2%A(i1+ioffset,j1+ioffset)
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

         do i=1,norb
           do j=1,norb
             !write(707,*)i,j,mat2%B(i,j)
           end do
         end do
        ! stop

         !write(*,*)"T U ==================___++++_+_)+"

!          do i=1,norb
!            do j=1,norb
!              write(704,*)i,j,mat2%T(i,j),mat2%B(i,j)
!            end do
!          end do

         UB=0.0d0
         BU=0.0d0
         AA=0.0d0
         AdR=0.0d0
         dRA=0.0d0  
         mat2%Ederi=0.0d0

         call MdXM(norb,mat2%U,mat2%B,UB)
         call MdXM(norb,mat2%B,mat2%U,BU)

         do i=1,norb
           do j=1,norb
             AA(i,j)=mat2%A(i,j)+mat2%A(j,i)
           end do
         end do

         call MXM(norb,AA,mat2%deltR,AdR)
         call MXM(norb,mat2%deltR,AA,dRA)
         !write(*,*)"UBBU1"
        

         mat2%Ederi=UB-BU-0.5d0*(AdR+dRA)

!         do i=1,norb
!           do j=1,norb
!            write(712,*)i,j,mat2%Ederi(i,j)
!           end do
!         end do
!         stop

! mat1%B is decleared in PreOneSCF
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
             !write(*,*)i1+ioffset,j1+joffset,mat2%A(i1+ioffset,j1+joffset),mat2%B(i1+ioffset,j1+joffset) 
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
               mat1(i)%Hdiag(i1,j1)=mat2%Hdiag(i1+ioffset,j1+ioffset)
               mat1(i)%Ederi(i1,j1)=mat2%Ederi(i1+ioffset,j1+ioffset)
               write(707,*)i,i1,j1,mat1(i)%Ederi(i1,j1)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

      end subroutine Derivaties2

!=======================================================================

      Subroutine Derivaties3(icycle,iloop)

        use global_control
        use matrix
 
        integer::icycle,iloop
        double precision BB(norb,norb)
        double precision AA(norb,norb)
        double precision AdR(norb,norb)
        double precision dRA(norb,norb)

        mat2%B=0.0d0
!  Generate the B-matrix
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%total(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%total(k)
!      induces for rdm
                  nr=n1+noffset1
                  mr=m1+moffset1
                  jr=j1+joffset1
                  kr=k1+koffset1
!      induce  for integrals
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset
!              B   part
                  mat2%B(ni,mi)=mat2%B(ni,mi)&
                             +mat2%G(ni,mi,ji,ki)*mat2%R(ki,ji) ! Mar. 16. 2015 (not yet)
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
             end do
             joffset=joffset+orb%total(j)
             joffset1=joffset1+orb%occ(j)
             end do
             !write(*,*)n1+noffset,m1+moffset,mat2%B(n1+noffset,m1+moffset)
           end do
           moffset=moffset+orb%total(m)
           moffset1=moffset1+orb%occ(m)
           end do
         end do
         noffset=noffset+orb%total(n)
         noffset1=noffset1+orb%occ(n)
         end do

! The Newton-Raphson B'
         mat2%B=mat2%B+mat2%A

         if(icycle.eq.0.or.mod(iloop,20).eq.0)then
           call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
           do i=1,norb
             do j=1,norb
               write(700,*)"Hdiag",i,j,mat2%Hdiag(i,j)   
             end do
           end do
         end if
          
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
               mat1(i)%A(i1,j1)=mat2%A(i1+ioffset,j1+ioffset)
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

         !write(*,*)"T U ==================___++++_+_)+"

          do i=1,norb
            do j=1,norb
              write(704,*)i,j,mat2%T(i,j),mat2%B(i,j)
            end do
          end do

         BB=0.0d0
         mat2%Ederi=0.0d0

         do i=1,norb
           do j=1,norb
             BB(i,j)=mat2%B(i,j)-mat2%B(j,i)
             AA(i,j)=mat2%A(i,j)+mat2%A(j,i)
           end do
         end do

         call MXM(norb,AA,mat2%R,AdR)
         call MXM(norb,mat2%R,AA,dRA)

         mat2%Ederi=BB-0.5d0*(AdR+dRA)

!         do i=1,norb
!           do j=1,norb
!            write(712,*)i,j,mat2%Ederi(i,j)
!           end do
!         end do
!         stop

! mat1%B is decleared in PreOneSCF
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
             !write(*,*)i1+ioffset,j1+joffset,mat2%A(i1+ioffset,j1+joffset),mat2%B(i1+ioffset,j1+joffset) 
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
               mat1(i)%Hdiag(i1,j1)=mat2%Hdiag(i1+ioffset,j1+ioffset)
               mat1(i)%Ederi(i1,j1)=mat2%Ederi(i1+ioffset,j1+ioffset)
               write(707,*)i,i1,j1,mat1(i)%Ederi(i1,j1)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

      end subroutine Derivaties3

!=======================================================================

      Subroutine Derivaties4(icycle,iloop)

        use global_control
        use matrix
 
        integer::icycle,iloop
        double precision Bp(norb,norb)
        double precision Bm(norb,norb)
        double precision BR(norb,norb)
        double precision RB(norb,norb)
        double precision R2B(norb,norb)
        double precision RBR(norb,norb)
        double precision BR2(norb,norb)
        double precision T1(norb,norb)

        mat2%B=0.0d0
          !  Generate the B-matrix
          noffset=0
          noffset1=0
        do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%total(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%total(k)
            !      induces for rdm
                      nr=n1+noffset1
                      mr=m1+moffset1
                      jr=j1+joffset1
                      kr=k1+koffset1
            !      induce  for integrals
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset
              !              B   part
                  mat2%B(ni,mi)=mat2%B(ni,mi)&
                             +mat2%G(ni,mi,ji,ki)*mat2%R(ki,ji) ! Mar. 16. 2015 (not yet)
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
             end do
             joffset=joffset+orb%total(j)
             joffset1=joffset1+orb%occ(j)
             end do
             !write(*,*)n1+noffset,m1+moffset,mat2%B(n1+noffset,m1+moffset)
           end do
           moffset=moffset+orb%total(m)
           moffset1=moffset1+orb%occ(m)
           end do
         end do
         noffset=noffset+orb%total(n)
         noffset1=noffset1+orb%occ(n)
        end do

         mat2%B=mat2%B+mat2%A

         if(icycle.eq.0.or.mod(iloop,20).eq.0)then
           call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
           do i=1,norb
             do j=1,norb
               write(700,*)"Hdiag",i,j,mat2%Hdiag(i,j)   
             end do
           end do
           !write(*,*)"iloops"
         end if
          
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
               mat1(i)%A(i1,j1)=mat2%A(i1+ioffset,j1+ioffset)
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

         !write(*,*)"T U ==================___++++_+_)+"

          do i=1,norb
            do j=1,norb
              write(704,*)i,j,mat2%T(i,j),mat2%B(i,j)
            end do
          end do

         Bp=0.0d0
         Bm=0.0d0
         BR=0.0d0
         RB=0.0d0
         R2B=0.0d0
         BR2=0.0d0
         RBR=0.0d0
         T1=0.0d0
         mat2%Ederi=0.0d0

         do i=1,norb
           do j=1,norb
             Bp(i,j)=mat2%B(i,j)+mat2%B(j,i)
             Bm(i,j)=mat2%B(i,j)-mat2%B(j,i)
           end do
         end do

         call MXM(norb,Bp,mat2%R,BR)
         call MXM(norb,mat2%R,Bp,RB)

         call MXM(norb,Bm,mat2%R,T1)
         call MXM(norb,T1,mat2%R,BR2)

         call MXM(norb,mat2%R,Bm,T1)
         call MXM(norb,T1,mat2%R,RBR)
         
         call MXM(norb,mat2%R,mat2%R,T1)
         call MXM(norb,T1,Bm,R2B)

         mat2%Ederi=Bm-0.5d0*(BR+RB)+1.0d0/6.0d0*(BR2+RBR+R2B)

            !         do i=1,norb
            !           do j=1,norb
            !            write(712,*)i,j,mat2%Ederi(i,j)
            !           end do
            !         end do
            !         stop

            ! mat1%B is decleared in PreOneSCF
         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
             !write(*,*)i1+ioffset,j1+joffset,mat2%A(i1+ioffset,j1+joffset),mat2%B(i1+ioffset,j1+joffset) 
               mat1(i)%B(i1,j1)=mat2%B(i1+ioffset,j1+ioffset)
               mat1(i)%Hdiag(i1,j1)=mat2%Hdiag(i1+ioffset,j1+ioffset)
               mat1(i)%Ederi(i1,j1)=mat2%Ederi(i1+ioffset,j1+ioffset)
               write(707,*)i,i1,j1,mat1(i)%Ederi(i1,j1)
             end do
           end do
           ioffset=ioffset+orb%total(i)
         end do

      end subroutine Derivaties4

! Maodou, best wishes
