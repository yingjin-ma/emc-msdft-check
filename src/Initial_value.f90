      Subroutine Initial_value()

! Generate the initial T matrix

        use global_control
        use matrix

        allocate(mat2%Ederi(norb,norb))
        allocate(mat2%deltR(norb,norb))
        allocate(mat2%deltT(norb,norb))

        mat2%Ederi=0.0d0
        mat2%deltR=0.0d0
        mat2%deltT=0.0d0

         do i=1,orb%nsub
           allocate(mat1(i)%R(orb%total(i),orb%total(i)))
           mat1(i)%R=0.0d0 
           allocate(mat1(i)%T(orb%total(i),orb%total(i)))
           mat1(i)%T=0.0d0
           allocate(mat1(i)%U(orb%total(i),orb%total(i)))
           mat1(i)%U=0.0d0
           !allocate(mat1(i)%B(orb%total(i),orb%total(i)))
           !mat1(i)%B=0.0d0
           allocate(mat1(i)%Hdiag(orb%total(i),orb%total(i)))
           mat1(i)%Hdiag=0.0d0
           allocate(mat1(i)%Ederi(orb%total(i),orb%total(i)))
           mat1(i)%Ederi=0.0d0
         end do

         ioffset=0
         do i=1,orb%nsub
           do i1=1,orb%total(i)
             do j1=1,orb%total(i)
!             mat2%T(i1+ioffset,j1+ioffset)=mat1(i)%T(i1,j1)
             end do
           end do
           ioffset=ioffset+orb%total(i) 
         end do 

         do i=1,norb
           do j=1,norb
             mat2%U(i,j)=mat2%T(i,j) 
           end do
           mat2%U(i,i)=mat2%U(i,i)+1.0d0
         end do 

!         call print_mat(norb,norb,mat1(1)%R)
!         write(*,*)"+++++++++++"

      end Subroutine Initial_value
