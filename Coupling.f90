      Subroutine coupling(icycle)
        
        use global_control
        use coupling_terms
        use matrix

        integer iloop 
        iloop=0

        call coupling_prepare(iloop)
! to Stefan : Eq.56
        call MPSci_gradient() 
! to Stefan : Eq.53
        call coupling_CIORB()
! to Stefan : Eq.54
        call coupling_CICI()

! Which order to be coupled using the extended Hessian
        if(ith_HESS.eq.0)then
          

        else if(ith_HESS.eq.1)then
          write(2,*)"Effective Hessian for orbital Hessian"
          call coupling_solver2(icycle) !In Coupled_solver.f90
         !call coupling_solver3(icycle) !In Coupled_solver.f90
        else if(ith_HESS.eq.2)then
          write(2,*)" extended Hessian (full matrix)"
          call coupling_solver(icycle) !In Coupled_solver.f90
        end if

        call coupling_finish()  

      end subroutine coupling

! This gradients part need to be modified for gradient case
! For orb-opt case, I will use use half-transformation way
!     instead of this. 
      Subroutine MPSci_gradient()
        
        use global_control
        use coupling_terms
        use matrix

! Used in the 2 particle grad part
        double precision tmp_KL, E_reformed, E2, T1,T2

        E2=0.0d0
        E_reformed=0.0d0

        do iC=1,nC
          Ci(iC)%grad=0.0d0

! This is actually too fussy!! 
!   RDM1 -> RDM1%sym
          j0=0; k0=0
          do i=1,orb%nsub
            !allocate(mat1(i)%D(orb%act(i),orb%act(i)))
            mat1(i)%D=0.0d0
            do j=1,orb%act(i)
              do k=1,orb%act(i)
                mat1(i)%D(j,k)=Ci(iC)%d(j0+j,k0+k)
!                mat1(i)%D(j,k)=mat2%D(j0+j,k0+k)
                !mat1(i)%D(j,k)=Ci(iC)%d(j0+j,k0+k)*Ci(iC)%dv
              end do
            end do
            j0=j0+orb%act(i)
            k0=k0+orb%act(i)
          end do

! MPSci gradient from one particle part
          ioffset=0
          do i=1,orb%nsub
            do j=1,orb%occ(i)
              koffset=0
              do k=1,orb%nsub
                if(i.eq.k)then
                  do l=1,orb%occ(k)
                    Ci(iC)%grad=Ci(iC)%grad+T(ioffset+j,koffset+l)*mat1(i)%D(j,l)
!                    write(*,*)j,l,T(ioffset+j,koffset+l),mat1(i)%D(j,l)
!                    if(mat1(i)%D(j,l).ne.0)then
!                      print *,"T",ioffset+j,koffset+l,T(ioffset+j,koffset+l)
!                    end if   
                  end do
                end if
                koffset=koffset+orb%total(k)
              end do
            end do
            ioffset=ioffset+orb%total(i)
          end do
!          print *,"Gradient (one) for ",iC,"-th vector is ",Ci(iC)%grad

! MPSci gradient from two particle part
          ioffset=0
          ioffset1=0
!         RDM_IJKL=0.0d0
          do i=1,orb%nsub; do i1=1,orb%occ(i)
            joffset=0
            joffset1=0
            do j=1,orb%nsub; do j1=1,orb%occ(j)
              tmp_KL=0.0d0
              koffset=0
              koffset1=0
              do k=1,orb%nsub; do k1=1,orb%occ(k)
                loffset=0
                loffset1=0
                do l=1,orb%nsub; do l1=1,orb%occ(l)
                  tmp_KL=tmp_KL+&
                    0.5d0*U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)&
              *Ci(iC)%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
!                  write(*,*)"U(ijkl)",U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1), &
!                            "P(iklj)",Ci(iC)%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
!                  if(Ci(iC)%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1).ne.0)then
!                    print *,"U",ioffset+i1,joffset+j1,koffset+k1,loffset+l1,U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)
!                  end if   
                end do
                loffset=loffset+orb%total(l)
                loffset1=loffset1+orb%occ(l)
                end do
              end do
              koffset=koffset+orb%total(k)
              koffset1=koffset1+orb%occ(k)
              end do
              Ci(iC)%grad=Ci(iC)%grad+tmp_KL
            end do
            joffset=joffset+orb%total(j)
            joffset1=joffset1+orb%occ(j)
            end do
          end do
          ioffset=ioffset+orb%total(i)
          ioffset1=ioffset1+orb%occ(i)
          end do

          E_reformed=E_reformed+Ci(iC)%grad*Ci(iC)%dv

          write(1,*)Ci(iC)%dv,"RDM_E", Ci(iC)%grad
!          write(*,*)" Hc",RDM_ENERGY
          write(1,*)Ci(iC)%dv,"  E*c",RDM_ENERGY*Ci(iC)%dv

          Ci(iC)%grad=Ci(iC)%grad-RDM_ENERGY*Ci(iC)%dv
          Ci(iC)%grad=Ci(iC)%grad*2.0d0

          write(1,*) "Gradient for ",iC,"-th vector is ",Ci(iC)%grad

        end do

        write(2,*) "Reformed Energy (all MPS_ci) is ",E_reformed
        write(2,*) "              The RDM Energy is ",RDM_ENERGY
        call flush(6)

        T2=0.0d0
        T1=0.0d0
        do iC=1,nC
          T2=Ci(iC)%grad*Ci(iC)%dv
          T1=T1+Ci(iC)%dv**2 
        end do
        E2=T2/T1
!        print *," The residual energy ",E2

!        do iC=2,nC
!          Ci(iC)%grad=Ci(iC)%grad-E2*Ci(iC)%dv
!          print *,"Gradient for ",iC,"-th vector (new) is ",Ci(iC)%grad
!        end do

!        stop
        
      End Subroutine MPSci_gradient   

! The orb-CI/orb-MPS coupling (i.e. orb-ci/orb-MPS Hessian part)
      Subroutine coupling_CIORB()

        use global_control
        use coupling_terms
        use matrix

        double precision,allocatable::T1(:,:),T2(:,:),T3(:,:)

!        character*72 oneRDMr 
!        character*72 twoRDMr 

        do iC=1,nC           
! Read in coefficients of Ci & DI, PI

!         The use of RDM derivative in the full form  

          ! The 1-order gradient for each element of Ci(davidson elems) in MPS
          if(.not.allocated(Ci(iC)%A))then
            allocate(Ci(iC)%A(norb,norb))
          end if
          if(.not.allocated(Ci(iC)%AL))then
            allocate(Ci(iC)%AL(norb,norb))
          end if
          if(.not.allocated(Ci(iC)%AR))then
            allocate(Ci(iC)%AR(norb,norb))
          end if
          Ci(iC)%A=0.0d0
          Ci(iC)%AL=0.0d0
          Ci(iC)%AR=0.0d0

! The hD part of Fock^i matrix 
!(very old code, should be updated later if necessary)
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%occ(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
! symmetried A
                Ci(iC)%A(n1+noffset,m1+moffset)=Ci(iC)%A(n1+noffset,m1+moffset)&
                +2.0d0*T(n1+noffset,j1+joffset)*Ci(iC)%D(m1+moffset1,j1+joffset1)
! AL
                Ci(iC)%AL(n1+noffset,m1+moffset)=Ci(iC)%AL(n1+noffset,m1+moffset)&
                +2.0d0*T(n1+noffset,j1+joffset)*Ci(iC)%DL(m1+moffset1,j1+joffset1)
! AR
                Ci(iC)%AR(n1+noffset,m1+moffset)=Ci(iC)%AR(n1+noffset,m1+moffset)&
                +2.0d0*T(n1+noffset,j1+joffset)*Ci(iC)%DR(m1+moffset1,j1+joffset1)
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%occ(j)
              end do
            end do
            moffset=moffset+orb%total(m)
            moffset1=moffset1+orb%occ(m)
            end do
          end do
          noffset=noffset+orb%total(n)
          noffset1=noffset1+orb%occ(n)
          end do

!         write(*,*)"The A operator only 1-partical"
!         call print_mat(norb,nact,Ci(iC)%A,6)

! The sum_kl[J^(kl)P^(lk)] part of Fock^i matrix 
! (very old code, should be updated later if necessary)         
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0
            do m=1,orb%nsub; do m1=1,orb%occ(m)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%occ(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%occ(k)
                  loffset=0
                  loffset1=0
                  do l=1,orb%nsub; do l1=1,orb%occ(l)
! Averaged A
                    Ci(iC)%A(n1+noffset,m1+moffset)=Ci(iC)%A(n1+noffset,m1+moffset)&
             +2.0d0*Ci(iC)%P(m1+moffset1,k1+koffset1,l1+loffset1,j1+joffset1)&
                     *U(n1+noffset,j1+joffset,k1+koffset,l1+loffset)
! AL 
                    Ci(iC)%AL(n1+noffset,m1+moffset)=Ci(iC)%AL(n1+noffset,m1+moffset)&
             +2.0d0*Ci(iC)%PL(m1+moffset1,k1+koffset1,l1+loffset1,j1+joffset1)&
                     *U(n1+noffset,j1+joffset,k1+koffset,l1+loffset)
! AR
                    Ci(iC)%AR(n1+noffset,m1+moffset)=Ci(iC)%AR(n1+noffset,m1+moffset)&
             +2.0d0*Ci(iC)%PR(m1+moffset1,k1+koffset1,l1+loffset1,j1+joffset1)&
                     *U(n1+noffset,j1+joffset,k1+koffset,l1+loffset)
                  end do
                  loffset=loffset+orb%total(l)
                  loffset1=loffset1+orb%occ(l)
                  end do
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%occ(j)
              end do
            end do
            moffset=moffset+orb%total(m)
            moffset1=moffset1+orb%occ(m)
            end do
          end do
          noffset=noffset+orb%total(n)
          noffset1=noffset1+orb%occ(n)
          end do           

!          deallocate(Ci(iC)%D)
!          deallocate(Ci(iC)%P)
!          Ci(iC)%A=Ci(iC)%A-Ci(iC)%dv*mat2%A
!          write(*,*)"The MPS_",ic,"'s A operator"
!          call print_mat(nact,nact,Ci(iC)%D,6)
!          call print_mat(norb,norb,Ci(iC)%A,6)

        end do
 
!        write(6,*)"The orbital's A matrix " 
!        call print_mat(norb,norb,mat2%A,6)
 
        write(1,*)"Finishing the contructions of the Fock-derivative" 
        call flush(1) 
!        stop

!  I) One way by constructing the full MPS-orb Hessian
        do iC=1,nC 
          do i=1,norb
            do j=1,norb
! This means add the abs(Dv) in constructing A^I
             ! Actually same as below
!              HCR(iC,i,j)=Ci(iC)%AL(i,j)-Ci(iC)%AL(j,i)&
!                         +Ci(iC)%AR(i,j)-Ci(iC)%AR(j,i)
!              HCR(iC,i,j)=1.0d0*HCR(iC,i,j)*0.5d0 
 
              KCR(iC,i,j)=Ci(iC)%A(i,j)-Ci(iC)%A(j,i) ! For property check

              HCR(iC,i,j)=Ci(iC)%A(i,j)-Ci(iC)%A(j,i)&
              -Ci(iC)%dv*mat2%A(i,j)+Ci(iC)%dv*mat2%A(j,i) ! <-- second line ()

!              HCR(iC,i,j)=Ci(iC)%A(i,j)-Ci(iC)%A(j,i)&
!              -Ci(iC)%dv*mat2%A(i,j)+Ci(iC)%dv*mat2%A(j,i)

!              HCR(iC,i,j)=  dabs(Ci(iC)%dv)*Ci(iC)%A(i,j)&
!                           -dabs(CI(iC)%dv)*Ci(iC)%A(j,i)&
!                                   -Ci(iC)%dv*mat2%A(i,j)&
!                                   +Ci(iC)%dv*mat2%A(j,i)
            end do
          end do
        end do
        HCR=1.0d0*HCR*2.0d0 
!        call print_TAT(nC,norb,norb,HCR,6)
!        stop


! =====================================================
       ! Check < H_cR | c > == 2 * Grad(orb)

        ndim1=norb
        allocate(T1(ndim1,ndim1))
        allocate(T2(ndim1,ndim1))
        allocate(T3(ndim1,ndim1))
        T1=0.0d0
        T2=0.0d0
        T3=0.0d0
 
        do i=1,ndim1
          do j=1,ndim1
            T2(i,j)=1.0d0*(mat2%A(i,j)-mat2%A(j,i))
            if(dabs(T2(i,j)).lt.1.0d-9)then
              T2(i,j)=0.0d0
            end if 
          end do
        end do

        do i=1,ndim1
          do j=1,ndim1
            do iC=1,nC
              T1(i,j)=T1(i,j)+KCR(iC,i,j)*Ci(iC)%dv
              T3(i,j)=T3(i,j)+HCR(iC,i,j)*Ci(iC)%dv
            end do
            if(dabs(T1(i,j)).lt.1.0d-9)then
              T1(i,j)=0.0d0
            end if 
            if(dabs(T3(i,j)).lt.1.0d-9)then
              T3(i,j)=0.0d0
            end if 
          end do
        end do    
     
!        write(*,*)"Check 1/2 < H_cR (K_cR) | c > == Grad(orb)"
!        call print_mat(ndim1,ndim1,T2,6) 
!        call print_mat(ndim1,ndim1,T1,6) 
!        write(*,*)"Check < H_cR | c > == 0, if convergenced c " 
!        call print_mat(ndim1,ndim1,T3,6) 
 
        deallocate(T1,T2)  
! ====================================================

! II) Other way is use the iterating only regarding to the H diagional & T*A^Ci 
!     (It can be added in the solver)

        write(1,*)"Finishing the coupling_MPS-ORB subroutine" 
        call flush(1)

      end Subroutine coupling_CIORB

! The CI/CI or MPS/MPS coupling 

      Subroutine coupling_CICI()

        use global_control
        use coupling_terms
        use matrix

        double precision Emcscf 
        double precision,allocatable::T1(:),Hami_all(:,:)
        double precision V1(nC)

! Here should be the read in the the full Hamiltonian
        ! read in code 
        allocate(Hami_all(nC_all,nC_all))
        Hami_all=0.0d0
        write(2,*)"nC_all in coupling_CICI",nC_all 
        open(unit=100,file="Hami.0.txt")
          do iC=1,nC_all
            read(100,*)(Hami_all(iC,j),j=1,nC_all) 
          end do
        close(100)
        do iC=1,nC
          do jC=1,nC
            Hami(iC,jC)=Hami_all(Ci(iC)%num,Ci(jC)%num)
          end do
        end do

! Set vary small value to 0
        do i=1,nC
          do j=1,nC
            if(dabs(Hami(i,j)).lt.dthrs)then
              Hami(i,j)=0.0d0  
            end if
          end do
        end do

! Also Energy (Assuming Erdm.eq.Edav ) 
        Emcscf=RDM_ENERGY

!        write(6,*)" Emcscf ", Emcscf 
!        write(6,*)" ====== Hamiltonian  start ====== " 
!        call print_mat(nc,nc,Hami,6)
!        write(6,*)" ====== Hamiltonian finish ====== " 
!        call flush(6) 

        do iC=1,nC            
          do jC=1,nC            
            HCC(iC,jC)=Hami(iC,jC)
          end do
          HCC(iC,iC)=Hami(iC,iC)-Emcscf ! for all
        end do
        HCC=HCC*2.0d0

! =====================================================
       ! Check < H_cc | c > ==  Grad(mps_ci)
        allocate(T1(nC))
        T1=0.0d0
        do iC=1,nC 
          do jC=1,nC
            T1(iC)=T1(iC)+HCC(iC,jC)*Ci(jC)%dv
          end do
        end do

!        write(*,*)"The H_cc * C (should always 0)"
!        call print_mat(1,nC,T1,6)
 
        deallocate(T1)
! =====================================================

!        write(6,*)" ====== MPS-MPS Hessian check ====== " 
!        call print_mat(nc,nc,Hami,6)
        call MXV(nC,nC,Hami,Ci(1:nC)%dv,V1) 
!        call print_mat(1,nC,V1,6)
        V1=0.0d0
        do iC=1,nC
          V1(iC)=Emcscf*Ci(iC)%dv
        end do 
!        call print_mat(1,nC,V1,6)
!        write(6,*)" ====== MPS-MPS Hessian check ====== " 
!        call flush(6)

        Hvec=Hami
        call EIGTQL2(nC,Hvec,Hval) 
! Unity the sign (to QCMaquis) if necessary
        if(Ci(1)%dv/Hvec(1,1).lt.0)then
          Hvec=-1.0d0*Hvec
        end if

! Set very small value to 0
        do i=1,nC
          if(dabs(Hval(i)).lt.dthrs)then
            Hval(i)=0.0d0
          end if
          do j=1,nC
            if(dabs(Hvec(i,j)).lt.dthrs)then
              Hvec(i,j)=0.0d0  
            end if
          end do
        end do

!        write(6,*)" ====== diagionalization ====== " 
!        call print_mat(nc,nc,Hvec,6)
!        write(6,*)" ====== diagionalization ====== " 
!        write(6,*)"diagionalized Eval(1)", Hval
!        call flush(6)
!        stop

      end Subroutine coupling_CICI


      Subroutine coupling_prepare(iloop)

        use global_control
        use coupling_terms
        use matrix
        
        integer::iloop  

! For the test case; check the completeness of RDM-deri       
        character*72 oneRDMr 
        character*72 twoRDMr 
! Temporary matrix for two-RDM-deri scratch
        double precision dv
! Whole 1' & 2' RDMs for checking correctness
        double precision oRDMchk(nact,nact)
        double precision tRDMchk(nact,nact,nact,nact)
! Temparory character input
        character ctmp

! Here should be read in the number of davidson elements, 
       ! Should from text file
       ! List all the 1'-RDM-deri, for counting the davidson elements 
        call system("ls oneRDM.0.deriL_* > info1.tmp")
        call system("ls twoRDM.0.deriL_* > info2.tmp")
        call system("ls oneRDM.0.deriR_* > info3.tmp")
        call system("ls twoRDM.0.deriR_* > info4.tmp")
        call system("wc -l info1.tmp > info_line.tmp")
        open(unit=100,file="info_line.tmp")
          read(100,*)nC
        close(100)
        write(2,*)"Number of Davidson vectors that need to be coupled :"
        write(2,*) nC

        ! Ci is a class for handling RDM derivative
        if(.not.allocated(Ci))then
          allocate(Ci(nC))     
        end if
        iC=0     

        if(iloop.eq.0)then ! only read in the coefficients  

          if(.true.)then
            write(2,*)"!!! > Readin cofficient from MPSCi.info"
            open(unit=100,file="MPSCi.0.info")
              read(100,*)ntype
              inum=0 
              do i=1,ntype
                read(100,*)j0,k0 
                do j=1,j0
                  do k=1,k0
                    inum=inum+1 
                    read(100,*)dv
                    if(dabs(dv).gt.1.0d-12)then
                      iC=iC+1
                      Ci(iC)%dv=dv
                      Ci(iC)%num=inum
                    end if 
                  end do
                end do
              end do
              nC_all=inum 
            close(100) 
          else 
            write(2,*)"!!!> Readin cofficient from MPSCi_value.info"
            open(unit=100,file="MPSCi_value.0.info")
              inum=0
              iC=0
              do 
                read(100,*,iostat=ierr)dv,i,j,k
                if(ierr.ne.0)exit 
                inum=inum+1
                if(dabs(dv).lt.dthrs)then
                else 
                  iC=iC+1
                  Ci(iC)%dv=dv
                  Ci(iC)%num=inum
                  !write(*,*)dv,inum
                end if 
              end do
              nC_all=inum 
            close(100)
          end if
          write(2,*)"Total Davidson vectors",nC_all
   
          ! read in the sign info. from overlap 
          open(unit=100,file="Loverlap.0.txt")
            do iC=1,nC
              read(100,*)ctmp,ctmp,ctmp,ctmp, dv
              read(100,*)
              Ci(iC)%dsL = sign(1.0d0,dv)
            end do 
          close(100)
          open(unit=100,file="Roverlap.0.txt")
            do iC=1,nC
              read(100,*)ctmp,ctmp,ctmp,ctmp, dv
              read(100,*)
              Ci(iC)%dsR = sign(1.0d0,dv)
            end do 
          close(100)

        else   
!          Ci(:)%dv=-1.0d0*Ci(:)%dv

        end if ! if now in micro-iters

!        write(6,*)" Read in sign of the Ci(:)%ds "
!        write(6,*)Ci(:)%ds 
!        stop

! For read in 1-RDM-deri (<i|o|mps> part) 
        open(unit=100,file="info1.tmp")
          do iC=1,nC 
            oneRDMr=""
! Averaged 1-RDM-deri
            if(.not.allocated(Ci(iC)%D))then
              allocate(Ci(iC)%D(nact,nact))
              Ci(iC)%D=0.0d0
            else
              Ci(iC)%D=0.0d0             
            end if
! Left side part
            if(.not.allocated(Ci(iC)%DL))then
              allocate(Ci(iC)%DL(nact,nact))           
              Ci(iC)%DL=0.0d0
            else
              Ci(iC)%DL=0.0d0
            end if
            read(100,*,iostat=ierr)oneRDMr 
            if(ierr.ne.0)exit 
            !write(6,*)"File of RDM-deri : ", oneRDMr
            call flush(6)
            open(unit=101,file=oneRDMr)
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  ! if NAN set to 0
                  if (isnan(dv).or.dabs(dv).lt.dthrs)then
                    dv = 0.0d0
                  end if 
                  Ci(iC)%DL(iorder(ij+1),iorder(ji+1))=dv
                end do
              end do              
            close(101)
          end do 
        close(100)
        write(2,*)"Finish importing 1-RDM-deri L"
        call flush(2)

! For read in 1-RDM-deri (<mps|o|i> part) 
        open(unit=100,file="info3.tmp")
          do iC=1,nC
            oneRDMr=""
! Right side part
            if(.not.allocated(Ci(iC)%DR))then
              allocate(Ci(iC)%DR(nact,nact))
              Ci(iC)%DR=0.0d0
            else
              Ci(iC)%DR=0.0d0
            end if
            read(100,*,iostat=ierr)oneRDMr
            if(ierr.ne.0)exit
            !write(6,*)"File of RDM-deri : ", oneRDMr
            call flush(2)
            open(unit=101,file=oneRDMr)
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  ! if NAN set to 0
                  if(isnan(dv).or.dabs(dv).lt.dthrs)then
                    dv = 0.0d0
                  end if
                  Ci(iC)%DR(iorder(ij+1),iorder(ji+1))=dv
                end do
              end do
            close(101)
          end do
        close(100)
        write(2,*)"Finish importing 1-RDM-deri R"
        call flush(2)

! For read in 2-RDM-deri (L)
        open(unit=100,file="info2.tmp")
          do iC=1,nC 
            twoRDMr=""
! Avergaed 
            if(.not.allocated(Ci(iC)%P))then
              allocate(Ci(iC)%P(nact,nact,nact,nact))
              Ci(iC)%P=0.0d0
            else
              Ci(iC)%P=0.0d0
            end if
! Left side part
            if(.not.allocated(Ci(iC)%PL))then            
              allocate(Ci(iC)%PL(nact,nact,nact,nact))
              Ci(iC)%PL=0.0d0
            else 
              Ci(iC)%PL=0.0d0
            end if
            read(100,*,iostat=ierr)twoRDMr 
            if(ierr.ne.0)exit
            ! read in the compressed form

            open(unit=101,file=twoRDMr)
              read(101,*)ijkl
              ijkl=0  ! Just a counter
              do
                read(101,*,iostat=ierr)ij,jk,kl,li,dv
                if(ierr.ne.0)exit
                ijkl=ijkl+1
!                      write(*,*)ij+1,jk+1,kl+1,li+1,rdm2P(ij+1,jk+1,kl+1,li+1)
                Ci(iC)%PL(iorder(ij+1),&
                          iorder(jk+1),&
                          iorder(kl+1),&
                          iorder(li+1))=dv
                Ci(iC)%PL(iorder(kl+1),&
                          iorder(li+1),&
                          iorder(ij+1),&
                          iorder(jk+1))=dv
                Ci(iC)%PL(iorder(jk+1),&
                          iorder(ij+1),&
                          iorder(li+1),&
                          iorder(kl+1))=dv
                Ci(iC)%PL(iorder(li+1),&
                          iorder(kl+1),&
                          iorder(jk+1),&
                          iorder(ij+1))=dv
              end do
            close(101)

            Ci(iC)%PL=Ci(iC)%PL*2.0d0
          end do 
        close(100)
        write(2,*)"Finish importing 2-RDM-deri L"
        call flush(2)

! For read in 2-RDM-deri (R)
        open(unit=100,file="info4.tmp")
          do iC=1,nC
            twoRDMr=""
            if(.not.allocated(Ci(iC)%PR))then 
              allocate(Ci(iC)%PR(nact,nact,nact,nact))
              Ci(iC)%PR=0.0d0
            else
              Ci(iC)%PR=0.0d0       
            end if
            read(100,*,iostat=ierr)twoRDMr
            if(ierr.ne.0)exit
            ! read in the compressed form

            open(unit=101,file=twoRDMr)
              read(101,*)ijkl
              ijkl=0  ! Just a counter
              do
                read(101,*,iostat=ierr)ij,jk,kl,li,dv
                if(ierr.ne.0)exit
                ijkl=ijkl+1
                Ci(iC)%PR(iorder(ij+1),&
                          iorder(jk+1),&
                          iorder(kl+1),&
                          iorder(li+1))=dv
                Ci(iC)%PR(iorder(kl+1),&
                          iorder(li+1),&
                          iorder(ij+1),&
                          iorder(jk+1))=dv
                Ci(iC)%PR(iorder(jk+1),&
                          iorder(ij+1),&
                          iorder(li+1),&
                          iorder(kl+1))=dv
                Ci(iC)%PR(iorder(li+1),&
                          iorder(kl+1),&
                          iorder(jk+1),&
                          iorder(ij+1))=dv
              end do
            close(101)

            Ci(iC)%PR=Ci(iC)%PR*2.0d0
          end do
        close(100)
        write(2,*)"Finish importing 2-RDM-deri R"
        call flush(2)

! Testing for the correctness
        oRDMchk=0.0d0        
        tRDMchk=0.0d0
        do iC=1,nC
          write(2,*)"Davidson C_",iC," : ", Ci(iC)%dv,& 
          dabs(Ci(iC)%dv)*Ci(iC)%dsL, dabs(Ci(iC)%dv)*Ci(iC)%dsR
!          Ci(iC)%dv=dabs(Ci(iC)%dv)*Ci(iC)%ds
!          Ci(iC)%D=Ci(iC)%D*Ci(iC)%ds*dabs(Ci(iC)%dv)
!          Ci(iC)%P=Ci(iC)%P*Ci(iC)%ds*dabs(Ci(iC)%dv)
!          oRDMchk=oRDMchk+Ci(iC)%D
!          tRDMchk=tRDMchk+Ci(iC)%P                    
 
          write(32,*)iC,"-iC",Ci(iC)%dv,&
                      "Ci(iC)%dv",Ci(iC)%dsL,Ci(iC)%dsR,"L and R"
          write(32,*)"L-deri"
          call print_mat(nact,nact,Ci(iC)%DL,32)
          write(32,*)"R-deri"
          call print_mat(nact,nact,Ci(iC)%DR,32)

       ! Fix the sign (only for the start iterations)
          !if(iloop.eq.0)then
            if(Ci(iC)%dv.eq.dabs(Ci(iC)%dv)*Ci(iC)%dsL)then
            else
              Ci(iC)%DL=-1.0d0*Ci(iC)%DL
              Ci(iC)%PL=-1.0d0*Ci(iC)%PL
            end if
            if(Ci(iC)%dv.eq.dabs(Ci(iC)%dv)*Ci(iC)%dsR)then
            else
              Ci(iC)%DR=-1.0d0*Ci(iC)%DR
              Ci(iC)%PR=-1.0d0*Ci(iC)%PR
            end if
            write(1,*)"The signs are forced changed before micro-it 0 ?"
          !end if

          if(iloop.gt.0)then 
            Ci(iC)%DL=-Ci(iC)%DL
            Ci(iC)%PL=-Ci(iC)%PL
            Ci(iC)%DR=-Ci(iC)%DR
            Ci(iC)%PR=-Ci(iC)%PR
          end if

! Sum up L and R, then average
          Ci(iC)%D=Ci(iC)%DL+Ci(iC)%DR
          Ci(iC)%P=Ci(iC)%PL+Ci(iC)%PR

          Ci(iC)%D=Ci(iC)%D/2.0d0
          Ci(iC)%P=Ci(iC)%P/2.0d0

!          Ci(iC)%D=Ci(iC)%D-Ci(iC)%Dv*mat2%D
!          Ci(iC)%P=Ci(iC)%P-Ci(iC)%Dv*mat2%P
!          oRDMchk=oRDMchk+Ci(iC)%D*dabs(Ci(iC)%dv)*Ci(iC)%ds
!          tRDMchk=tRDMchk+Ci(iC)%P*dabs(Ci(iC)%dv)*Ci(iC)%ds
          oRDMchk=oRDMchk+Ci(iC)%D*Ci(iC)%dv
          tRDMchk=tRDMchk+Ci(iC)%P*Ci(iC)%dv
!          deallocate(Ci(iC)%DL)
!          deallocate(Ci(iC)%PL)
        end do        
        write(2,*)"Notice that the corrected Davidson vector c is used."
 
        open(unit=100,file="oneRDM.reformed")
          do i=1,nact
            do j=1,nact
              write(100,*)"oRDMchk",i,j,oRDMchk(i,j)
            end do
          end do  
        close(100) 
        open(unit=100,file="twoRDM.reformed")
          do i=1,nact
            do j=1,nact
              do k=1,nact
                do l=1,nact
                  write(100,*)"tRDMchk",i,j,k,l,tRDMchk(i,j,k,l)/2.0d0
                end do 
              end do
            end do
          end do  
        close(100)

        write(2,*)"DMRG Calculated RDM : " 
        call Print_mat(nact,nact,mat2%d,2)
        write(2,*)"      resambled RDM : " 
        call Print_mat(nact,nact,oRDMchk,2)

! nC is the dim (number of) of the Davidson vecotrs
!        nC=16 ! Use 10 as a example

! Allocate the orb-CI and CI-CI Hessian
        if(.not.allocated(HCR))then
          allocate(HCR(nC,norb,norb)); HCR=0.0d0
        else
          HCR=0.0d0
        end if
        if(.not.allocated(KCR))then
          allocate(KCR(nC,norb,norb)); KCR=0.0d0
        else
          KCR=0.0d0
        end if  
        if(.not.allocated(HCC))then
          allocate(HCC(nC,nC)); HCC=0.0d0
        else
          HCC=0.0d0
        end if
! Allocate for the full Hamiltonian
        if(.not.allocated(Hami))then
          allocate(Hami(nC,nC)) ; Hami=0.0d0
        else
          Hami=0.0d0
        end if 
        if(.not.allocated(Hvec))then
          allocate(Hvec(nC,nC)) ; Hvec=0.0d0
        else
          Hvec=0.0d0
        end if
        if(.not.allocated(Hval))then
          allocate(Hval(nC))    ; Hval=0.0d0
        else
          Hval=0.0d0
        end if

      end Subroutine coupling_prepare

      Subroutine coupling_finish()  

        use global_control
        use coupling_terms
        use matrix

        deallocate(HCC)        
        deallocate(HCR)
        deallocate(KCR)
        deallocate(Ci)
        deallocate(Hami,Hvec,Hval)

      End Subroutine coupling_finish
