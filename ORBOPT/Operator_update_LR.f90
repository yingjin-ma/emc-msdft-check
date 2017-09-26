      Subroutine Operator_update_LR(iloop)

        use global_control
!        use coupling_terms
        use matrix

        integer::iloop

! For the 2-index transformed integrals
        type::transform
          double precision,allocatable::U(:,:)
          double precision,allocatable::T(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision U_updated(norb,norb)
        double precision U_tmp(norb,norb)
        double precision T_tmp(norb,norb)

        double precision dv,dtmp,RDM_energy_updated

        double precision Tact(nact,nact)
        double precision Uact(nact,nact,nact,nact) 

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)        

        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

!        write(21,*)"  =====  In the loop 0000  ", iloop, " ===== " 
!        call print_mat(nact,nact,mat2%D,21)
!        write(22,*)"  =====  In the loop 0000  ", iloop, " ===== "
!        call print_gat(nact,nact,nact,nact,mat2%P,22)

        do iroot=0,dmrg_nstates-1

          i=iroot+1
          write(21,*)"  =====  In the loop 0000  ", iloop, " ===== " 
          call print_mat(nact,nact,MPS(i)%D,21)
          write(22,*)"  =====  In the loop 0000  ", iloop, " ===== "
          call print_gat(nact,nact,nact,nact,MPS(i)%P,22)

          call RDM_update_MPSci_LR(MPS(i)%nC,        &
                                   MPS(i)%CI(:)%dv,  & 
                                   MPS(i)%CI(:)%num, & 
                                   MPS(i)%nC_all,    &
                                   nact,             &
                                   MPS(i)%D,         &
                                   MPS(i)%P,         & 
                                   iroot )

          write(21,*)"  =====  In the loop 1111 ",iloop, " ===== " 
          call print_mat(nact,nact,MPS(i)%D,21)
          write(22,*)"  =====  In the loop 1111 ",iloop, " ===== "
          call print_gat(nact,nact,nact,nact,MPS(i)%P,22)

          write(1,*)" RDMs are updated for ",iroot,"-th root "
          call flush(1)     

        end do

!        call RDM_update_MPSci_LR&
!             (nC,CI(1:nC)%dv,CI(1:nC)%num,nC_all,nact,mat2%D,mat2%P)

!        stop 
     
        call MXM(norb,UMAT_FINAL,mat2%U,U_updated)
        UMAT_FINAL=U_updated

        write(1,*)"In micro-iterations",iloop 
        call flush(1)
!        call print_mat(norb,norb,UMAT_FINAL,1)
!        call print_mat(norb,norb,U_updated,1)

! For quick test, I use 2-index transformation
! It can be easily back to 1-index transform, later together with the other 1-index and act only on 2-dim matrix (i.e. A and B)
        if(.true.)then
          call transform_2index_ALL&
               (orb%nsub,orb%total,norb,UMAT_FINAL,&
                T_origin,U_origin,T,U,orb%grouptable) 
        else
          ! to be added later
        end if
        
        mat2%U=0.0d0
        mat2%T=0.0d0
        do i=1,norb
          mat2%U(i,i)=1.0d0
        end do
        mat2%deltR=0.0d0

        write(1,*)"after 2/1-index transform" 
        call flush(1)

! 2index assumes it is the multipled form 
! retreat to 1-index form when combining Molcas later
        Tact=0.0d0; Uact=0.0d0
        call transform_2index&
            (orb%nsub,orb%act,orb%total,nact,norb,&
             MAT2%U,T,U,Tact,Uact,orb%grouptable)

!        goto 1111

! Distribute the transform matrix base on group (not used)
        allocate(UM(orb%nsub))
        i0=0
        do i=1,orb%nsub
          itot=orb%total(i)
          allocate(UM(i)%U(itot,itot))
          allocate(UM(i)%T(itot,itot))
          UM(i)%U  = 0.0d0
          UM(i)%T  = 0.0d0
          UM(i)%U  = mat2%U(i0+1:i0+itot,i0+1:i0+itot)
          UM(i)%T  = UM(i)%U
!         T=U-1 in sub-irreps
          do j=1,itot
            UM(i)%T(j,j)=UM(i)%U(j,j)-1.0d0
          end do
!          call print_mat(tot(i),tot(i),UM(i)%U)
          i0=i0+orb%total(i)
        end do

        if(.NOT.ALLOCATED(mat2%A))then
          allocate(mat2%A(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%B))then
          allocate(mat2%B(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%G))then
        ! still use the old form
          allocate(mat2%G(norb,norb,norb,norb)) 
        end if
!        if(.NOT.ALLOCATED(mat2%H))then
!          allocate(mat2%H(norb,norb,norb,norb)) !  
!        end if

! calculate the E(2) during the micro-iterations      
        RDM_energy_updated=0.0d0
     ! 1-e part
        allocate(TM1(nact,nact));TM1=0.0d0
        allocate(TM2(nact,nact));TM2=0.0d0
        ia=0
        do i=1,orb%nsub
          ij=orb%occ(i)
          TM1(ia+1:ia+ij,ia+1:ia+ij)=Tact(ia+1:ia+ij,ia+1:ia+ij)
          ia=ia+orb%occ(i)
        end do
        call MXM(nact,TM1,mat2%D,TM2)
        call trace(nact,TM2,dv)

        RDM_energy_updated=dv
!        write(6,*)"energy 1-e",RDM_energy_updated
        deallocate(TM1)
        deallocate(TM2)

     ! 2-e part
        k0=0
        ka=0
        do k=1,orb%nsub; do kk=1,orb%occ(k)
          l0=0
          la=0
          do l=1,orb%nsub; do ll=1,orb%occ(l)

          ! \sum_trace(J^{kl}*P^{lk})
            allocate(TM1(nact,nact));TM1=0.0d0
            allocate(TM2(nact,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(nact,nact));TM4=0.0d0

            TM1=  Uact(:,:,kk+ka,ll+la)  ! bk
            TM3=mat2%P(:,kk+ka,ll+la,:)  ! bk

!           call print_mat(nact,nact,TM1,6)

            i0=0
            ia=0
            do i=1,orb%nsub
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jj=orb%occ(j)
               TM2(ia+1:ia+ij,ja+1:ja+jj)=TM1(ia+1:ia+ij,ja+1:ja+jj)

                ja=ja+orb%occ(j)
              end do
              ia=ia+orb%occ(i)
            end do
!                write(*,*)" ===  TM1/TM2  === "
!                call print_mat(nact,nact,TM2,6)

            call MXM(nact,TM2,TM3,TM4)
            call trace(nact,TM4,dv)

            RDM_energy_updated=RDM_energy_updated+dv*0.5d0

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4)

          end do
          la=la+orb%occ(l)
          end do
        end do
        ka=ka+orb%occ(k)
        end do

!       write(*,*)"updated RDM_energy(1+2)",RDM_ENERGY_updated
        write(1,*)iloop,"-moit updated RDM_energy ",&
                  RDM_ENERGY_updated+fronre0

        E0_micro=RDM_ENERGY_updated

!       Construct the A matrix base on the transformed 1-integrals
        call fock_gen(orb%nsub,orb%act,orb%total,&
                  nact,norb,T,U,mat2%D,mat2%P,mat2%A)

        mat2%B=mat2%A

!        write(6,*)" mat2%A, updated in micro-iterations "
!        call print_mat(norb,norb,mat2%A,6)

        deallocate(TM1)
        deallocate(TM5)

        deallocate(UM)

!1111    continue
!  ==============
!     G="hD"+
!  ==============
         allocate(GM1(norb,norb,norb,norb)) ! the k,ij,l | The i,j can be further reduced
         GM1=0.0d0
         ioffset=0
         ioffset1=0
         do i=1,orb%nsub; do ii=1,orb%occ(i)
           i0=ioffset
           i1=orb%total(i)
           i2=orb%occ(i)
           ix=ioffset1

           joffset=0
           joffset1=0
           do j=1,orb%nsub; do jj=1,orb%occ(j)
             j0=joffset
             j1=orb%total(j)
             j2=orb%occ(j)
             jx=joffset1

             GM1(:,i0+ii,j0+jj,:)=T*mat2%D(ix+ii,jx+jj)*2.0d0
           end do
           joffset=joffset+orb%total(j)
           joffset1=joffset1+orb%occ(j)
           end do
         end do
         ioffset=ioffset+orb%total(i)
         ioffset1=ioffset1+orb%occ(i)
         end do
         mat2%G=GM1
!         write(*,*)"mat2%G || GM1" 
!         call print_GAT(norb,norb,norb,norb,GM1,8) 
         deallocate(GM1)

!  =============================================
!   G=hD+"sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"
!  =============================================
         allocate(GM1(norb,norb,norb,norb))
         GM1=0.0d0

         ioffset=0
         ioffset1=0
         do i=1,orb%nsub; do ii=1,orb%occ(i)
           i0=ioffset
           i1=orb%total(i)
           i2=orb%occ(i)
           ix=ioffset1

           joffset=0
           joffset1=0
           do j=1,orb%nsub; do jj=1,orb%occ(j)
             j0=joffset
             j1=orb%total(j)
             j2=orb%occ(j)
             jx=joffset1

             allocate(TM1(norb,norb));TM1=0.0d0

             koffset=0
             koffset1=0
             do k=1,orb%nsub; do kk=1,orb%occ(k)
               k0=koffset
               k1=orb%total(k)
               k2=orb%occ(k)
               kx=koffset1

               loffset=0
               loffset1=0
               do l=1,orb%nsub
! this is for JP
                 if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
                   do ll=1,orb%occ(l)
                     l0=loffset
                     l1=orb%total(l)
                     l2=orb%occ(l)
                     lx=loffset1
                     dtmp=mat2%P(ix+ii,kx+kk,lx+ll,jx+jj)
!                     if(dtmp.ne.0)then
!                       write(9,*)"i,j,k,l",i,j,k,l,"P",dtmp
!                       write(9,*)"i,j,k,l",i,j,k,l,"U",U(1,1,k0+kk,l0+ll)
!                     end if
                     TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp*2.0d0
                   end do
                 end if
! this is for KQ
                 if(orb%grouptable(i,k).eq.orb%grouptable(j,l))then
                   do ll=1,orb%occ(l)
                     l0=loffset
                     l1=orb%total(l)
                     l2=orb%occ(l)
                     lx=loffset1
                     dtmp=mat2%P(ix+ii,jx+jj,lx+ll,kx+kk)
                     TM1=TM1+U(:,k0+kk,l0+ll,:)*dtmp*4.0d0
                   end do
                 end if
                 loffset=loffset+orb%total(l)
                 loffset1=loffset1+orb%occ(l)
               end do
             end do
             koffset=koffset+orb%total(k)
             koffset1=koffset1+orb%occ(k)
             end do
             GM1(:,i0+ii,j0+jj,:)=TM1
             deallocate(TM1)

           end do
           joffset=joffset+orb%total(j)
           joffset1=joffset1+orb%occ(j)
           end do
         end do
         ioffset=ioffset+orb%total(i)
         ioffset1=ioffset1+orb%occ(i)
         end do
         mat2%G=mat2%G+GM1
         deallocate(GM1)
         write(2,*)"G operator is also updated in micro-iters"
!        write(1,*)"Before finishing Operator_update_CPNR"
!        call flush(1)

      End Subroutine Operator_update_LR


      Subroutine RDM_update_MPSci_LR&
                 (ndimC,Vmpsci,Vindex,nCA,nact,D,P,iroot)

         integer::ndimC,nCA,nact,iroot 
         integer::Vindex(ndimC)
         double precision::Vmpsci(ndimC)    
         double precision VmpsciA(nCA) 
         character*192 string1,oneRDMfile,twoRDMfile 
 
         character    ctmp1
         character*2  ctmp2

         double precision::D(nact,nact)
         double precision::P(nact,nact,nact,nact)


           ctmp2=""
           if(iroot.ge.10)then
             write(ctmp2,"(I2)")iroot
           else
             write(ctmp1,"(I1)")iroot
             ctmp2(1:1)=ctmp1
           end if          

           string1="updated_davidson.vec"
           
           VmpsciA=0.0d0
           do i=1,ndimC
             VmpsciA(Vindex(i))=Vmpsci(i)
           end do

           open(unit=100,file=trim(string1))
           do i=1,nCA
             write(100,*)1.0d0*VmpsciA(i)
           end do
           close(100)

           call run_mpsci_update(iroot) 
      
           open(unit=200,file="rdms_update_moit")
             string1="mv oneRDM_1PT oneRDM."//trim(ctmp2)//".moit"
             write(200,*)trim(string1)
             string1="mv twoRDM_1PT twoRDM."//trim(ctmp2)//".moit"
             write(200,*)trim(string1)
           close(200)
           call system("chmod +x rdms_update_moit")
           call system("./rdms_update_moit")

           oneRDMfile="oneRDM."//trim(ctmp2)//"."//"moit"
           twoRDMfile="twoRDM."//trim(ctmp2)//"."//"moit"

           D=0.0d0
           open(unit=110,file=trim(oneRDMfile))
             read(110,*)ij
             do i=1,nact
               do j=1,nact
                 read(110,*)ij,ji,D(i,j)
               end do
             end do
           close(110)
           !TM1=TM1+mat2%d*dmrg_weight(iroot+1)
           ! write(2,*)"print_mat D"
           ! call print_mat(nact,nact,D,2) 

           P=0.0d0
           open(unit=110,file=trim(twoRDMfile))
             READ(110,*)ijkl
             ijkl=0  ! Just a counter
             do
               read(110,*,iostat=ierr)ij,jk,kl,li,dv
               if(ierr.ne.0)exit
               ijkl=ijkl+1
               P(ij+1,jk+1,kl+1,li+1)=dv
               P(kl+1,li+1,ij+1,jk+1)=dv
               P(jk+1,ij+1,li+1,kl+1)=dv
               P(li+1,kl+1,jk+1,ij+1)=dv
             end do
           close(110)
           P=P*2.0d0
           !GM1=GM1+mat2%p*dmrg_weight(iroot+1)

      end Subroutine RDM_update_MPSci_LR 



