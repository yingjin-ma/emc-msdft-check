      Subroutine DMRG_CASSCF()
  
        use global_control
        use matrix

! Read in the orbital informations

! order in FCIDUMP
        allocate(iorder(nact)) ; iorder=0
        open(unit=110,file=dmrg_reorder_name)
          read(110,*)(iorder(i),i=1,nact)
!          write(*,*)(iorder(i),i=1,nact)
        close(110)
!write(*,*)dmrg_reorder_name
!stop
        write(*,*)""
        write(*,*)"MCSCF Part (orbital optimization):"
        write(*,*)"--------------------------------------------"
        write(*,"(A20)",advance='no')"    frozen orbitals:"
        do i=1,orb%nsub 
          write(*,"(I3)",advance='no')orb%frozen(i)
        end do
        write(*,*)""
        write(*,"(A20)",advance='no')"    closed orbitals:"
        do i=1,orb%nsub
          write(*,"(I3)",advance='no')orb%closed(i)
        end do
        write(*,*)""
        write(*,"(A20)",advance='no')"  occupied orbitals:"
        do i=1,orb%nsub
          write(*,"(I3)",advance='no')orb%   occ(i)
        end do
        write(*,*)""
        write(*,*)"--------------------------------------------"
        write(*,"(A20)",advance='no')"    active orbitals:"
        do i=1,orb%nsub
          write(*,"(I3)",advance='no')orb%act(i)
        end do
        write(*,*)""
        write(*,*)"--------------------------------------------"
        write(*,"(A20)",advance='no')"  external orbitals:"
        do i=1,orb%nsub
          write(*,"(I3)",advance='no')orb%extnl(i)
        end do
        write(*,*)""
        write(*,*)"--------------------------------------------"
        write(*,*)""
        !write(*,*)" DMRG-CASSCF strategy : ",method
        if(method(1).ne."WMKUBAR")then  
          write(*,*)" Iter | -- RDM_ENERGY -- |  sum|R| |  Method " 
        else
          write(*,*)" Iter |   DMRG/CI-E(0)  |     E(2)_SCF    |    &
              E(2)_VAL    |   sqrt(\sum(T**2)) "
        end if  

! Globle control in DMRG-CASSCF iterations

        allocate(threshold(ncycle))
        allocate(energies(ncycle))
        allocate(Rabs(ncycle))
        allocate(Udiis(ncycle))   

        ! initial in whole SCF
        allocate(T_ini(norb,norb)) ;           T_ini=0.0d0            !define
        allocate(U_ini(norb,norb,norb,norb)) ; U_ini=0.0d0

        threshold=0.0d0
        energies=0.0d0
        Rabs=0.0d0

! back up the orbital information for the later usage (cls <-> fro)
        orbbk%frozen=orb%frozen !backup
        orbbk%closed=orb%closed
        orbbk%   occ=orb%   occ
        orbbk%   act=orb%   act
        orbbk%  act2=orb%  act2
        orbbk% extnl=orb% extnl
        orbbk% total=orb% total
        orbbk%   ele=orb%   ele 


      end Subroutine DMRG_CASSCF

! =====================================================================

      Subroutine PreOneSCF(icycle,DFT_Fc,ifDFT)
  
        use global_control
        use matrix

        integer::icycle
        character*72 oneRDMfile,twoRDMfile
        double precision T1(norb,norb)
        double precision T2(norb,norb)
        double precision H_temp(norb,norb,norb,norb)
        double precision Hdiag_temp(norb,norb)
        double precision Gc(norb,norb)

        double precision dtmp,dtmp2,dij,dik,dil,djl,dv
        double precision RDM_IJ,RDM_KL,RDM_IJKL,RDM_J,RDM_J1
        double precision RDM_ENERGY_MAT
        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:) 
        double precision,allocatable::TM4(:,:)

        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

        character    ctmp1
        character*2  ctmp2         
        character*192 string0,string1,string2,string3,string4

!       exchange-correlation term
        double precision::DFT_Fc(norb,norb)
        logical::ifdft

!        write(*,*)"In Preone SCF ",iorder,nact
        call flush(6)
!        stop

!  ==========================================================================
!  ==============     Active the closed shell approximatation    ============
!  ==============        Modifying the RDMs (can be omited)      ============
!  ==========================================================================

!        call print_mat(nocc,nocc,mat2%d)
!        call print_GAT(nocc,nocc,nocc,nocc,mat2%p,11)
!        stop 

        allocate(mat1(orb%nsub)) 
        !mat1=0.0d0;
        j0=0; k0=0
        do i=1,orb%nsub
!          write(6,*)"==> orb%act(i)",orb%act(i)
          allocate(mat1(i)%D(orb%act(i),orb%act(i)))
          mat1(i)%D=0.0d0
          do j=1,orb%act(i)  
            do k=1,orb%act(i)
              mat1(i)%D(j,k)=mat2%d(j0+j,k0+k)
            end do
          end do
          j0=j0+orb%act(i)
          k0=k0+orb%act(i)
        end do
 
        if(.NOT.ALLOCATED(mat2%T))then
          allocate(mat2%T(norb,norb))
          mat2%T=0.0d0;
        end if
        allocate(mat2%U(norb,norb))
        allocate(mat2%R(norb,norb))
        allocate(mat2%A(norb,norb))
        allocate(mat2%B(norb,norb))
        allocate(mat2%A0(norb,norb))
        allocate(mat2%B0(norb,norb))
        allocate(mat2%Hdiag(norb,norb))
!        if(old_code)then
        allocate(mat2%G(norb,norb,norb,norb)) ! still use the old form
!        else
!        allocate(mat2%G(norb,nact,nact,norb)) ! partly reduced on 2015.3.16 (later)
!        end if        
        allocate(mat2%H(norb,norb,norb,norb)) !  

        allocate(mat2%Rocc(nact,nact))
        allocate(mat2%Aocc(nact,nact))
        allocate(mat2%Hocc(nact,nact,nact,nact)) !  

        mat2%U=0.0d0; mat2%R=0.0d0
        mat2%B=0.0d0; mat2%G=0.0d0; mat2%Hdiag=0.0d0
        mat2%H=0.0d0; mat2%A=0.0d0; mat2%Hocc=0.0d0
        mat2%Aocc=0.0d0; mat2%Rocc=0.0d0

! MATRIX A=hD+sum_kl[J^(kl)P^(lk)]
        do i=1,orb%nsub
          allocate(mat1(i)%A(orb%total(i),orb%total(i)))
          allocate(mat1(i)%B(orb%total(i),orb%total(i)))
          allocate(mat1(i)%hD(orb%total(i),orb%total(i)))
          allocate(mat1(i)%XC(orb%total(i),orb%total(i)))
!          allocate(mat1(i)%JP_KL(orb%total(i),orb%total(i)))
          mat1(i)%A=0.0d0
          mat1(i)%B=0.0d0
          mat1(i)%hD=0.0d0
          mat1(i)%XC=0.0d0
!          mat1(i)%JP=0.0d0
        end do

        write(6,*)" === Start to calculate RDM_energy === "
        call flush(6)
!        stop

        RDM_ENERGY=0.0d0      ! initial version
!        RDM_ENERGY_MAT=0.0d0  ! MKL (dgemm) version
        if(.not. ifcore) then

        if(old_code)then ! === old code ===
        else  ! === new code ===
          ioffset=0
          do i=1,orb%nsub 
            N=orb%occ(i)
            if(N.ne.0)then
              M=ioffset 
!              write(*,*)"loop of RDM_ENERGY "
!              call print_mat(N,N,mat1(i)%D(1:N,1:N),6)
!              call flush(6)
!              call print_mat(N,N,T(M+1:M+N,M+1:M+N),6)              
!              call flush(6) 
              call MXM(N,T(M+1:M+N,M+1:M+N),mat1(i)%D,mat1(i)%hD(1:N,1:N))
!              call print_mat(N,N,mat1(i)%hD(1:N,1:N),6)              
!              call flush(6) 
              call trace(N,mat1(i)%hD(1:N,1:N),dtmp)
              RDM_ENERGY=RDM_ENERGY+dtmp
!              write(*,*)"RDM_ENERGY in step",i,"is",RDM_ENERGY
!              write(*,*)"RDM_ENERGY_MAT in step",i,"is",RDM_ENERGY_MAT
            end if 
            ioffset=ioffset+orb%total(i)
          end do 

          if(ifDFT)then
            dtmp=0.0d0
            call print_mat(norb,norb,DFT_Fc,10001)
            write(6,*)"XC correlation part is involved"
            ioffset=0
            do i=1,orb%nsub 
              N=orb%occ(i)
              if(N.ne.0)then
                M=ioffset
                do i1=1,N
                  do j1=1,N 
                    dtmp=dtmp+DFT_Fc(M+i1,M+j1)
                    RDM_ENERGY=RDM_ENERGY+DFT_Fc(M+i1,M+j1)
!                    write(6,*)"M+i1,M+j1",M+i1,M+j1,DFT_Fc(M+i1,M+j1)
                  end do
                end do
              end if 
              ioffset=ioffset+orb%total(i)
            end do 
            write(6,*)"XC correlation energy is", dtmp 
          end if


        end if ! === 
         
          !write(6,*)"RDM_ENERGY 1-body ",RDM_ENERGY
          !write(6,*)"RDM_ENERGY_MAT",RDM_ENERGY_MAT
          !call flush(6)

         !call print_mat(norb,norb,T,6)       
         !call print_mat(nact,nact,mat2%d,6)

         ! write(6,*)"RDM_ENERGY 1-body ",RDM_ENERGY
         ! write(6,*)"RDM_ENERGY_MAT",RDM_ENERGY_MAT
         !stop
          !RDM_ENERGY=0.0d0

         if(old_code)then !== old code ===
 
         else ! === new code === 

! Will be used later
          !dtmp2=0.0d0
          ioffset=0
          ioffset1=0
          do i=1,orb%nsub; do i1=1,orb%occ(i)
            joffset=0
            joffset1=0
            do j=1,orb%nsub; do j1=1,orb%occ(j)
              Mk=0  ! MK is the offset of total orbital
              Nk=0  ! Nk is the offset of occupied orbital
              do k=1,orb%nsub
                Ml=0  ! Ml is the offset of total orbital
                Nl=0  ! Nl is the offset of occupied orbital
                do l=1,orb%nsub
                  if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
                    N1=orb%occ(k) 
                    N2=orb%occ(l)
                    if(N1.ne.0.and.N2.ne.0)then 
                      allocate(TM1(N1,N2));TM1=0.0d0 
                      allocate(TM2(N1,N2));TM2=0.0d0
                      allocate(TM3(N1,N2));TM3=0.0d0
                      TM1=U(ioffset+i1,joffset+j1,Mk+1:Mk+N1,Ml+1:Ml+N2)
                      TM2=mat2%P(ioffset1+i1,Nk+1:Nk+N1,Nl+1:Nl+N2,joffset1+j1)
!                     if(N1.eq.N2)then
!                       call MXM(N1,TM1,TM2,TM3)
!                       call trace(N1,TM3,dtmp) 
!                     else
                      dtmp=0
                      do k1=1,N1
                        do l1=1,N2
                          dtmp=dtmp+TM1(k1,l1)*TM2(k1,l1)
                        end do
                      end do 
!                     end if
!                     if(j.eq.7)then
!                       write(*,*)"i,j,k,l",i,j,k,l
!                       call print_MAT(N1,N2,mat2%P(ioffset1+i1,Nk+1:Nk+N1,Nl+1:Nl+N2,joffset1+j1))
!                       write(*,*)"i,k,l,j",i,k,l,j
!                       call print_MAT(N1,N2,TM1)
!                       call print_MAT(N1,N2,TM2)        
!                     end if 
                      RDM_ENERGY=RDM_ENERGY+dtmp*0.5d0
                      deallocate(TM1)
                      deallocate(TM2)
                      deallocate(TM3)
                    end if
                  end if
                  Ml=Ml+orb%total(l)
                  Nl=Nl+orb%occ(l)
                end do ! for l
                Mk=Mk+orb%total(k)
                NK=NK+orb%occ(k)
              end do ! for k 
!              write(*,*)i,j,"RDM_KL ",i1,j1,dtmp2/2.0d0 
            end do 
            joffset=joffset+orb%total(j)
            joffset1=joffset1+orb%occ(j)
            end do
          end do 
          ioffset=ioffset+orb%total(i)
          ioffset1=ioffset1+orb%occ(i)
          end do
          !call print_gat(nact,nact,nact,nact,mat2%P,32)
 
         end if ! === 

!         call print_GAT(norb,norb,norb,norb,U,11)
!         call print_GAT(nact,nact,nact,nact,mat2%p,12)

!          write(*,*)"2-RDM RDM_J   energy : ",RDM_J
!          write(*,*)"2-RDM RDM_J1 (total) : ",RDM_J1

!          write(*,*)"2-RDM RDM_IJKL energy : ",RDM_IJKL
!          write(*,*)"2-RDM energy (total) : ",dtmp2
            
!          stop

          write(*,*)"FRONRE0,FRONRE",FRONRE0,FRONRE
          write(*,*)"RDM2",RDM_ENERGY+FRONRE0
          write(*,*)"RDM2",RDM_ENERGY+FRONRE
!          write(*,*)"RDM 2-e",RDM_ENERGY
          !RDM_ENERGY=0.0d0
!          stop

2211      format(A7,I3,I3,I3,I3,G21.12)          

! The hD part
         if(old_code)then ! === old code ===

        else ! === new code ===

! The hD part -- MKL form
          noffset=0 
          do n=1,orb%nsub
            n1=orb%total(n)
            n2=orb%occ(n)
            if(n2.ne.0)then
              n0=noffset 
              allocate(TM1(N1,N2));TM1=0.0d0 
              allocate(TM2(N2,N2));TM2=0.0d0
              allocate(TM3(N1,N2));TM3=0.0d0
              mat1(n)%hD=0.0d0
              TM1=T(n0+1:n0+n1,n0+1:n0+n2)
              TM2=mat1(n)%D
              call MXMG(n1,n2,n2,TM1,TM2,TM3,"NN")
              TM3=TM3*2.0d0
!              mat1(n)%hD(1:n1,1:n2)=TM3
!              call print_mat(n1,n1,mat1(n)%hD)
              mat2%A(n0+1:n0+n1,n0+1:n0+n2)=TM3
              ! IF active MSDFT
              if(ifDFT)then
                write(6,*)"The xc term is covered in orb-grad"
                mat2%A(n0+1:n0+n1,n0+1:n0+n2)&
               =mat2%A(n0+1:n0+n1,n0+1:n0+n2)&
               +DFT_Fc(n0+1:n0+n1,n0+1:n0+n2)
              end if
              deallocate(TM1)
              deallocate(TM2)
              deallocate(TM3)
            end if 
            noffset=noffset+orb%total(n)
          end do

        end if ! ===

!         call print_mat(norb,norb,mat2%A,6)
!          mat2%A=0.0d0

!          write(*,*)"RDM_ENERGY",RDM_ENERGY
!          write(*,*)"RDM_ENERGY_MAT",RDM_ENERGY_MAT
          !RDM_ENERGY=0.0d0
!          ioffset=0
!          ioffset1=0
!          RDM_IJKL=0.0d0
!          do i=1,orb%nsub; do i1=1,orb%occ(i)
            !dtemp=0.0d0
            !write(*,*)"RDM2-I",RDM_ENERGY+FRONRE
            !write(2,*)"RDM2-I",RDM_ENERGY+FRONRE

!          call Print_MAT(norb,norb,mat2%A)

!          do i=1,norb
!            do j=1,norb
!              write(*,*)i,j,mat2%A(i,j)
!            end do 
!          end do
!          stop 
 
!          write(202,*)mat2%A
!          stop

!          mat2%A=0.0d0 ! for testing

! The sum_kl[J^(kl)P^(lk)] part         
       
         if(old_code)then ! === old code === 

! The sum_kl[J^(kl)P^(lk)] part --- mkl part,done, will be used later
         else ! === new code ===
        
          allocate(TM4(norb,norb));TM4=0.0d0
          koffset=0
          koffset1=0
          do k=1,orb%nsub; do k1=1,orb%occ(k) 
            loffset=0
            loffset1=0
            do l=1,orb%nsub; do l1=1,orb%occ(l)
              noffset=0 
              do n=1,orb%nsub
                n0=noffset 
                n1=orb%total(n)
                n2=orb%occ(n)
                moffset=0
                moffset1=0
                do m=1,orb%nsub
                  m0=moffset
                  m1=orb%total(m)
                  m2=orb%occ(m)
                  mx=moffset1
!                  if(orb%grouptable(n,m).eq.orb%grouptable(k,l))then
                  joffset=0
                  joffset1=0
                  do j=1,orb%nsub
                    if(orb%occ(j).ne.0.and.m2.ne.0)then
                      if(orb%grouptable(n,j).eq.orb%grouptable(k,l))then
                        j0=joffset
                        j1=orb%total(j) 
                        j2=orb%occ(j) 
                        jx=joffset1
                        allocate(TM1(N1,J2));TM1=0.0d0 
                        allocate(TM2(J2,M2));TM2=0.0d0
                        allocate(TM3(N1,M2));TM3=0.0d0               
                        TM1=U(n0+1:n0+n1,j0+1:j0+j2,k1+koffset,l1+loffset) 
                        TM2=mat2%p(mx+1:mx+m2,k1+koffset1,l1+loffset1,jx+1:jx+j2)
!                          TM2=mat2%p(mx+1:mx+m2,l1+loffset1,k1+koffset1,jx+1:jx+j2)
!                        if(n.eq.1.and.m.eq.1)then
!                          write(*,*)"========k,l,n,m,j==========",k,l,n,m,j
!                          call print_MAT(n1,j2,TM1)  
!                          write(*,*)"        ---       "
!                          call print_MAT(j2,m2,TM2) 
!                        end if
                        call MXMG(n1,m2,j2,TM1,TM2,TM3,"NT") 
                        TM3=TM3*2.0d0 
                        TM4(N0+1:N0+N1,M0+1:M0+M2)=TM4(N0+1:N0+N1,M0+1:M0+M2)+TM3
                        deallocate(TM1)
                        deallocate(TM2)
                        deallocate(TM3)                    
                      end if
                    end if
                    joffset=joffset+orb%total(j)
                    joffset1=joffset1+orb%occ(j)
                  end do 
                  !end if
                  moffset=moffset+orb%total(m)
                  moffset1=moffset1+orb%occ(m)
                end do
                noffset=noffset+orb%total(n)
              end do
            end do
            loffset=loffset+orb%total(l)
            loffset1=loffset1+orb%occ(l)
            end do
          end do 
          koffset=koffset+orb%total(k)
          koffset1=koffset1+orb%occ(k)
          end do
          mat2%A=mat2%A+TM4
          deallocate(TM4)
         end if ! ===
         !write(*,*)" ====A(1)/A(2)== "
         !call print_MAT(norb,norb,MAT2%A,6)

        ! Whole space to subspace
          ioffset=0 
          do i=1,orb%nsub
            do k=1,orb%total(i)
              do l=1,orb%total(i)
                if(dabs(mat2%A(k+ioffset,l+ioffset)).lt.1.0e-9)then
                   mat2%A(k+ioffset,l+ioffset)=0.0d0
                end if
                mat1(i)%A(k,l)=mat2%A(k+ioffset,l+ioffset)
              end do
            end do
            ioffset=ioffset+orb%total(i)
          end do

!          do i=1,orb%nsub
!            call print_MAT(orb%total(i),orb%total(i),MAT1(i)%A)
!          end do
!          stop 

          do i=1,norb
            do j=1,norb
              if(dabs(mat2%A(i,j)).lt.1.0d-9) mat2%A(i,j)=0.0d0
              if(mat2%A(i,j)-mat2%A(j,i).gt.1.0d-6) then
                write(1,*)i,j,"i,j",mat2%A(i,j)-mat2%A(j,i)
              end if
            end do 
          end do

!=====================================================================
! Also the sum_kl[J^(kl)P^(lk)] part, but not success, omit

          ! This is the mcscf Fcok matrix 
!          do i=1,orb%nsub
!            do k=1,orb%total(i)
!              do l=1,orb%total(i)           
!                write(*,*)"Fock",i,l,k,mat1(i)%A(l,k)/2.0d0
!              end do
!            end do
!          end do
!          stop


! The B=A+hTD+sum_(kl)[J^(kl)TP^(lk)+2K^(kl)TQ^(lk)]
! The G=hD+sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]

!     G="hD"+
!     B="hTD"+
        if(old_code)then ! === old code ===

! The G=hD+sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)] --  MKL & point-group
! only hD part, i.e. h*D_ij
        else ! === new code ===

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
        end if ! ===

!         call print_GAT(norb,nact,nact,norb,mat2%G,7)
!         stop

!         do i=1,norb
!           do j=1,norb
!              do k=1,norb
!                do l=1,norb 
!             write(*,*)i,j,mat2%B(i,j)
!                  write(4444,*)i,j,k,l,mat2%G(i,j,k,l)
!                end do
!              end do
!           end do
!         end do
!         stop
!       mat2%G=0.0d0  ! for testing
!      G=hD+"sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"
!      B=A+hTD+"sum_(kl)[J^(kl)TP^(lk)+2K^(kl)TQ^(lk)]"
        if(old_code)then ! === old code ===

!      G=hD+"sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"
!      sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"  part ---- mkl version  
       ! stage-1 :  matrix form ready, partly point-group 
       ! stage-2 :  dgemm and fully point group (later) since it will be too time-consuming
        else ! === new code ===      
 
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
!         call print_GAT(norb,nact,nact,norb,GM1,7)
         mat2%G=mat2%G+GM1
         deallocate(GM1)
        end if ! ===

!        call print_GAT(norb,norb,norb,norb,mat2%G,8)
!        stop

         !write(201,*)mat2%P
         !write(202,*)U
!         do i=1,norb
!           do j=1,norb
!             do k=1,norb
!               do l=1,norb
!                 write(10,*)i,j,k,l,mat2%G(i,j,k,l)
!               end do
!             end do 
!           end do
!         end do

!         stop  
 
         if(method(icycle).eq."microit")then
           call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
!           call print_mat(norb,norb,mat2%Hdiag)
!           stop
!           do i=1,norb
!             do j=1,norb
!               write(700,*)"Hdiag",i,j,mat2%Hdiag(i,j)   
!               write(3,*)i,j,mat2%A(i,j)
!             end do
!           end do
         end if

! Too lazy, actually only diagional should be enough in this case (also microit)
         if(method(icycle).eq."cp-inte")then
           call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag) 
         end if

         !call Hessian_test(norb,A,mat2%G,mat2%H,mat2%Hdiag)
         !stop

        else  ! This is a very old version!!!!! need to be rewritten 


        end if

          old_code=.false.   ! temporary active it 

! the rdm energy, inaccurate if RDM (old version of Block)
! QCMaquis is always fine
          write(*,2433,advance='no')"  ",icycle
          if(ifcore2)then
            write(*,2434,advance='no')"    ",rdm_energy+fronre0
          else
            write(*,2434,advance='no')"    ",rdm_energy+fronre
          end if
2433      format(a2,i3)
2434      format(a4,f14.8)
!          stop
 

          energies(icycle)=rdm_energy+fronre0
!          Echeck=rdm_energy+fronre0
          E0_micro=rdm_energy
          if(icycle.gt.1)then
            Echeck=dabs(energies(icycle)-energies(icycle-1))
          else
            Echeck=dabs(energies(icycle))
          end if

!          call print_mat(norb,norb,mat2%A,6)
! allocate the Umatrix for the next macro-iteration
          allocate(UMAT_FINAL(norb,norb)); UMAT_FINAL=0.0d0
          do i=1,norb
            UMAT_FINAL(i,i)=1.0d0
          end do 
!          call print_mat(norb,norb,UMAT_FINAL,6)

!          write(*,*)"energies : ",energies(1:icycle)

      end subroutine preonescf

! -------------------------------------------------
!  Use it to deallocate the allocated parameters
      Subroutine deallocate_preonescf()

        use global_control
        use matrix

        deallocate(mat1)  
        deallocate(mat2%d)  
        deallocate(mat2%p)  
!        deallocate(mat2%T)
        deallocate(mat2%U)
        deallocate(mat2%R)
        deallocate(mat2%A)
        deallocate(mat2%B)
        deallocate(mat2%A0)
        deallocate(mat2%B0)
        deallocate(mat2%Hdiag)
        deallocate(mat2%G)
        deallocate(mat2%H)
        deallocate(mat2%Rocc)
        deallocate(mat2%Aocc)
        deallocate(mat2%Hocc)
        deallocate(mat2%Ederi)
        deallocate(mat2%deltR)
        deallocate(mat2%deltT)

      End Subroutine deallocate_preonescf


