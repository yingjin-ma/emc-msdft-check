
      Subroutine Operator_update(icycle,CIupdate,H2update,iloop)
  
        use global_control
        use matrix
        use date_time

        integer::icycle
        logical::CIupdate,H2update
        character*72 oneRDMfile,twoRDMfile
        double precision T1(norb,norb)
        double precision T2(norb,norb)
        double precision H_temp(norb,norb,norb,norb)
        double precision Hdiag_temp(norb,norb)
        double precision Gc(norb,norb)

        double precision dtmp,dtmp2,dij,dik,dil,djl,dv
!        double precision RDM_IJ,RDM_KL,RDM_IJKL,RDM_J,RDM_J1
        double precision RDM_ENERGY_updated
!        double precision,allocatable::TM0(:,:)
        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)

        double precision,allocatable::TN1(:,:),TN2(:,:),TN3(:,:)
        double precision,allocatable::TN4(:,:),TN5(:,:),TN6(:,:)

        double precision,allocatable::TA1(:,:),TA2(:,:),TA3(:,:)
        double precision,allocatable::TA4(:,:),TA5(:,:),TA6(:,:)

        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

! For the 2-index transformed integrals
        type::transform
          double precision,allocatable::U(:,:)
          double precision,allocatable::T(:,:)
          double precision,allocatable::dR(:,:)
          double precision,allocatable::UdR(:,:)
        end type transform
        type(transform),allocatable::UM(:)

        double precision Tact(nact,nact)
        double precision Uact(nact,nact,nact,nact)

        double precision,allocatable::dT1dx(:,:)
        double precision,allocatable::dJ1dx(:,:,:,:)
        double precision,allocatable::dK1dx(:,:,:,:)
        double precision UdR(norb,norb)
        double precision deltaB(norb,norb)

        double precision U_updated(norb,norb)
        double precision U_tmp(norb,norb)
        double precision T_tmp(norb,norb)

        character    ctmp1
        character*2  ctmp2         
        character*192 string0,string1,string2,string3,string4

        allocate(dT1dx(norb,norb)) ; dT1dx=0.0d0
        allocate(dJ1dx(norb,norb,norb,norb)) ; dJ1dx=0.0d0
        allocate(dK1dx(norb,norb,norb,norb)) ; dK1dx=0.0d0

        deltaB=0.0d0
        Tact=0.0d0
        Uact=0.0d0

        U_tmp=mat2%U
        T_tmp=mat2%T

! check the current tr[T^+(A+B)]
!        allocate(TM1(norb,norb));TM1=0.0d0
!        allocate(TM2(norb,norb));TM2=0.0d0
!        TM1=0.0d0 
!        TM1=mat2%A+mat2%B
!        call MdXM(norb,mat2%T,TM1,TM2)
!        dv=0.0d0
!        call trace(norb,TM2,dv)
!        dv=dv*0.5d0
!        print *,"The Tr[A+B]) value : ",dv
!        print *,"The E2(T)    value : ",dv+RDM_energy+FRONRE0 
!        deallocate(TM1)
!        deallocate(TM2) 

!       allocate(TM(norb,norb));   TM=0.0d0
! Active Integrals for T and U(V) 
!       allocate(Tact(nact,nact)); Tact=0.0d0     
!       allocate(Uact(nact,nact,nact,nact));  Uact=0.0d0   
!       call print_mat(norb,norb,TRANS,6)


        if(H2update.or.CIupdate)then
          write(2,*) " ==============================================="
          write(2,*) "  == Hami will be updated in this micro-iter =="
          write(2,*) " ==============================================="
!          write(2,*) "          allow the H(2) update"
          call MXM(norb,UMAT_FINAL,mat2%U,U_updated)
          UMAT_FINAL=U_updated
!          print *,"updated the UMAT_FINAL"
!          call print_mat(norb,norb,UMAT_FINAL,6)
          call transform_2index_ALL&
               (orb%nsub,orb%total,norb,UMAT_FINAL,&
                T_origin,U_origin,T2H,U2H,orb%grouptable)
!          print *,"Finished the transform_2index_ALL"
!          call print_mat(norb,norb,T_origin,101)
!          call print_gat(norb,norb,norb,norb,U_origin,102)
!          call print_mat(norb,norb,T2H,6)
!          call print_gat(norb,norb,norb,norb,U2H,112)
          !stop
          T=T2H 
          U=U2H
          mat2%U=0.0d0
          mat2%T=0.0d0
          do i=1,norb 
            mat2%U(i,i)=1.0d0
          end do
          mat2%deltR=0.0d0
        else
          T=T2H 
          U=U2H          
          !call MXM(norb,UMAT_FINAL,mat2%U,U_updated)
          !!UMAT_FINAL=mat2%U
        end if
!        write(6,*)" iorder ",iorder; call flush(6)
!        call print_mat(norb,norb,T,111)
!        call print_gat(norb,norb,norb,norb,U,112)

!       Do the DMRG calculation during the micro-iter.
        if(CIupdate)then 
 
!          write(6,*)"CIupdated is used"
!          call flush(6)
 
          if(CP_integrals)then
            if(.false.)then ! ?? 
             ! Tact is already transformed before
              T2nd=Tact
            ! Uact (not yet done)
              call U2nd_integrals(orb%nsub,orb%occ,orb%total,nocc,norb,&
                                  U,J2dx,K2dx,Uact,orb%grouptable)
              U2nd=Uact 
            else
              call transform_2index&
                (orb%nsub,orb%occ,orb%total,nocc,norb,&
                MAT2%U,T,U,T2nd,U2nd,orb%grouptable) 
            end if
!            write(6,*)T2nd 
!            call flush(6)
            call print_mat(nact,nact,T2nd,121)
            call print_gat(nact,nact,nact,nact,U2nd,122)
!            stop
!            write(6,*) "Befoer 2-ap int"; call flush(6)
            call CS_dim_fro2cls()   ! CCCC
!            write(6,*) "After 2-ap int"; call flush(6)

!            write(111,*)"The next output"
!            write(112,*)"The next output"
!            call print_mat(norb,norb,T,111)
!            call print_gat(norb,norb,norb,norb,U,112)
!            call print_gat(nact,nact,nact,nact,U2nd,6)
!            call flush(6) 

! to Stefan : Eq.59 Eq.60
            call integrals_active("FCIDUMP_ACTIVE_2nd",.true.)
            call system("mv FCIDUMP_ACTIVE_2nd FCIDUMP_ACTIVE")     
!            call system("cat FCIDUMP_ACTIVE")     
!            write(6,*) "After 2-ap int "; call flush(6)

! to Stefan : Eq.57 Eq.61
!             The DMRG solver during micro-iterations
!             The Full DMRG sweeps can be reduced to 1) partly sweep 2) edge MPO 3) 1st-order PT
!             (explained in P11 bottom and P12 top) 

            walltime(19) = wtime()
            if(ith_INTE.eq.2)then

              write(3,"(A)")&
              "  micro-iters is extended by updated RDMs "
              write(30,"(A)")&
              "  micro-iters is extended by updated RDMs "

              !call run_dmrg_direct(icycle)
              call run_dmrg(icycle,iloop)

              walltime(20) = wtime()
              call timing_reporter(3,"MPS update ext-micro (strat 2)",walltime(20)-walltime(19))

              write(2,*)" "
              write(2,*)"============================================"
              write(2,*)"  MPSci update by the DMRG in micro-iter"
              write(2,*)"============================================"
            else if(ith_INTE.eq.3)then
              write(3,"(A)") &
              "  micro-iters is extended by updated RDMs "
              write(30,"(A)") &
              "  micro-iters is extended by updated RDMs "

              call run_dmrg_direct(icycle,iloop)
              !call run_dmrg(icycle,iloop)

              walltime(20) = wtime()
              call timing_reporter(3,"MPS update ext-micro (strat 3)",walltime(20)-walltime(19))

              write(2,*)" "
              write(2,*)"============================================"
              write(2,*)" MPSci update by the edge-DMRG in micro-iter"
              write(2,*)"============================================"
            else if(ith_INTE.eq.1)then
              write(3,"(A)") &
              "  micro-iters is extended by updated RDMs "
              write(30,"(A)") &
              "  micro-iters is extended by updated RDMs "

              write(2,*)
              write(2,*)"============================================"
              write(2,*)"     MPSci update by the 1st-order PT"
              write(2,*)"============================================"
              write(2,*) 
              call run_dmrg_diag(icycle)
              !call run_dmrg(icycle)
              call mpsci_update(icycle)
              call run_mpsci_update(icycle)        
              call system("mv oneRDM_1PT oneRDM.0.0")
              call system("mv twoRDM_1PT twoRDM.0.0")

              walltime(20) = wtime()
              call timing_reporter(3,"MPS update ext-micro (strat 1)",walltime(20)-walltime(19))

            else if(ith_INTE.eq.0)then

              write(2,*)
              write(2,*)"============================================"
              write(2,*)"      NO MPSci update is implemented"
              write(2,*)"============================================"
              write(2,*)
            end if
            write(2,*)&
            "Extended micro-iter (DMRG/CI side) is implemented"
            call RDMS_read() !CCCC
            ! write(6,*)"after RDMs_read in Operator_update"; call flush(6)

            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
            call closed_shell_update(mat2%d,mat2%P,TM1,GM1,1.0d0)  ! CCCC
            deallocate(mat2%d,mat2%P)
            allocate(mat2%d(nocc,nocc))
            allocate(mat2%p(nocc,nocc,nocc,nocc))
            mat2%d=TM1
            mat2%p=GM1
            deallocate(TM1,GM1)
            
          end if       
        end if

!        write(6,*)" iorder ",iorder; call flush(6)
!        write(6,*)"after fro2cls",nocc,nact
!        call flush(6) 

! Distribute the transform matrix base on group
        allocate(UM(orb%nsub))
        i0=0
        call MXM(norb,mat2%U,mat2%deltR,UdR) 
        do i=1,orb%nsub
          itot=orb%total(i)
          allocate(UM(i)%U(itot,itot))
          allocate(UM(i)%T(itot,itot))
          UM(i)%U  = 0.0d0
          UM(i)%T  = 0.0d0
          UM(i)%U  = mat2%U(i0+1:i0+itot,i0+1:i0+itot)
          UM(i)%UdR=    UdR(i0+1:i0+itot,i0+1:i0+itot)
          UM(i)%T =UM(i)%U
!         T=U-1 in sub-irreps
          do j=1,itot
            UM(i)%T(j,j)=UM(i)%U(j,j)-1.0d0
          end do
!         call print_mat(tot(i),tot(i),UM(i)%U)
          i0=i0+orb%total(i)
        end do

! The N-index transform for h,J,K
        !   hU, hUdR 
        call MXM(norb,T,mat2%U,T1dx)
        call MXM(norb,T,UdR,dT1dx)
        !U^+hU
        allocate(TM1(norb,norb));TM1=0.0d0
        call MdXM(norb,mat2%U,T1dx,TM1)
        i0=0
        ia=0
        do i=1,orb%nsub
          ir=orb%total(i)
          ij=orb%occ(i)
          Tact(ia+1:ia+ij,ia+1:ia+ij)=TM1(i0+1:i0+ij,i0+1:i0+ij)
          i0=i0+orb%total(i)
          ia=ia+orb%occ(i)
        end do
        deallocate(TM1)
!        write(*,*)"hU and hUdR"

!        call print_gat(norb,norb,norb,norb,U,115)
!        stop

        ! JU & U^+JU     KT & T^+KT 
        k0=0
        ka=0
        do k=1,orb%nsub; do kk=1,orb%occ(k)
          l0=0
          la=0
          do l=1,orb%nsub; do ll=1,orb%occ(l)
            i0=0
            ia=0
            do i=1,orb%nsub
              j0=0
              ja=0
              do j=1,orb%nsub
                if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
                  Ni=orb%total(i)
                  Nj=orb%total(j)
                  if(Ni.ne.0.and.Nj.ne.0)then
    
                    !write(6,*)i,j,k,l

                    allocate(TM1(Ni,Nj));TM1=0.0d0
                    allocate(TM2(Ni,Nj));TM2=0.0d0
                    allocate(TM3(Ni,Nj));TM3=0.0d0
                    TM1=U(i0+1:i0+ni,j0+1:j0+nj,kk+k0,ll+l0)
                    ! 1-index for J  
                    call MXMG(ni,nj,nj,TM1,UM(j)%U,TM2,"NN")
                     J1dx(i0+1:i0+ni,j0+1:j0+nj,kk+k0,ll+l0)=TM2
                    ! 2-index for J
                    call MXMG(ni,nj,ni,UM(i)%U,TM2,TM3,"TN")
                     J2dx(i0+1:i0+ni,j0+1:j0+nj,kk+k0,ll+l0)=TM3
                    ! 1-index for delta J 
                    call MXMG(ni,nj,nj,TM1,UM(j)%UdR,TM2,"NN")
                    dJ1dx(i0+1:i0+ni,j0+1:j0+nj,kk+k0,ll+l0)=TM2
                    deallocate(TM1)
                    deallocate(TM2)
                    deallocate(TM3)

                    !write(6,*)"J transformed pass"
                    !call flush(6)

                    allocate(TM1(Ni,Nj));TM1=0.0d0
                    allocate(TM2(Ni,Nj));TM2=0.0d0
                    allocate(TM3(Ni,Nj));TM3=0.0d0
                    TM1=U(i0+1:i0+ni,kk+k0,ll+l0,j0+1:j0+nj)
                    ! 1-index for K  
                    call MXMG(ni,nj,nj,TM1,UM(j)%T,TM2,"NN")
                     K1dx(i0+1:i0+ni,kk+k0,ll+l0,j0+1:j0+nj)=TM2
                    ! 2-index for K
                    call MXMG(ni,nj,ni,UM(i)%T,TM2,TM3,"TN")
                     K2dx(i0+1:i0+ni,kk+k0,ll+l0,j0+1:j0+nj)=TM3
                    ! 1-index for delta K  
                    call MXMG(ni,nj,nj,TM1,UM(j)%UdR,TM2,"NN")
                    dK1dx(i0+1:i0+ni,kk+k0,ll+l0,j0+1:j0+nj)=TM2
                    deallocate(TM1)
                    deallocate(TM2)
                    deallocate(TM3)
 
                    !write(6,*)"K transformed pass"
                    !call flush(6)                    

                  end if
                end if
                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do
              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do
          end do
          l0=l0+orb%total(l)
          la=la+orb%occ(l)
          end do
        end do
        k0=k0+orb%total(k)
        ka=ka+orb%occ(k)
        end do
!        write(6,*)"J,K, transformed pass"
!        call flush(6)
!        stop 

!        call flush(6) 
!        call print_mat(nact,nact,mat2%D,21)
!        call print_gat(nact,nact,nact,nact,mat2%P,22)
!        stop

        if(.not.allocated(mat1))then
          allocate(mat1(orb%nsub)) 
        end if
        !mat1=0.0d0;
        j0=0; k0=0
        do i=1,orb%nsub
          if(.not.allocated(mat1(i)%D))then
            allocate(mat1(i)%D(orb%act(i),orb%act(i)))
            mat1(i)%D=0.0d0
          else
            mat1(i)%D=0.0d0
          end if
          do j=1,orb%act(i)  
            do k=1,orb%act(i)
              mat1(i)%D(j,k)=mat2%d(j0+j,k0+k)
            end do
          end do
          j0=j0+orb%act(i)
          k0=k0+orb%act(i)
        end do
 
        if(.NOT.ALLOCATED(mat2%A))then
          allocate(mat2%A(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%B))then
          allocate(mat2%B(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%A0))then
          allocate(mat2%A0(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%B0))then
          allocate(mat2%B0(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%Hdiag))then
          allocate(mat2%Hdiag(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%G))then
!         if(old_code)then
          allocate(mat2%G(norb,norb,norb,norb)) ! still use the old form
        end if
!        else
!        allocate(mat2%G(norb,nact,nact,norb)) ! partly reduced on 2015.3.16 (later)
!        end if        
        if(.NOT.ALLOCATED(mat2%H))then
          allocate(mat2%H(norb,norb,norb,norb)) !  
        end if

!        write(6,*)"after allocated"
!        call flush(6)


        if(CIupdate)then

          if(CP_integrals)then

!           Check the new RDM energy
            ! 1-ele
            RDM_energy_updated=0.0d0         
            allocate(TM1(nact,nact));TM1=0.0d0
            allocate(TM2(nact,nact));TM2=0.0d0
            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)
              TM1(ia+1:ia+ij,ia+1:ia+ij)=T2nd(ia+1:ia+ij,ia+1:ia+ij)
              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do       
            call MXM(nact,TM1,mat2%D,TM2)
            call trace(nact,TM2,dv)  
      
            RDM_energy_updated=dv
!            write(6,*)"energy 1-e",RDM_energy_updated
!            call flush(6)
     
            deallocate(TM1)
            deallocate(TM2)
     
!            write(6,*)"nact,nocc",nact,nocc
!           call flush(6) 

            ! 2-ele
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

                TM1=  U2nd(:,:,kk+ka,ll+la)  ! bk
                TM3=mat2%P(:,kk+ka,ll+la,:)  ! bk
 
!                call print_mat(nact,nact,TM1,6)
!                call flush(6)

                i0=0
                ia=0
                do i=1,orb%nsub
                  ir=orb%total(i)
                  ij=orb%occ(i)

                  j0=0
                  ja=0
                  do j=1,orb%nsub
                    jr=orb%total(j)
                    jj=orb%occ(j)                                       
                   TM2(ia+1:ia+ij,ja+1:ja+jj)=TM1(ia+1:ia+ij,ja+1:ja+jj)
 
                    j0=j0+orb%total(j)
                    ja=ja+orb%occ(j)                      
                  end do 
                  i0=i0+orb%total(i)
                  ia=ia+orb%occ(i)
                end do
!                write(*,*)" ===  TM1/TM2  === "
!                call print_mat(nact,nact,TM2,6)
!                call flush(6)
     
                call MXM(nact,TM2,TM3,TM4)     
                call trace(nact,TM4,dv) 

                RDM_energy_updated=RDM_energy_updated+dv*0.5d0
     
                deallocate(TM1)
                deallocate(TM2)
                deallocate(TM3)
                deallocate(TM4)
     
              end do
              l0=l0+orb%total(l)
              la=la+orb%occ(l)
              end do
            end do
            k0=k0+orb%total(k)
            ka=ka+orb%occ(k)
            end do

!            write(6,*)"updated RDM_energy(1+2)",RDM_ENERGY_updated
!            call flush(6)
            write(2,*)"updated RDM_energy(NRE)",RDM_ENERGY_updated+fronre0
            E0_micro=RDM_ENERGY_updated
!            call print_mat(nact,nact,T2nd,21)
!            call print_gat(nact,nact,nact,nact,U2nd,22)
!            call print_mat(nact,nact,mat2%D,23)
!            call print_gat(nact,nact,nact,nact,mat2%P/2.0d0,24)
!            stop
!            if(Echeck.gt.RDM_ENERGY_updated+fronre0)then
!              Echeck=RDM_ENERGY_updated+fronre0
!            else
!              Lexit=.true. 
!            end if
          end if
          !if(icycle.eq.1)stop
          !stop
        end if 

!       Construct the B matrix base on 1dx integrals

        mat2%B =0.0d0     
        mat2%A0=0.0d0     
        mat2%B0=0.0d0     
 
        allocate(TM1(norb,nact));TM1=0.0d0
        allocate(TM2(norb,nact));TM2=0.0d0
        allocate(TN1(norb,nact));TN1=0.0d0
        allocate(TN2(norb,nact));TN2=0.0d0
        allocate(TA1(norb,nact));TA1=0.0d0
        allocate(TA2(norb,nact));TA2=0.0d0
        i0=0
        ia=0 
        do i=1,orb%nsub
          ir=orb%total(i)
          ij=orb%occ(i)
          TM1(i0+1:i0+ir,ia+1:ia+ij)= T1dx(i0+1:i0+ir,i0+1:i0+ij)
          TN1(i0+1:i0+ir,ia+1:ia+ij)=dT1dx(i0+1:i0+ir,i0+1:i0+ij)
          TA1(i0+1:i0+ir,ia+1:ia+ij)=    T(i0+1:i0+ir,i0+1:i0+ij)
          i0=i0+orb%total(i)
          ia=ia+orb%occ(i)         
        end do

        call MXMG(norb,nact,nact,TM1,mat2%D,TM2,"NN")        
        call MXMG(norb,nact,nact,TN1,mat2%D,TN2,"NN")        
        call MXMG(norb,nact,nact,TA1,mat2%D,TA2,"NN")        

        i0=0
        ia=0 
        do i=1,orb%nsub
          ir=orb%total(i)
          ij=orb%occ(i)
           mat2%B(i0+1:i0+ir,i0+1:i0+ij)=TM2(i0+1:i0+ir,ia+1:ia+ij)
           deltaB(i0+1:i0+ir,i0+1:i0+ij)=TN2(i0+1:i0+ir,ia+1:ia+ij)
          mat2%A0(i0+1:i0+ir,i0+1:i0+ij)=TA2(i0+1:i0+ir,ia+1:ia+ij)
          i0=i0+orb%total(i)
          ia=ia+orb%occ(i)         
        end do
!        mat2%B=mat2%B*2.0d0  ! in order to match my old formula
        
!        call print_mat(norb,norb,Mat2%B,6)  ! checked
!        mat2%B=0.0d0

        deallocate(TM1)
        deallocate(TM2)
        deallocate(TN1)
        deallocate(TN2)
        deallocate(TA1)
        deallocate(TA2)

!       Continue to construct the B matrix base on J1dx and K1dx integrals

        allocate(TM5(norb,nact)); TM5=0.0d0
        allocate(TM6(norb,nact)); TM6=0.0d0
        allocate(TN5(norb,nact)); TN5=0.0d0
        allocate(TN6(norb,nact)); TN6=0.0d0
        allocate(TA5(norb,nact)); TA5=0.0d0

        k0=0
        ka=0
        do k=1,orb%nsub; do kk=1,orb%occ(k)
          l0=0
          la=0
          do l=1,orb%nsub; do ll=1,orb%occ(l)


! J^{kl}*P^{lk} 
            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=     U(:,:,kk+k0,ll+l0)
            TM3=mat2%P(:,kk+ka,ll+la,:)

            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jr=orb%total(j)
                jj=orb%occ(j)

              TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)

                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do

              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")

            TA5=TA5+TM4

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4)


! J^{kl}U*P^{lk} 
            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=J1dx(:,:,kk+k0,ll+l0)
            TM3=mat2%P(:,kk+ka,ll+la,:)

            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jr=orb%total(j)
                jj=orb%occ(j)

              TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)

                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do

              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")
 
            TM5=TM5+TM4 
 
            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)                   
            deallocate(TM4)                 

! J^{kl}*UdR*P^{lk}

            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=dJ1dx(:,:,kk+k0,ll+l0)
            TM3=mat2%P(:,kk+ka,ll+la,:)

            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jr=orb%total(j)
                jj=orb%occ(j)
                TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)
                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do

              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")

            TN5=TN5+TM4

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4)
 

! K^{kl}*Q^{lk} 
            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=K1dx(:,kk+k0,ll+l0,:)
            TM3=mat2%P(:,:,ll+la,kk+ka)

            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jr=orb%total(j)
                jj=orb%occ(j)

                TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)

                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do

              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")

            TM6=TM6+TM4 

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4)

! K^{kl}*UdR*Q^{lk} 
            allocate(TM1(norb,norb));TM1=0.0d0
            allocate(TM2(norb,nact));TM2=0.0d0
            allocate(TM3(nact,nact));TM3=0.0d0
            allocate(TM4(norb,nact));TM4=0.0d0
            TM1=dK1dx(:,kk+k0,ll+l0,:)
            TM3=mat2%P(:,:,ll+la,kk+ka)

            i0=0
            ia=0
            do i=1,orb%nsub
              ir=orb%total(i)
              ij=orb%occ(i)

              j0=0
              ja=0
              do j=1,orb%nsub
                jr=orb%total(j)
                jj=orb%occ(j)

                TM2(i0+1:i0+ir,ja+1:ja+jj)=TM1(i0+1:i0+ir,j0+1:j0+jj)

                j0=j0+orb%total(j)
                ja=ja+orb%occ(j)
              end do

              i0=i0+orb%total(i)
              ia=ia+orb%occ(i)
            end do

            call MXMG(norb,nact,nact,TM2,TM3,TM4,"NN")

            TN6=TN6+TM4

            deallocate(TM1)
            deallocate(TM2)
            deallocate(TM3)
            deallocate(TM4) 
 
          end do
          l0=l0+orb%total(l)
          la=la+orb%occ(l)
          end do
        end do
        k0=k0+orb%total(k)
        ka=ka+orb%occ(k)
        end do

        allocate(TM1(norb,norb));TM1=0.0d0
        allocate(TM2(norb,norb));TM2=0.0d0
        allocate(TN1(norb,norb));TN1=0.0d0
        allocate(TN2(norb,norb));TN2=0.0d0
        allocate(TA1(norb,norb));TA1=0.0d0

        i0=0
        ia=0
        do i=1,orb%nsub
          ir=orb%total(i)
          ij=orb%occ(i)
          TM1(i0+1:i0+ir,i0+1:i0+ij)=TM5(i0+1:i0+ir,ia+1:ia+ij)
          TM2(i0+1:i0+ir,i0+1:i0+ij)=TM6(i0+1:i0+ir,ia+1:ia+ij)
          TN1(i0+1:i0+ir,i0+1:i0+ij)=TN5(i0+1:i0+ir,ia+1:ia+ij)
          TN2(i0+1:i0+ir,i0+1:i0+ij)=TN6(i0+1:i0+ir,ia+1:ia+ij)
          TA1(i0+1:i0+ir,i0+1:i0+ij)=TA5(i0+1:i0+ir,ia+1:ia+ij)
          i0=i0+orb%total(i)
          ia=ia+orb%occ(i)
        end do       

        mat2%A0=mat2%A0+TA1
        mat2%A0=2.0d0*mat2%A0        

!        call print_mat(norb,norb,TM1+TM2,6)
        mat2%B=mat2%B+TM1+2.0d0*TM2   
        deltaB=deltaB+TN1+2.0d0*TN2   

        mat2%B=mat2%B*2.0d0  ! in order to match my old formula
!        deltaB=deltaB 
        mat2%B0=mat2%B 

!        mat2%B=Aini
 
!        mat2%U=U_tmp
!        mat2%T=T_tmp

!        write(*,*)" op ====A(1)/A(2)== op "
!        call print_mat(norb,norb,mat2%B,6)

        call MdXM(norb,mat2%U,mat2%B,mat2%A)

!        call print_mat(norb,norb,mat2%U,6)
!        write(*,*)" op ==== U // A = UB == op "
!        call print_mat(norb,norb,mat2%A,6)
    
        mat2%B=mat2%B+deltaB

!        call print_mat(norb,norb,deltaB,6)
 
        deallocate(TM1)
        deallocate(TM2)
        deallocate(TM5)
        deallocate(TM6)
        deallocate(TA5)

        deallocate(TN1)
        deallocate(TN2)
        deallocate(TN5)
        deallocate(TN6)

!        goto 2222

!  ==============
! The G=hD+sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)] --  MKL & point-group
       if(CIupdate.or.H2update)then 
         ! The G need to be updated if new RDMs or integrals(??) are obtained

          call Gmat_gen(orb%nsub,orb%act,orb%total,&
               nact,norb,T,U,mat2%D,mat2%P,mat2%G,orb%grouptable) 
               write(2,*)"G operator is also updated in DMRG/CI steps"

       end if 

!        call print_GAT(norb,norb,norb,norb,mat2%G,8)
!        write(6,*)"Finishing operator_update"
!        call flush(6) 

!2222     continue 

        deallocate(UM)
        deallocate(dT1dx)
        deallocate(dJ1dx)
        deallocate(dK1dx)

      end subroutine operator_update

! -------------------------------------------------
!  Use it to deallocate the allocated parameters
      Subroutine deallocate_operator_update()

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
        deallocate(mat2%Hdiag)
        deallocate(mat2%G)
        deallocate(mat2%H)
        deallocate(mat2%Rocc)
        deallocate(mat2%Aocc)
        deallocate(mat2%Hocc)
        deallocate(mat2%Ederi)
        deallocate(mat2%deltR)
        deallocate(mat2%deltT)

      End Subroutine deallocate_operator_update


