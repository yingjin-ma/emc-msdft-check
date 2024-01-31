! =====================================================
! Severl ways to use the coupling terms
! 1) origonal way without approximation
! 2) treat the MPS changes as the 1st order response 
! 3) extended iterative way 
! 4) 2nd Hamiltanion treatment way
! =====================================================

! 1) origonal way without approximation 
!   It should be used for gradients, but not for coupled orb-opt.
!     Since the \DeltaR also involved in orb-opt
      Subroutine coupling_solver(icycle)

        use global_control
        use coupling_terms
        use matrix

        integer icycle 
! This solver should be very similar to the Solver_Hessian
        double precision,allocatable::H(:,:)
        integer offset 

        character    ctmp1
        character*2  ctmp2
        character*2  ctmp3
        character*192 string0,string1,string2,string3,string4

! This solver should be very similar to the Solver_Hessian
        integer,allocatable::porb(:)
        double precision d0,dv,dv0,redundant
        double precision d0_sum0,d2_sum0
        double precision,allocatable::D(:),X(:),Y(:)
        double precision,allocatable::HP(:,:),HPA(:,:),XP(:),YP(:)
        double precision,allocatable::T1(:,:),T2(:,:)
        double precision,allocatable::  L1(:),  L2(:)

        integer,allocatable::valid(:,:)
        Logical Lvalid
        double precision lamda
        double precision grad_digit

        lamda=1.0d0 ! "1" means the standard AH

! Assuming the whole Hessian is needed by ndim1 = norb
        ndim1=norb
        ndim2=ndim1**2
        ndimH=ndim2+nC
        allocate(H(ndimH,ndimH))
        H=0.0d0

! Check the redundant rotations
        ! initialize of Hessian matrix that will be used in this solver
          allocate(X(ndimH)); X=0.0d0
          allocate(Y(ndimH)); Y=0.0d0

        ! MCSCF Fock/gradient matrix together with MPS_ci gradient part 

        ! 1) the orbital gradients
          redundant=1.0e-9
          ij=0
          nij=0
          grad_digit=0.0d0
          do i=1,norb
            do j=i+1,norb ! only half is needed
              ij=ij+1
              Y(ij)=mat2%A(i,j)-mat2%A(j,i)
              if(dabs(Y(ij)).lt.redundant)then
              else
                if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                  ! Check for degeneracy irreps
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(Y(ij)).gt.grad_digit)then
                      grad_digit=dabs(Y(ij))
                    end if
                  end if
                  ! end Check for degeneracy irreps
                end if
              end if
            end do
          end do

          write(1,*)"Checked for the no-redundant orbital-rotations 1/2"
          call flush(1)

          do i=1,norb
            do j=1,i ! only half is needed
              ij=ij+1
              Y(ij)=mat2%A(i,j)-mat2%A(j,i)
              if(dabs(Y(ij)).lt.redundant)then
              else
                if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                  ! Check for degeneracy irreps
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(Y(ij)).gt.grad_digit)then
                      grad_digit=dabs(Y(ij))
                    end if
                  end if
                  ! end Check for degeneracy irreps
                end if
              end if
            end do
          end do

          write(1,*)"Checked for the no-redundant orbital-rotations 2/2"
          call flush(1)

          nij0=nij
        ! re-run it to get the valid rotations index
          allocate(valid(nij,2)); valid=0
          ij=0
          kl=0
          do i=1,norb
            do j=i+1,norb  !half
              ij=ij+1
              if(dabs(Y(ij)).lt.redundant)then
              else
                IF(LMCSCF)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                  end if
                else
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      kl=kl+1
                      valid(kl,1)=i
                      valid(kl,2)=j
                    end if
                  end if
                end if
              end if
            end do
          end do
          do i=1,norb
            do j=1,i  !half
              ij=ij+1
              if(dabs(Y(ij)).lt.redundant)then
              else
                IF(LMCSCF)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                  end if
                else
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      kl=kl+1
                      valid(kl,1)=i
                      valid(kl,2)=j
                    end if
                  end if
                end if
              end if
            end do
          end do

          write(1,*)"valid rotations in coupled solver",nij
          call flush(1)
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,1)
          end do
          write(1,*)" "
          do i=1,nij
            write(1,"(I3)",advance='no')valid(i,2)
          end do
          call flush(1)

          nij0=nij
           nij=nij/2
          nijO=nij ! number of valid orbital rotation 

          ! 2) the MPS_ci gradients
          do iC=1,nC
            ij=ij+1
            !Y(ij)=2.0d0*Ci(iC)%grad
            Y(ij)=1.0d0*Ci(iC)%grad
            if(iC.eq.1)then ! Assuming it is the ground state
              Y(ij)=0.0  ! No rotation for its own
            end if

!            if(dabs(Y(ij)).lt.redundant)then
!           It should be activated (two-step WNR)
            if(dabs(Y(ij)).lt.0)then
            else
              nij=nij+1
              write(1,*) "MPSci",ic,"gradient",Y(ij)
!             print *, "MPSci",ic,"gradient set to 0"
!             Y(ij)=0.00000001
            end if
          end do
          nijC=nij-nijO ! number of valid MPS_ci rotation

          write(1,*)"nij0,nij",nij0,nij
          call flush(1)

          allocate(porb(nijO));   porb=0
          allocate(YP(nij));     YP=0.0d0
          allocate(HP(nij,nij)); HP=0.0d0

          ij=0
          kl=0
          do i=1,norb
            do j=1,norb ! only half is needed
              ij=ij+1
              Y(ij)=mat2%A(i,j)-mat2%A(j,i)
              !Check the valid rotations  
              Lvalid=.false.
              do k=1,nijO
                if(i.eq.valid(k,1).and.j.eq.valid(k,2))then
                  Lvalid=.true.
                  exit
                end if
              end do
              if(Lvalid)then  ! can be deleted
                kl=kl+1
                porb(kl)=ij
              end if
            end do
          end do

          write(1,*)"porb : "
          write(1,*) porb
          call flush(1)


        iloop=0 
        iloopmax=max_iter 
        do 
          iloop=iloop+1
          write(1,*)"The ",iloop,"-th micro-iteration "

          ! updated E(0) = v'Hv
          allocate(L1(nC))  ; L1=0.0d0
          call MXMG(1,nC,nC,Ci(:)%dv,Hami,L1,"NN")
          call MXMG(1,1,nC,L1,Ci(:)%dv,dv0,"NN")
!          write(1,*)" Debug uncomment the MXMG()above "
!          write(1,*)"The E(0) or v'Hv value is ",dv0+fronre0
          deallocate(L1)
!          call print_mat(nC,nC,Hami,1) 
!          call print_mat(1,nC,Ci(:)%dv,1) 
          call flush(1)
          ! updated E(2)=E(0)+1/2*Tr(A+B)
 
          ! Here should be 1-index updated operators
          ! Add later
          if(iloop.gt.1)then 
!            write(1,*)"Before the Operator_update_CPNR" 
!            call flush(1)

            call Operator_update_CPNR(iloop)

            call coupling_update(iloop)

            write(1,*)"After coupling_update" ,iloop
            call flush(1)

!            allocate(T1(norb,norb)); T1=0.0d0              
!            T1=mat2%A+mat2%B
!            call trace(norb,T1,dv)
!            write(1,*)"1/2*tr(A+B)",0.5d0*dv 
!            write(1,*)"E(2)",0.5d0*dv+dv0 
!            deallocate(T1)

! This is the delta-orbiatl-rotation gradients,
!       i.e. re-calculate the orbital residual
            ij=0
            kl=0
            do i=1,ndim1
              do j=1,ndim1 ! only half is needed
              !do j=1,ndim1 ! all
                 ij=ij+1
                 Y(ij)=(mat2%A(i,j)-mat2%A(j,i)) ! gradients
                 !Check the valid rotations  
              end do
            end do           
!       i.e. re-calculate the MPSci residual
!            to be added later



          end if 
!          if(iloop.eq.3) stop

          ! calculate the orb-orb Hessian
          call Hessian(ndim1,mat2%A,mat2%G,mat2%H,mat2%Hdiag)        

          write(1,*)"after generate the orb-orb Hessian"
          call flush(1)

          ! 1) Put orb-orb coupling
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

        ! offset value, used later
          offset=norb*norb

        ! 2) Put orb-MPS coupling
          do iC=1,nC 
            ij=0
            do i=1,ndim1
              do j=1,ndim1 
                ij=ij+1
                H(offset+iC,ij)=HCR(iC,i,j) 
                H(ij,offset+iC)=HCR(iC,i,j) 
              end do
            end do
          end do

        ! 3) Put the MPS-MPS coupling
          do i=1,nC
            do j=1,nC
              H(offset+i,offset+j)=HCC(i,j) 
            end do
          end do 
 
          call print_mat(ndimH,ndimH,H,716)
          call flush(716)

          write(1,*)"after print the full Hessian"
          call flush(1)

       ! valid orb-orb part
          grad_digit=0.0d0
          do i=1,nijO
            YP(i)=Y(porb(i))
            if(dabs(YP(i)).gt.grad_digit)then ! orbital residual
              grad_digit=dabs(YP(i))
            end if
!            if(i.gt.nijO)then ! treat Grad(mps_ci) as 0
!            YP(i)=0.0d0
!          end if 
            do j=1,nijO
              HP(i,j)=H(porb(i),porb(j))
            end do
          end do

        ! valid mps-orb coupling
          do i=1,nijC
            do j=1,nijO
              HP(nijO+i,j)=H(ndim2+i,porb(j))
              HP(j,nijO+i)=H(porb(j),ndim2+i)
            end do
          end do      

        ! MPSci gradients
          do iC=1,nC
            YP(iC+nijO)=1.0d0*Ci(iC)%grad
          end do
 
          write(1,*)"ndim2 should be 100, and it is ", ndim2
          write(1,*)" nij0/nijC/nij should be ", nijO,nijC,nij
          write(1,*)" Current residual ", YP     
          call flush(1)

        ! valid mps-mps part 
          do i=1,nijC
            do j=1,nijC
              HP(nijO+i,nijO+j)=H(ndim2+i,ndim2+j)
            end do
          end do 
 
          if(.true.)then
            allocate(T1(nij,nij))
            allocate(L1(nij))
            T1=HP
            write(2,*)"The valid Hessian"
            call flush(2)
            call print_mat(nij,nij,T1,2)
            call eigtql2(nij,T1,L1)
            write(2,*)"The  eigen-values"
            call flush(2)
            write(2,*)L1
            write(2,*)"The eigen-vectors"
            call print_mat(nij,nij,T1,2)
            call flush(2) 
            deallocate(T1)
            deallocate(L1)

            write(*,*)
            write(*,*)"==============================="
            write(*,*)" Stop just after Hessian Check"
            write(*,*)"==============================="
            stop

          end if
  
          if(.true.)then
          ! Augmented Hessian Solver
            allocate(HPA(nij+1,nij+1))
            HPA=0.0d0

          ! The new augmented Hessian matrix {except HPA(1,1)}
            HPA(2:nij+1,1)=YP
            HPA(1,2:nij+1)=YP
            HPA(2:nij+1,2:nij+1)=HP
          ! And the augmented x-vector
            allocate(XP(nij+1))
            allocate(D(nij+1))
            XP=0.0d0
             D=0.0d0

            write(1,*)"The effective coupled Hessian : "
            call print_mat(nij+1,nij+1,HPA,1)

            if(.false.)then
              call Davidson(nij+1,HPA,XP,D,lamda,.true.,nij+1,1.0d-12)
            else
! to Stefan : Eq.51
!             Solver equation for coupled NR formula (Augmented Hessian algorithm)
              call eigtql2(nij+1,HPA,XP)
              write(2,*)"The AH  eigen-values"
              write(2,*)XP
              write(2,*)"The AH eigen-vectors"
              call print_mat(nij+1,nij+1,HPA,2)
              XP=HPA(:,1)
              call flush(2)
            end if

            if(XP(1).lt.0)then
              XP=-1.0d0*XP
            end if

            write(1,*) "XP", XP
            call flush(1)

!          do i=1,nij
!            X(porb(i))=XP(i+1)
!          end do

            write(1,*) "  parameter  :  (1 is done) "
            call print_mat(1,1,XP(1),1)
            call flush(1)
            write(1,*) " orbital rot :  "
            call print_mat(1,nijO,XP(2:nijO+1),1)
            call flush(1)
            write(1,*) "  MPS_ci rot :  "
            call print_mat(1,nijC,XP(nijO+2:nij+1),1)
            call flush(1)
!          stop

            !call print_mat(ndim1,ndim1,mat2%A,6)
!            call print_mat(ndim1,ndim1,mat2%B,6)
!         update the A and B matrix        

!          do iC=1,nC
!            mat2%A=mat2%A+2.0d0*X(iC+ndim2)*Ci(iC)%A
!            mat2%B=mat2%B+2.0d0*X(iC+ndim2)*Ci(iC)%A
!            write(6,*)iC,X(iC+ndim2)
!          end do 

            !call print_mat(ndim1,ndim1,mat2%A,6)
!            call print_mat(ndim1,ndim1,mat2%B,6)

! ===================================================================
!             recover all orbital rotation
! ===================================================================
              ij=0;
              d0=0.0d0; d2=0.0d0;
              d0_sum=0.0d0; d2_sum=0.0d0
              do i=1,nijO  ! Notice here should only orbital index
                i1=valid(i,1)
                j1=valid(i,2)
! recover the anti-symmetric rotations
                mat2%deltR(i1,j1)=  XP(i+1)/2.0d0
                mat2%deltR(j1,i1)= -XP(i+1)/2.0d0
                write(1,*)"deltaR in",iloop,"micro-iter",&
                          i1,j1,mat2%deltR(i1,j1)
                if(dabs(mat2%deltR(i1,j1)).gt.0.10)then
                  write(2,*)"notice : large rotation between",i1,j1,&
                            " : ",mat2%deltR(i1,j1)
                end if
                d0=d0+dabs(XP(i+1))
                d2=d2+(XP(i+1))**2
              end do
              d0_sum = d0_sum + d0
              d2_sum = d2_sum + d2

! recover the degenrate rotations
              is0=0
              do isym=1,orb%nsub
                js0=0
                do jsym=1,orb%nsub
                  if(redu%irreps(isym).eq.redu%irreps(jsym)&
                    .and.isym.lt.jsym)then
                    do i=1,orb%total(isym)
                      do j=1,orb%total(isym)
                        ii=i+is0; jj=j+is0
                        kk=i+js0; ll=j+js0
                ! avoid possible over-shooting by 1/2 ()
                        mat2%deltR(ii,jj)=&
                        mat2%deltR(ii,jj)/2.0d0
                        mat2%deltR(kk,ll)=& 
                        mat2%deltR(ii,jj)*redu%sign_mat(kk,ll)
                        write(1,*)ii,jj,mat2%deltR(ii,jj),&
                                  kk,ll,mat2%deltR(kk,ll)
                      end do
                    end do
                  end if
                  js0=js0+orb%total(jsym)
                end do
                is0=is0+orb%total(isym)
              end do
! ===================================================================

         ! 1) Orb modification
            allocate(T1(norb,norb)); T1=0.0d0
            allocate(T2(norb,norb)); T2=0.0d0
            call CAL_T(norb,mat2%deltR,1.0d0,150,T1)
            call T_TO_U(norb,T1,T2)
            call MXM(norb,mat2%U,T2,T1)
            mat2%U=T1
            mat2%T=mat2%U
            do i=1,norb
              do j=1,norb
                if(dabs(mat2%U(i,j)).lt.1.0d-16)then
                  mat2%U(i,j)=0.0d0
                end if
              end do
              mat2%T(i,i)=mat2%U(i,i)-1.0d0
            end do
            deallocate(T1,T2)

!           mat2%R=mat2%deltR

         ! 2) MPSci modification
            if(iloop.gt.3)then
              call renormalize_deltaCi&
                (nC,CI(1:nC)%dv,XP(nijO+2:nij+1),CI(1:nC)%dvp) 
            end if 

!          ij=0; d0=0.0d0
!          do i=1,ndim1
!            do j=1,ndim1
!               ij=ij+1
!               mat2%R(i,j)=X(ij)
!               d0=d0+dabs(X(ij))
!            end do
!          end do
!          Rabs(icycle)=d0

            call print_mat(1,nC,CI(:)%dv,1)
            call flush(1)           
 
            if(iloop.gt.iloopmax)then
              write(2,*)"reach the maximum loop",iloopmax
              exit
            end if
            write(2,*)"Finish the ",iloop," iteration"
            call flush(2)
            if(grad_digit.lt.redundant)then
              write(2,*)" negligible residual "
              exit
            end if
            if(dabs(XP(1)).gt.0.99999999d0)then
              write(2,*)" XP(1) flag close to 1.0d0 "
              exit
            end if
            write(1,*)XP(1) 

            if(ith_HESS.ne.2)then
              UMAT_FINAL=mat2%U
            end if

            deallocate(HPA)
            deallocate(XP)
            deallocate(D)

          else
  
            allocate(D(nij));       D=0.0d0
            allocate(XP(nij));     XP=0.0d0
            allocate(T1(nij,nij)); T1=0.0d0
            allocate(T2(nij,nij)); T2=0.0d0

            T1=HP
            call eigtql2(nij,T1,D)
  
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
              write(6,*)"Value of X",i,XP(i) 
              X(porb(i))=XP(i)
            end do
            Rabs(icycle)=d0
 
            deallocate(XP)
            deallocate(D)
            deallocate(T1,T2)
 
          end if ! advanced solver or not


        end do ! loop for the micro-iterations

        d0_sum0=0.0d0
        d2_sum0=0.0d0
        do i=1,norb
          do j=1,norb
            if(i.ne.j)then
              d0_sum0=d2_sum0+dabs(UMAT_FINAL(i,j))
              d2_sum0=d2_sum0+UMAT_FINAL(i,j)**2
            end if
          end do
        end do

        write(*,2434)E0_micro+fronre0, dsqrt(d2_sum0)
2434    format(4X,f14.8,4X,f8.4)

        Rabs(icycle)=d0_sum0
        write(2,*)icycle,d0_sum0,d2_sum0

        deallocate(H)
        deallocate(HP)
        deallocate(YP)
        deallocate(X,Y)
        deallocate(porb)
        deallocate(valid)

!        write(6,*)"Finishing the coupling_solver"
!        stop

        
      ! recover the checkpoint file, in case it is changed in micro-iters
        do iroot=0,dmrg_nstates-1

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if

          open(unit=200,file="rdm_ckp_operate_mi")        
            string0=""
            string0="checkpoint_state."//trim(ctmp2)          ! string0  (ori_)  
            string4=""
            string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_) 

            string2="rm -rf "//trim(string0)
            write(200,*)trim(string2)
            string2="cp -r "//trim(string4)//" "//trim(string0)
            write(200,*)trim(string2)
          close(200)          
          call system("chmod +x rdm_ckp_operate_mi")
          call system("./rdm_ckp_operate_mi")

        end do


      end Subroutine coupling_solver


      Subroutine renormalize_deltaCi(nC,C0,C1,deltaC)

        integer::nC
        double precision::C0(nC),C1(nC),deltaC(nC)
        double precision dL

!        write(1,*) "C0", C0 
!        write(1,*) "C1", C1 
!        write(1,*) "dC", deltaC 

        deltaC=C0+C1 
        !deltaC=C0-C1 
        dL=0.0d0 
        do i=1,nC
          dL=dL+deltaC(i)**2
        end do 
        dL=dsqrt(dL) 
 
        deltaC=deltaC/dL

        C0=deltaC

      end Subroutine renormalize_deltaCi


      Subroutine coupling_solver3(icycle)

        use global_control
        use coupling_terms
        use matrix

        integer::icycle
 

      end Subroutine coupling_solver3

! 2) treat the MPS change as the qusi-second order response 
!          In this case, couplings ralated to MPS are 
!     added to orb-orb Hessian using the perturbation
!     way up to first order;
!     (Actually, many people treat it as second order)
!     Advantage is that only the orbitals rotation 
!     need to be solved.

      Subroutine coupling_solver2(icycle)

        use global_control
        use coupling_terms
        use matrix

        integer::icycle
! This solver should be very similar to the Solver_Hessian
        integer,allocatable::porb(:)
        double precision d0,redundant
        double precision,allocatable::D(:),X(:),Y(:),XP(:),YP(:)
        double precision,allocatable::H(:,:),HCCinv(:,:)
        double precision,allocatable::H2(:,:),H2vec(:,:),H2val(:)
        double precision,allocatable::HP(:,:),HPA(:,:)
! Temporary value
        double precision dv
        double precision lamda,levelshift
! Temporary matrix
        double precision,allocatable::T1(:,:),T2(:,:),T3(:,:),T4(:,:)
        double precision,allocatable::Hcorr(:,:),HcorrP(:,:)

        integer offset 

        lamda=1.00d0 ! "1" means the standard AH 
        levelshift=0.0 

! Assuming the whole Hessian is needed by ndim1 = norb
        ndim1=norb
        ndim2=ndim1**2

! initialize of Hessian matrix that will be used in this solver
        allocate(X(ndim2)); X=0.0d0
        allocate(Y(ndim2)); Y=0.0d0
        allocate(H(ndim2,ndim2)); H=0.0d0
        allocate(H2(nC-1,nC-1));  H2=0.0d0
        allocate(H2val(nC-1));    H2val=0.0d0
        allocate(H2vec(nC-1,nC-1));  H2vec=0.0d0
        allocate(HCCinv(nC-1,nC-1)); HCCinv=0.0d0
       


 
! 1) The MKL way to inverse        
! Well, it is not stable due to the diagional is not dominant 
! (need consider the k!=0 case)
!        goto 2233 ! Not use this way due to 1/0 exist
        if(.true.)then
          H2=HCC(2:nC,2:nC) 
          HCCinv=HCC(2:nC,2:nC) 
!          do i=1,nC
!            HCCinv(i,i)=HCCinv(i,i)
!          end do
          call EIGTQL2(nC-1,H2vec,H2val)

!          print *,"The HCC inv (MKL) before inv : " 
          !call print_mat(nC-1,nC-1,HCCinv,6)       
          call inv(nC-1,HCCinv)
          !print *,"The HCC inv after inv : " 
          !call print_mat(nC-1,nC-1,HCCinv,6)

          !stop
! 2) Other way as D^-1_IJ=\sum_k c^k_I (E_k-E_0) c^k_J
! Which means using the existing E-diag and E_vec 
!        2233 continue ! Use this way instead       
        else 
          do i=1,nC
            do j=1,nC
              do k=1,nC ! means avoid the k=0 case
              if(k.eq.1)then
                dv=0.0d0
              else             
                dv=1.0d0/(2.0*(Hval(k)-RDM_ENERGY))
               !dv=1.0d0/(Hval(k))
              end if
!                write(*,*)i,k,"Hvec(k,i)",Hvec(k,i)
                HCCinv(i,j)=HCCinv(i,j)+Hvec(k,i)*dv*Hvec(k,j)
              end do
            end do
          end do         
        end if 

        print *,"--> Check the HCCinv * Hcc (inv = 1) : " 
        !allocate(T1(nC,nC)); T1=0.0d0
        allocate(T1(nC-1,nC-1)); T1=0.0d0
        !call print_mat(nC,nC,HCC,6) 
        !call print_mat(nC,nC,HCCinv,6) 
        !call MXM(nC-1,HCC(2:nC,2:nC),HCCinv(2:nC,2:nC),T1)
        call MXM(nC-1,HCC(2:nC,2:nC),HCCinv,T1)
        call print_mat(nC-1,nC-1,T1,6) 
        deallocate(T1)
        !stop 

! calculate the orb-orb Hessian
          mat2%H=0.0d0 
          call Hessian(ndim1,mat2%A,mat2%G,mat2%H,mat2%Hdiag)        
!        call PRINT_GAT(ndim1,ndim1,ndim1,ndim1,mat2%H,702)

        allocate(Hcorr(ndim2,ndim2)) ! The Ci_orb correction to Hessian (orb)
        Hcorr=0.0d0

        if(.true.)then
          allocate(T1(nC-1,ndim2))   ; T1=0.0d0
          allocate(T2(nC-1,ndim2))   ; T2=0.0d0
          allocate(T3(ndim2,ndim2)); T3=0.0d0 
          allocate(T4(ndim2,nC-1))   ; T4=0.0d0
! Get the orb-MPS coupling matrix
          do iC=1,nC-1  ! except the dominant one
            ij=0
            do i=1,ndim1
              do j=1,ndim1 
                ij=ij+1
                T1(iC,ij)=HCR(iC+1,i,j) 
                T2(iC,ij)=HCR(iC+1,i,j) 
              end do
            end do
          end do
          !T2=T1
        !call print_tat(nC,ndim1,ndim1,HCR,703)
        !call print_mat(ndim2,nC,T1,704)
        !call print_mat(nC,ndim2,T2,705)
        !call print_mat(nC,nC,HCC,706)
        !call print_mat(nC,nC,HCCinv,707)
!        stop

! Coupling the MPS change into orbital in a perturbation way
! This also work, but currently use the old way for reading 
!  (Coupled at the first order with MPS)
         call MXMG(ndim2,nC-1,nC-1,T1,HCCinv,T4,"TN")  ! H o-m * D^-1 m-m = T
!        call print_mat(ndim2,nC,T4,708)
!        stop
         call MXMG(ndim2,ndim2,nC-1,T4,T2,T3,"NN")   ! T o-m * H m-o = New Hm-m
         !print *,"--> The deltaH : " 
         !call print_mat(ndim2,ndim2,T3,6)

          ij=0
          do i=1,ndim1
            do j=1,ndim1
              ij=ij+1
              kl=0
              do k=1,ndim1
                do l=1,ndim1
                  kl=kl+1
! Actually, this is "plus" rather than minus, but works (?, check later)
                  H(ij,kl)=mat2%H(i,j,k,l)
                  H(ij,kl)=H(ij,kl)-T3(ij,kl)
                end do
              end do
            end do
          end do
 
          Hcorr=T3
          deallocate(T1,T2,T3,T4)

        else 

          allocate(deltaH(ndim1,ndim1,ndim1,ndim1))
          deltaH=0.0d0

          do i=1,ndim1
            do j=1,ndim1
              do k=1,ndim1
                do l=1,ndim1
                  dv=0.0d0 
                  do iC=1,nC-1
                    do jC=1,nC-1
                      dv=dv+HCR(iC+1,i,j)*HCCinv(iC,jC)*HCR(jC+1,k,l)
                    end do
                  end do
                  deltaH(i,j,k,l)=dv
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
! Actually, this is "plus" work rather than minus (??????)
                   H(ij,kl)=mat2%H(i,j,k,l)-deltaH(i,j,k,l)
                  !H(ij,kl)=mat2%H(i,j,k,l)+deltaH(i,j,k,l)
                  !H(ij,kl)=mat2%H(i,j,k,l)
                  
                  Hcorr(ij,kl)=deltaH(i,j,k,l)
                end do
              end do
            end do
          end do
          !print *,"The CI_orb correction to Hessian(orb) part : "
          !call PRINT_MAT(ndim2,ndim2,T1,6)

!          deallocate(T4)
          deallocate(deltaH)

        end if  

!          call PRINT_GAT(ndim1,ndim1,ndim1,ndim1,mat2%H,702)
!          call PRINT_GAT(ndim1,ndim1,ndim1,ndim1,deltaH,715)

!        call PRINT_MAT(ndim2,ndim2,H,712)
!        call PRINT_MAT(ndim2,ndim2,T3,713)
!        stop  

!        call PRINT_GAT(ndim1,ndim1,ndim1,ndim1,H,712)
!        call PRINT_GAT(ndim1,ndim1,ndim1,ndim1,T3,713)

!        stop

        !call PRINT_MAT(ndim1,ndim1,mat2%A,6)
        !stop 

        ! mcscf Fock/gradient matrix 
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
               print *, "orb",i,j,"gradient",Y(ij)
             end if
!             H(ij,ij)=100+ij ! when testing
          end do
        end do

! only consider the non-redundant rotation (nij)
        allocate(porb(nij));   porb=0
        allocate(YP(nij));     YP=0.0d0
        allocate(HP(nij,nij)); HP=0.0d0
        allocate(HcorrP(nij,nij)); HcorrP=0.0d0

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
!            H(ij,ij)=100+ij ! when testing
          end do
        end do

        do i=1,nij
          YP(i)=Y(porb(i))
          do j=1,nij
            HP(i,j)=H(porb(i),porb(j))
            HcorrP(i,j)=Hcorr(porb(i),porb(j))
          end do
        end do

!        print *,"The effective Hessian" 
!        call print_mat(nij,nij,HP,6)

         print *,"The Hessian correction from (MPS)CI-orb part :" 
         call print_mat(nij,nij,HcorrP,6)
        !stop

!        do i=1,nij
!          do j=1,nij
!            write(321,*)i,j,HP(i,j) 
!          end do
!        end do
!        call flush(6)

        if(.true.)then  ! Here use modified Davidson (AH)

          YP=-1.0d0*YP
          allocate(HPA(nij+1,nij+1))
          HPA=0.0d0

! The new augmented Hessian matrix {except HPA(1,1)}
          HPA(2:nij+1,1)=YP
          HPA(1,2:nij+1)=YP
          HPA(2:nij+1,2:nij+1)=HP
! And the augmented x-vector
          allocate(XP(nij+1))
          allocate(D(nij+1))
          XP=0.0d0
           D=0.0d0
          call Davidson(nij+1,HPA,XP,D,lamda,.true.,100,1.0d-12)

          write(*,*)"AH_inital"
          write(*,*) XP

          if(XP(1).lt.0)then
            XP=-1.0d0*XP
          end if

          if(.true.)then
            do
              call MODIFIED_GRAM_SCHMIDT(nij+1,1,XP)
              if(XP(1).lt.0)then
                XP=-1.0d0*XP
              end if
              if(dabs(XP(1)).lt.levelshift)then
                XP(1)=XP(1)*1.25d0
                XP(2:nij+1)=XP(2:nij+1)/1.25d0
              else
               !levelshift=dabs(XP(1)) ! reset levelshift
                exit
              end if
            end do
          end if
          write(*,*)"AH_restricted"
          write(*,*) XP

!          stop
          do i=1,nij
            X(porb(i))=XP(i+1)
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

          deallocate(HPA)
          deallocate(XP)
          deallocate(D)

        else

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
            if(dabs(D(i)).lt.1.0e-12)then
              D(i)=0.0d0
            end if
            write(322,*)"hessian eigval in PH solver",i,D(i) 
          end do

          XP=0.0d0
          do i=1,nij
            do k=1,nij
              do l=1,nij
                if(D(l).gt.1.0e-9)then
                  dv=T1(i,l)*T1(k,l)*YP(k)/D(l) 
!                  if(dabs(dv).gt.0.1)dv=0.1
                  XP(i)=XP(i)+dv
                else if(D(l).lt.-1.0e-9)then
                  dv=T1(i,l)*T1(k,l)*YP(k)/D(l)
!                  if(dabs(dv).gt.0.1)dv=0.1
                  XP(i)=XP(i)-dv
                end if
              end do
            end do
            write(323,*)"Value of X",i,XP(i) 
            X(porb(i))=XP(i)
          end do

          deallocate(D,XP,T1,T2)

        end if 

        ij=0; d0=0.0d0
        do i=1,ndim1
          do j=1,ndim1
             ij=ij+1
             mat2%R(i,j)=X(ij)
             d0=d0+dabs(X(ij))
          end do
        end do
        Rabs(icycle)=d0       
        call print_mat(ndim1,ndim1,mat2%R,6) 
        print *,"Rabs(icycle)",Rabs(icycle)
!        stop 

!       Now, deal with the remaining rotations

        call CAL_T(norb,mat2%R,d3,150,mat2%T)
!        print *," === comment the CAL_T(above)  === "
        call T_TO_U(norb,mat2%T,mat2%U)
        call GRAM_SCHMIDT(norb,mat2%U)

!        print *," === The transformed matrix U === "
!        call print_mat(norb,norb,mat2%U,6)
!        print *," ================================ "

        write(*,2435,advance='no')"    ",Rabs(icycle)
        write(*,2436)"   ",method(icycle)
2435    format(A4,f8.4)
2436    format(A3,A7)

!        deallocate(T1,T2)
        deallocate(X)
        deallocate(Y,YP)
        deallocate(Hcorr)
        deallocate(H,HP)
        deallocate(H2,H2vec,H2val)

        print *,"finishing coupling_solver2"        

      end Subroutine coupling_solver2

