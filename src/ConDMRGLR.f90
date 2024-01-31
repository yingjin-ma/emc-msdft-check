      Subroutine ConDMRGLR

        use global_control
        use matrix

        double precision,allocatable::TM1(:,:),TM2(:,:)
        double precision,allocatable::TM3(:,:),TM4(:,:)
        double precision,allocatable::TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
        double precision,allocatable::GM3(:,:,:,:),GM4(:,:,:,:)
        double precision,allocatable::TL1(:),TL2(:),TL3(:)

        character    ctmp1
        character*2  ctmp2
        character*192 string0,string1,string2,string3,string4

        character*72 oneRDMfile,twoRDMfile

        double precision dv
        double precision grad_digit
        double precision redundant

        double precision,allocatable::Y(:),G(:)
        double precision,allocatable::V1(:),V2(:),V3(:)

        ! specific state & valid rotations
        integer rlx,offset
        integer,allocatable::valid(:,:)
        integer,allocatable::valid_add(:,:)
        logical Lvalid

        ! orthogonal SA space
        integer,allocatable::conf(:)

        ! effective orbital index & residual
        integer,allocatable::porb(:),porb_redu(:,:),porb_tmp(:)
        double precision,allocatable::YP(:),YP0(:),YP00(:)
        double precision,allocatable::LP(:),LPx(:,:),LP0(:)
        double precision,allocatable::Dia(:,:),Dsa(:,:)

        double precision,allocatable::H(:,:),HP(:,:),HP2(:,:)
        double precision,allocatable::SAvec(:),HPM(:,:)

        ! Parameters that used to control PCG/CG process
        ! (It need to be separately considered, caused by closed-shell approximation)
        double precision rAlphaORB,rdeltaORB,resMPS,rsigmaORB
        double precision rAlphaMPS,rdeltaMPS,resORB,rsigmaMPS
        double precision rAlpha,rdelta,res,rsigma,rbeta
        double precision,allocatable::sigma(:),xvec(:),rvec(:)
        double precision dv1,dv2,dv3

        ! inactive-Fock and active-Fock
        double precision FockI(norb,norb)
        double precision FockA(norb,norb)

        ! number of non-redundant rotation for each orb (with correction from closed shell)
        integer nprei(norb)

!        ! flag for checking the phase
!        logical ifcorrectedRHS

        nprei=0
        ncls=nocc-nact
        nvir=norb-nocc

! 1-RDMs
        if(.not.allocated(mat2%d))then
           allocate(mat2%d(nact,nact))           ; mat2%d=0.0d0
        else
           mat2%d=0.0d0
        end if
        allocate(TM1(nact,nact))                 ;    TM1=0.0d0

! 2-RDMs
        if(.not.allocated(mat2%p))then
          allocate(mat2%p(nact,nact,nact,nact))  ; mat2%p=0.0d0
        else
           mat2%p=0.0d0
        end if
        allocate(GM1(nact,nact,nact,nact))       ;    GM1=0.0d0

        do iroot=0,dmrg_nstates-1

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if

          oneRDMfile="oneRDM."//trim(ctmp2)//"."//trim(ctmp2)
          twoRDMfile="twoRDM."//trim(ctmp2)//"."//trim(ctmp2)

          mat2%d=0.0d0
          open(unit=110,file=trim(oneRDMfile))
            read(110,*)ij
            do i=1,nact
              do j=1,nact
                read(110,*)ij,ji,mat2%d(iorder(i),iorder(j))
              end do
            end do
          close(110)
          write(2,*)"The DMRG calculated 1-RDM for state",iroot
          call print_mat(nact,nact,mat2%d,2)
          TM1=TM1+mat2%d*dmrg_weight(iroot+1)

          mat2%p=0.0d0
          open(unit=110,file=trim(twoRDMfile))
            READ(110,*)ijkl
            ijkl=0  ! Just a counter
            do
              read(110,*,iostat=ierr)ij,jk,kl,li,dv
              if(ierr.ne.0)exit
              ijkl=ijkl+1
              mat2%P(iorder(ij+1),&
                     iorder(jk+1),&
                     iorder(kl+1),&
                     iorder(li+1))=dv
              mat2%P(iorder(kl+1),&
                     iorder(li+1),&
                     iorder(ij+1),&
                     iorder(jk+1))=dv
              mat2%P(iorder(jk+1),&
                     iorder(ij+1),&
                     iorder(li+1),&
                     iorder(kl+1))=dv
              mat2%P(iorder(li+1),&
                     iorder(kl+1),&
                     iorder(jk+1),&
                     iorder(ij+1))=dv
!                    write(*,*)ij+1,jk+1,kl+1,li+1,rdm2P(ij+1,jk+1,kl+1,li+1)
            end do
          close(110)
          write(34,*)"The DMRG calculated 2-RDM for state",iroot
          call print_gat(nact,nact,nact,nact,mat2%p,34)
          GM1=GM1+mat2%p*dmrg_weight(iroot+1)

        end do

        ! return the weighted RDMs
        mat2%d=TM1
        mat2%p=GM1

        mat2%p=mat2%p*2.0d0
!       call symmetrize()
        deallocate(TM1,GM1)

        ! Closed shell update
        ! call closed_shell_update()  !CCCC

        write(2,*)"Re-read-in all RDMs"
        call flush(2)

        if(.NOT.ALLOCATED(mat2%A))then
          allocate(mat2%A(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%Hdiag))then
          allocate(mat2%Hdiag(norb,norb))
        end if
        if(.NOT.ALLOCATED(mat2%G))then
          allocate(mat2%G(norb,norb,norb,norb)) ! still use the old form
        end if
        if(.NOT.ALLOCATED(mat2%Y))then
          allocate(mat2%Y(norb,norb,norb,norb)) ! use the same form as G
        end if
        if(.NOT.ALLOCATED(mat2%Y1))then
          allocate(mat2%Y1(norb,norb)) ! use the same form as Y1
        end if
        if(.NOT.ALLOCATED(mat2%Y2))then
          allocate(mat2%Y2(norb,norb)) ! use the same form as Y2
        end if
        if(.NOT.ALLOCATED(mat2%Dh))then
          allocate(mat2%Dh(norb,norb,norb,norb)) ! use the same form as G
        end if
        if(.NOT.ALLOCATED(mat2%H))then
          allocate(mat2%H(norb,norb,norb,norb)) !
        end if
        if(.NOT.ALLOCATED(mat2%T))then
          allocate(mat2%T(norb,norb))
        end if

        write(2,*)"after allocate mat2"
        call flush(2)

        write(2,*)" mat2%D and mat2%P : mat2%D "
        call print_mat(nact,nact,mat2%D,2)
!        write(2,*)" mat2%D and mat2%P : mat2%P "
!        call print_gat(nact,nact,nact,nact,mat2%P,2)

! to Leon :!     Notice that the same redundancy check should be emplyed for LR part

        ! construct the orbital gradient for the specific state
        rlx=dmrg_rlxstate
        redundant=1.0d-9
        grad_digit=0.0d0
        ! Check for the no-redundant orbiatl-rotations
        allocate(Y(norb**2)); Y=0.0d0
        allocate(G(norb**2)); G=0.0d0
        ij=0
        nij=0
        nij_add=0

!1) option-1 : If use o-v, i-p to check the redundant
!        open(unit=110,file="orb_grad.txt")
!          ij=0
!          do i=1,ncls
!            do j=ncls+1,norb
!              ij=ij+1
!              write(110,*)i,j,MPS(rlx)%G(i,j),"i-p",ij
!            end do
!          end do
!          do i=ncls+1,nocc
!            do j=nocc+1,norb
!              ij=ij+1
!              write(110,*)i,j,MPS(rlx)%G(i,j),"o-v",ij
!            end do
!          end do
!        close(110)

!        stop

!2) option-2 : If use predefined orb-gradients to check the redundant

        do i=1,norb
          do j=i+1,norb ! only half is needed
            ij=ij+1
            Y(ij)=MPS(rlx)%G(i,j)  ! orbital gradients for the specific state
            if(dabs(Y(ij)).lt.redundant)then
            else
!                 write(1,*)"in else"
              IF(LMCSCF)then
                write(*,*)"LR can only handle CASSCF-type rotations "
                stop
                ! Check for degeneracy irreps
                if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                  write(1,*)"LR gradient check",i,j,Y(ij),ij
                  nij=nij+1
                  if(dabs(Y(ij)).gt.grad_digit)then
                    grad_digit=dabs(Y(ij))
                  end if
                end if
                ! end Check for degeneracy irreps
              else
                ! no closed shell
                if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    ! Check for degeneracy irreps
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      write(1,*)"LR gradient check",i,j,Y(ij),ij
                      nij=nij+1
                      if(dabs(Y(ij)).gt.grad_digit)then
                        grad_digit=dabs(Y(ij))
                      end if
                    end if
                    ! end Check for degeneracy irreps
                  end if
                ! i is closed shell
                else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"LR gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(G(ij)).gt.grad_digit)then
                      grad_digit=dabs(G(ij))
                    end if
                    ! Additional rotation caused by closed shell
                    if(redu%orbirp(j).lt.0)then
                      nij_add=nij_add+1
                    end if
                  end if
                ! j is closed shell
                else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"LR gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(G(ij)).gt.grad_digit)then
                      grad_digit=dabs(G(ij))
                    end if
                  end if
                end if ! close shell
              end if ! MCSCF or CASSCF
            end if
          end do
        end do

        write(1,*)"Checked for the no-redundant orbiatl-rotations 1/2"
        call flush(1)

        do i=1,norb
          do j=1,i ! only half is needed
            ij=ij+1
            Y(ij)=MPS(rlx)%G(i,j)
            if(dabs(Y(ij)).lt.redundant)then
            else
!                 write(1,*)"in else"
              IF(LMCSCF)then
                write(*,*)&
                "Currently, LR can only handle CASSCF-type rotations "
                stop
                ! Check for degeneracy irreps
                if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                  write(1,*)"LR gradient check",i,j,Y(ij),ij
                  nij=nij+1
                  if(dabs(Y(ij)).gt.grad_digit)then
                    grad_digit=dabs(Y(ij))
                  end if
                end if
                ! end Check for degeneracy irreps
              else
                ! no closed shell
                if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    ! Check for degeneracy irreps
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      write(1,*)"LR gradient check",i,j,Y(ij),ij
                      nij=nij+1
                      if(dabs(Y(ij)).gt.grad_digit)then
                        grad_digit=dabs(Y(ij))
                      end if
                    end if
                    ! end Check for degeneracy irreps
                  end if
                ! i is closed shell
                else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"LR gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(G(ij)).gt.grad_digit)then
                      grad_digit=dabs(G(ij))
                    end if
                  end if
                ! j is closed shell
                else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    write(1,*)"LR gradient check",i,j,Y(ij),ij
                    nij=nij+1
                    if(dabs(G(ij)).gt.grad_digit)then
                      grad_digit=dabs(G(ij))
                    end if
                    ! Additional rotation caused by closed shell
                    if(redu%orbirp(i).lt.0)then
                      nij_add=nij_add+1
                    end if
                  end if
                end if ! close shell
              end iF ! MCSCF-type or not
            end if
          end do
        end do

        write(1,*)"Checked for the no-redundant orbiatl-rotations 2/2"
        write(1,*)"nij and nij_add",nij,nij_add
        call flush(1)

        nij0=nij
! re-run it to get the valid rotations index
        allocate(valid(nij,2)); valid=0
        allocate(valid_add(nij_add,2)); valid_add=0
        ij=0
        kl=0
        kl_add=0
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
                ! no closed shell
                if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      kl=kl+1
                      valid(kl,1)=i
                      valid(kl,2)=j
                    end if
                  end if
                ! i is closed shell | i-a,i-v
                else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                    if(redu%orbirp(j).lt.0)then
                      kl_add=kl_add+1
                      !valid_add(kl_add,1)=i
                      !valid_add(kl_add,2)=j
                      ! Treat separately for G_ai and G_ia
                      valid_add(kl_add,1)=j
                      valid_add(kl_add,2)=i
                    end if
                  end if
                ! j is closed shell | j-a,j-v
                else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                  end if
                end if !
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
                ! no closed shell
                if(IABS(redu%orbirp(i)).lt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(redu%orbirp(i)*redu%orbirp(j).lt.0)then
                    if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                      kl=kl+1
                      valid(kl,1)=i
                      valid(kl,2)=j
                    end if
                  end if
                ! i is closed shell | i-a,i-v
                else if(IABS(redu%orbirp(i)).gt.10.and.IABS(redu%orbirp(j)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                  end if
                ! j is closed shell | j-a,j-v
                else if(IABS(redu%orbirp(j)).gt.10.and.IABS(redu%orbirp(i)).lt.10)then
                  if(dabs(redu%sign_mat(i,j)).lt.1.0d-8)then
                    kl=kl+1
                    valid(kl,1)=i
                    valid(kl,2)=j
                    if(redu%orbirp(i).lt.0)then
                      kl_add=kl_add+1
                      !valid_add(kl_add,1)=i
                      !valid_add(kl_add,2)=j
                      ! Treat separately for G_ai and G_ia
                      valid_add(kl_add,1)=j
                      valid_add(kl_add,2)=i
                    end if
                  end if
                end if !
              end if
            end if
          end do
        end do

        write(1,*)"valid rotations",nij
        call flush(1)
        do i=1,nij
          write(1,"(I3)",advance='no')valid(i,1)
        end do
        write(1,*)" "
        do i=1,nij
          write(1,"(I3)",advance='no')valid(i,2)
        end do
        call flush(1)

        write(1,*)
        write(1,*)"Additional valid rotations caused by closed shell",&
                   nij_add
        do i=1,nij_add
          write(1,"(I3)",advance='no')valid_add(i,1)
        end do
        write(1,*)" "
        do i=1,nij_add
          write(1,"(I3)",advance='no')valid_add(i,2)
        end do

        write(1,*)
        write(1,*)"Total non-redundant rotations :", (nij+nij_add)/2

        if(nij.eq.0)then
          write(*,*)"No orbital gradients large than ", redundant
          write(*,*)"Lagrange optimization is not need"
          stop
        end if
! to Leon :!     Always need the closed shell approximation for MCLR part

! Closed shell extention for SA
          allocate(TM1(nocc,nocc))
          allocate(GM1(nocc,nocc,nocc,nocc))
          !call print_mat(nact,nact,MPS(iroot+1)%d,6)
          !call print_gat(nact,nact,nact,nact,MPS(iroot+1)%p,6)
          call closed_shell_update&
          (mat2%d,mat2%P,TM1,GM1,1.0d0)  ! CCCC
          deallocate(mat2%d,mat2%P)
          allocate(mat2%d(nocc,nocc))
          allocate(mat2%P(nocc,nocc,nocc,nocc))
          mat2%d=TM1
          mat2%p=GM1
!          call print_mat(nocc,nocc,mat2%d,6)
!          call print_gat(nocc,nocc,nocc,nocc,mat2%p,6)
          deallocate(TM1,GM1)

! SA energy check
        call energy_gen(orb%nsub,orb%act,orb%total,&
                  nact,norb,T,U,mat2%D,mat2%P,SA_RDMenergy)
        write(2,*)"Checked the SA_RDMenergy in LR - AG",SA_RDMenergy
        write(2,*)"            SA_RDMenergy (-NRE)    ",SA_RDMenergy+FRONRE

        dv=0.0d0
        do iroot=0,dmrg_nstates-1
          dv=dv+MPS(iroot+1)%RDM_energy*dmrg_weight(iroot+1)
        end do
        write(2,*)"Checked the reformed RDMenergy",dv
        call flush(2)
        if(dabs(dv-SA_RDMenergy).gt.thrs%e)then
          write(*,*)"Warnning: SA_RDMenergy .ne. reformed RDMenergy"
          write(*,*)"          first consider to check RDM-deri"
        else
          write(2,*)"Matched for the reformed RDMenergy"
        end if

! Condtruct the A,G,H matrix
!       Construct the A matrix base on T/U integrals
        call fock_gen(orb%nsub,orb%act,orb%total,&
                  nact,norb,T,U,mat2%D,mat2%P,mat2%A)

        write(2,*)"The SA A-mat :"
        call print_mat(norb,norb,Mat2%A,2)  ! checked
        call flush(2)
!        mat2%A=0.0d0

! Construct the G matrix
        call Gmat_gen(orb%nsub,orb%act,orb%total,&
             nact,norb,T,U,mat2%D,mat2%P,mat2%G,orb%grouptable)

!        write(2,*)"second round"
!        call flush(2)
!        call print_GAT(norb,norb,norb,norb,mat2%G,9)

! to Leon :
!       construct the SA orbital-orbital Hessian
!       Many amendments were implemented in order to generate the same results as Molcas
!       differences emerge if closed shell elements (i.e. Hia,ia !== Hai,ai in Molcas)

        iHessian=1

        if(iHessian.eq.0)then
          call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)

          call print_gat(norb,norb,norb,norb,mat2%H,321)
          call flush(321)

        else if(iHessian.eq.1)then
          ! Same way as in "Molecular electronic-structure theory"

          write(2,*)"Dh_gen start"
          call flush(2)

          write(2,*)"orb%nsub" ,  orb%nsub
          call flush(2)
          write(2,*)"orb%act"  ,  orb%act
          call flush(2)
          write(2,*)"orb%total",  orb%total
          call flush(2)
          write(2,*)"nact"     ,     nact
          call flush(2)
          write(2,*)"norb"     ,     norb
          call flush(2)
!          write(2,*)"T"        ,        T
!          call flush(2)
!          write(2,*)"mat2%D"   ,   mat2%D
!          call flush(2)

          call Dh_gen(orb%nsub,orb%act,orb%total,&
                     nact,norb,T,mat2%D,mat2%Dh,orb%grouptable)

          write(2,*)"Dh_gen finish"
          call flush(2)

          call Y_gen(orb%nsub,orb%act,orb%total,&
                     nact,norb,U,mat2%P,mat2%Y,orb%grouptable)

          write(2,*)"Y_gen finish"
          call flush(2)

!          call HessianHF(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)
          call Hessian3(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)

          write(1,*)"Hessian3 finish, print the Hess(ii,jj)"
          call print_mat(norb,norb,mat2%Hdiag,1)
          call flush(1)

!          call Hessian2(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)

!          call print_gat(norb,norb,norb,norb,mat2%H,322)
!          call flush(322)

! to Leon :! Further Correct the diagional Hess_ia,ia
          call CS_dim_fro2cls()
          allocate(Dia(ncls,nact)); Dia=0.0d0
          allocate(Dsa(nvir,nact)); Dsa=0.0d0

!   inactive Fock matrix in MO basis
          FockI=0.0d0
          call fock_I(orb%nsub,orb%closed,orb%total,norb,T,U,FockI)
!          write(2,*)" ==== Fock_I ===="
!          call print_mat1(norb,norb,FockI,2)

!     active Fock matrix in MO basis
          FockA=0.0d0
          call fock_A(orb%nsub,orb%occ,orb%closed,orb%total,&
                      nocc,norb,U,mat2%D,FockA)
!          write(2,*)" ==== Fock_A ===="
!          call print_mat1(norb,norb,FockA,2)

!          write(6,*)"Fock MC"
!          call print_mat(norb,norb,mat2%A/2.0d0,6)

! Hessian for 2index ii, diagional elements are used
          call precaii(orb%nsub,orb%occ,orb%closed,orb%total,&
                       norb,nocc,ncls,mat2%D,mat2%p,&
                       FockI,FockA,mat2%A/2.0d0,U,Dia)
! Hessian for 2index ss, diagional elements are used
          call precabb(orb%nsub,orb%occ,orb%closed,orb%total,&
                       norb,nocc,ncls,mat2%D,mat2%p,&
                       FockI,FockA,mat2%A/2.0d0,U,Dsa)

          write(126,*)"Hess_ia,ia"
          call print_mat(ncls,nact,Dia,126)
          write(126,*)"Hess_sa,sa"
          call print_mat(nvir,nact,Dsa,126)

          nacttmp=nact

! closed shell update
          call CS_dim_cls2fro()

          write(1,*)"ncls",ncls

            do j=1+ncls,nacttmp+ncls
              do i=1,ncls
                write(1,*)"Corrected Hess(ji,ji)",j,i,j,i,mat2%H(j,i,j,i),Dia(i,j-ncls)
                mat2%H(j,i,j,i)= Dia(i,j-ncls)
!              mat2%H(j,i,i,j)=-Dia(i,j-ncls)
!              mat2%H(i,j,j,i)=-Dia(i,j-ncls)
              end do
            end do

!            do j=1+ncls,nacttmp+ncls
!              do i=1+nocc,nvir+nocc
!                write(1,*)"Corrected Hess(ji,ji)",j,i,j,i,mat2%H(j,i,j,i),Dsa(i-nocc,j-ncls)
!                mat2%H(j,i,j,i)= Dsa(i-nocc,j-ncls)
!!              mat2%H(j,i,i,j)=-Dia(i,j-ncls)
!!              mat2%H(i,j,j,i)=-Dia(i,j-ncls)
!              end do
!            end do

          deallocate(Dia,Dsa)

!         stop
        else if(iHessian.eq.2)then

          call Dh_gen(orb%nsub,orb%act,orb%total,&
                     nact,norb,T,mat2%D,mat2%Dh,orb%grouptable)

          write(2,*)" The 2*D*h part "
          call flush(2)
!          call print_gat(norb,norb,norb,norb,mat2%Dh,2)

          ! one-by-one Hessian element calculation
          i0=0
          do isym=1,orb%nsub
            do i=1,orb%act(isym)
              j0=0
              do jsym=1,orb%nsub
                do j=1,orb%act(jsym)
                  mat2%Y1=0.0d0
                  call PmultU(orb%nsub,orb%act,orb%total,nact,norb,&
                       U,mat2%P,mat2%Y1,orb%grouptable,isym,i,jsym,j,1)
                  mat2%Y2=mat2%Y1

                  write(2,*)"The squared PmultU type 1 outsides ===="
                  call print_mat(norb,norb,mat2%Y2,2)
                  write(2,*)"The squared PmultU type 1 outsides ===="

                  call PmultU(orb%nsub,orb%act,orb%total,nact,norb,&
                       U,mat2%P,mat2%Y2,orb%grouptable,isym,i,jsym,j,2)

                  write(2,*)"=========================================="
                  write(2,*)" Finished ",i,j,"-th PmultU calculations "
                  write(2,*)"=========================================="

                  write(2,*)"The squared PmultU ==== outsides ===="
                  call print_mat(norb,norb,mat2%Y2,2)
                  write(2,*)"The squared PmultU ==== outsides ===="
                  write(2,*)" "

                  mat2%Y(i0+i,:,j0+j,:)=mat2%Y2 ! origin
                  !mat2%Y(:,i0+i,:,j0+j)=mat2%Y2

                end do
                j0=j0+orb%total(jsym)
              end do
            end do
            i0=i0+orb%total(isym)
          end do

          write(2,*)"The mat2%Y  "
          call print_gat(norb,norb,norb,norb,mat2%Y,323)

          call Hessian3(norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)
!          call Hessian_risj&
!               (norb,mat2%Dh,mat2%A,mat2%Y,mat2%H,mat2%Hdiag)

          call print_gat(norb,norb,norb,norb,mat2%H,325)

        end if

! to Leon :  the following:
!         Construct the couling for every states

!        write(2,*)" ======== Print use the following format ======== "
!        write(2,*)"   iroot+1, iC, i, j, MPS(iroot+1)%HCR(iC,i,j)"
!        write(2,*)" ================================================ "

! Option-1: construct the coupling base on 1) MPS_Davidson_elements --> 2) orbital_rotations

!        write(6,*)dmrg_nstates,MPS(iroot+1)%nC,norb,norb
        call flush(6)

        if(.true.)then
          do iroot=0,dmrg_nstates-1
            do iC=1,MPS(iroot+1)%nC

!            write(2,*)"i-i TDMS (actually RDMs)"
!            call print_mat(norb,norb,MPSTDMs(iroot+1,iroot+1)%A,2)

!            write(2,*)iroot,"-th root ",iC,"-th RDMs deri)"
!            call print_mat(norb,norb,MPS(iroot+1)%Ci(iC)%A,2)
              allocate(TM3(norb,norb))
              allocate(TM4(norb,norb))
              TM3=0.0d0
              TM4=0.0d0

              do i=1,norb
                do j=1,norb

                  iipos=MPS(iroot+1)%Ci(iC)%ipos

!                  write(6,*)"iipos",iipos,"L_R%nC",L_R%nC
!                  write(6,*)"L_R%idx",L_R%idx
                  call get_position(iipos,L_R%nC,L_R%idx,L_R%pos,iCn,0)
!                  write(6,*)"iCn",iCn
!                  call flush(6)

                  !iCn=MPS(iroot+1)%Ci(iC)%n
                  MPS(iroot+1)%KCR(iCn,i,j)=MPS(iroot+1)%Ci(iC)%A(i,j)&
                             -MPS(iroot+1)%Ci(iC)%A(j,i) ! For property check

!                do jroot=0,dmrg_nstates-1 ! Asuming orthogonality

                  MPS(iroot+1)%HCR(iCn,i,j)=MPS(iroot+1)%Ci(iC)%A(i,j)&
                                           -MPS(iroot+1)%Ci(iC)%A(j,i)&
                -MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(i,j)&
                +MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(j,i) ! <-- second line ()

                  TM3(i,j)=&
                  MPS(iroot+1)%Ci(iC)%A(i,j)-MPS(iroot+1)%Ci(iC)%A(j,i)

                  TM4(i,j)=&
                 -MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(i,j)&
                 +MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(j,i)

!               end do
                  MPS(iroot+1)%HCR(iCn,i,j)= &
                  MPS(iroot+1)%HCR(iCn,i,j)*dmrg_weight(iroot+1)*(1.0d0)! sign problem

!               write(*,*)" does more than two states also need TDM? no"
!                 write(2,*)iroot+1,iC,i,j,MPS(iroot+1)%HCR(iC,i,j)

                end do
              end do

!              write(2,*)"The A^I-A^I^dagger for state",iroot,"and iC",ic
!              call print_mat1(norb,norb,TM3,2)
!              write(2,*)"The Ci(A-A^dagger) for state",iroot,"and iC",ic
!              call print_mat1(norb,norb,TM4,2)
              write(2,*)"The non-transtion part"
              call print_mat1(norb,norb,TM3+TM4,2)

              deallocate(TM3,TM4)

!              allocate(TM3(norb,norb))
!              allocate(TM4(norb,norb))
!              TM3=0.0d0
!              TM4=0.0d0

              allocate(TM5(norb,norb))
              allocate(TM6(norb,norb))
              TM5=0.0d0
              TM6=0.0d0

!              allocate(TM1(norb,norb))
!              allocate(TM2(norb,norb))
!              TM1=0.0d0
!              TM2=0.0d0

              do jroot=0,dmrg_nstates-1
                if(jroot.ne.iroot)then
!                  irt=iroot+1
!                  jrt=jroot+1

                  ! check if this MPS element exist
!                  jexist=0
!                  do jCx=1,MPS(jroot+1)%nC
!                    if(MPS(jrt)%Ci(jCx)%num.eq.MPS(irt)%Ci(iC)%num)then
!                      jexist=1
!                    end if
!                  end do
!                  if(jexist.eq.1)then

                    do i=1,norb
                      do j=1,norb

!                TM3(i,j)=MPSTDMs(iroot+1,jroot+1)%A(i,j) !  Ax
!
!                TM4(i,j)=&
!                MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(i,j)
!                      TM1(i,j)=&
!           -MPSTDMs(iroot+1,iroot+1)%A(i,j)&
!           +MPSTDMs(iroot+1,iroot+1)%A(j,i) ! <-- second line ()
!                      TM2(i,j)=&
!           -MPSTDMs(jroot+1,jroot+1)%A(i,j)&
!           +MPSTDMs(jroot+1,jroot+1)%A(j,i) ! <-- second line ()

!                      TM5(i,j)=&
!                  -MPSTDMs(iroot+1,jroot+1)%A(i,j)& !  Ax
!                  +MPSTDMs(iroot+1,jroot+1)%A(j,i)  !
                      TM6(i,j)=&
            MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(i,j)&
           -MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(j,i)
!           -MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(i,j)&
!           +MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,iroot+1)%A(j,i) ! <-- second line ()
!                  -MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(i,j)&
!                  +MPS(iroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(j,i)
!                    MPS(iroot+1)%HCR(iC,i,j)=MPS(iroot+1)%HCR(iC,i,j)&
!              -MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(i,j)&
!              +MPS(jroot+1)%Ci(iC)%dv*MPSTDMs(iroot+1,jroot+1)%A(j,i)

                 MPS(iroot+1)%HCR(iCn,i,j)= &
                 MPS(iroot+1)%HCR(iCn,i,j)-TM6(i,j)*dmrg_weight(jroot+1)

                      end do
                    end do
!                  end if  ! if jroot has the same MPS element as that of iroot

                end if
              end do

              write(2,*)"The cI(A-A^dagger(TDM)) for state",iroot
              call print_mat1(norb,norb,TM6,2)

              write(2,*)"The coupling for state",iroot,"element",iC
            call print_mat1(norb,norb,2.0d0*MPS(iroot+1)%HCR(iC,:,:),2)

!              deallocate(TM1,TM2)
!              deallocate(TM3,TM4)
              deallocate(TM5,TM6)

        write(2,*)" ================================================ "
        write(2,*)" ================================================ "
        write(2,*)" ================================================ "

!            if(ic.eq.1)then
!              stop
!            end if

            end do

!            MPS(iroot+1)%HCR = MPS(iroot+1)%HCR*(-2.0d0)
!     ??? minus(-) or plus(+)
            MPS(iroot+1)%HCR = MPS(iroot+1)%HCR*( 2.0d0)
            write(1,*)"Need check if need minus"
            write(1,*)"Need check if need minus"
            write(1,*)"Need check if need minus"
            write(1,*)"Need check if need minus"
            write(1,*)"Need check if need minus"
!          MPS(iroot+1)%HCR = MPS(iroot+1)%HCR*(-2.0d0)
!          call print_tat(MPS(iroot+1)%nC,norb,norb,MPS(iroot+1)%HCR,12)
          end do
! Option-2: construct the coupling base on 1) orbital_rotations --> 2) MPS_Davidson_elements
!        stop
       else

       end if

! to Leon :
! ====================================================================================
! Construct the full LR equations (left part, no consideration of redundant variables)
! ====================================================================================

! Assuming the whole Hessian is needed by ndim1 = norb
        ndim1=norb
        ndim2=ndim1**2
        ndimH=ndim2
        do iroot=0,dmrg_nstates-1
!          ndimH=ndimH+MPS(iroot+1)%nC
          ndimH=ndimH+L_R%nC
        end do
        write(2,*)"dim for the full Hessian", ndimH
        call flush(2)
        allocate(H(ndimH,ndimH))
        H=0.0d0

! 1) Put orb-orb Hessian
        ij=0
        do i=1,norb
          do j=1,norb
            ij=ij+1
            kl=0
            do k=1,norb
              do l=1,norb
                kl=kl+1
                H(ij,kl)=mat2%H(i,j,k,l)
              end do
            end do
          end do
        end do
        write(2,*)"After put orb-orb coupling in the full Hessian"
        call flush(2)

!        call print_mat(ndimH,ndimH,H,13)
!        call flush(13)

!        call print_gat(norb,norb,norb,norb,mat2%H,20)
!        call flush(20)

! 2) Put orb-MPS coupling
       ! offset value, use for cycle the orbital-orbital part
        offset=norb*norb

        do iroot=0,dmrg_nstates-1
         
          do iC=1,L_R%nC     !!!!!!!!   nMPS_LIST
            ij=0
            do i=1,ndim1
              do j=1,ndim1
                ij=ij+1
                H(offset+iC,ij)=MPS(iroot+1)%HCR(iC,i,j)
                H(ij,offset+iC)=MPS(iroot+1)%HCR(iC,i,j)
              end do
            end do
          end do
          offset=offset+MPS(iroot+1)%nC
        end do
        write(2,*)"After put orb-MPS coupling in the full Hessian"
        call flush(2)

! 3) Put the MPS-MPS coupling
       ! offset value, use for cycle the orbital-orbital part
        offset=norb*norb

        do iroot=0,dmrg_nstates-1
          !call print_mat(MPS(iroot+1)%nC,MPS(iroot+1)%nC,MPS(iroot+1)%HCC,1)
          !call flush(1)
          do i=1,L_R%nC
            do j=1,L_R%nC
              H(offset+i,offset+j)=MPS(iroot+1)%HCC(i,j)
            end do
          end do
          offset=offset+L_R%nC
        end do
        write(2,*)"After put MPS-MPS coupling in the full Hessian"
        call flush(2)

!        call print_mat(ndimH,ndimH,H,14)
!        call flush(14)
!        stop

        do i=1,ndimH
          do j=1,ndimH
!            write(14,*)"Full(modified) Hessian",i,j,H(i,j)
          end do
        end do

! =============================================================
!        Only the non-redundant parameters are needed
! =============================================================
! to Leon :!       The nij_add is caused by Hess_ia,ia != Hess_ai,ai


        ! Add solver selection later
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"
        write(6,*)&
        "Notice : nij_add is set to 0 in order to use the whole solver"

        nij_add=0
        nij=(nij0+nij_add)/2 ! remove the redundant from anti-symmetry

        ! non-redundant index
        allocate(porb(nij));   porb=0
        ! redundant but closed shell, !!-> index in porb
        if(nij_add.gt.0)then
          allocate(porb_redu(nij_add/2,2)); porb_redu=0
          allocate(porb_tmp(nij_add/2));    porb_tmp=0
        else ! for mem issue
          allocate(porb_redu(1,2)); porb_redu=0
          allocate(porb_tmp(1));    porb_tmp=0
        end if
        write(6,*)"porb_redu",porb_redu
        write(6,*)"porb_tmp ",porb_tmp
        call flush(6)

        nC0=0
        do iroot=0,dmrg_nstates-1
          nC0=nC0+L_R%nC
        end do

        allocate(L_R%Y(nij));        L_R%Y  =0.0d0
        allocate(L_R%HOO(nij,nij));  L_R%HOO=0.0d0
        allocate(L_R%HCO(nC0,nij));  L_R%HCO=0.0d0
        allocate(L_R%HOC(nij,nC0));  L_R%HOC=0.0d0
        allocate(L_R%HCC(nC0,nC0));  L_R%HCC=0.0d0

        ! The full orbital Lagrange
        allocate(L_R%R(norb,norb));  L_R%R  =0.0d0

        ! If .false., the gradients from Molcas will used
        Lcorrected_Phase=.false.
        if(.not.Lcorrected_Phase)then     !  1/2   IfcorrectedRHS
          !!! Leon: for debugging, read RHS calculated by MOLCAS
          !!! Make sure you're using this with the MOLCAS version
          !!! that produces this file!
          open(unit=100, file="OrbGrad.txt")
          ij=0
          do i=1,norb
            do j=1,norb
                ij=ij+1
              read (100,*) G(ij)
            end do
          end do
          close(100)
          !!!
          !!! Leon: full RHS for debugging
          write (1,*) "================================================"
          write (1,*) " Full orbital gradient is imported from Molcas"
          write (1,*) "================================================"
          ij=0
          do i=1,norb
            do j=1,norb
              ij=ij+1
              write (1,*) G(ij)
            end do
          end do
        else
          write (1,*)"================================================="
          write (1,*)"The phase for MPS-LAG and Orb-LAG will be altered"
          write (1,*)"================================================="
        end if

        ! This is the delta-orbiatl-rotation gradients
        ij=0
        kl=0
        iredu=0
        do i=1,norb
          do j=1,norb ! only half is needed
            ijkl=0
            ij=ij+1
            !  2/2   IfcorrectedRHS
            if(.not.Lcorrected_phase)then
            else
              G(ij)=MPS(rlx)%A(i,j)-MPS(rlx)%A(j,i)
            end if
            !Check the valid rotations
            Lvalid=.false.
            do k=1,nij0/2
              if(i.eq.valid(k,1).and.j.eq.valid(k,2))then
                Lvalid=.true.
                exit
              end if
            end do
            do k=nij_add/2+1,nij_add
              if(i.eq.valid_add(k,2).and.j.eq.valid_add(k,1))then
                Lvalid=.true.
                ijkl=1
                exit
              end if
            end do

            ! act
            if(Lvalid.and.ijkl.eq.0)then
              kl=kl+1
              porb(kl)=ij
            end if
            ! closed additional
            if(Lvalid.and.ijkl.eq.1)then
              kl=kl+1
              kl1=0
              do i0=1,norb
                do j0=1,norb
                  kl1=kl1+1
                  !if(i0.eq.valid_add(k,2).and.j0.eq.valid_add(k,1))then ! i,j re-counter
                  if(i0.eq.valid_add(k,2).and.j0.eq.valid_add(k,1))then ! i,j
                    ij1=kl1
                    !write(*,*)valid_add(k,1),valid_add(k,2),ij1,kl
                    iredu=iredu+1
                  end if
                  if(j0.eq.valid_add(k,2).and.i0.eq.valid_add(k,1))then ! j,i
                    ij2=kl1
                  end if
                end do
              end do
              porb(kl)=ij1
              porb_tmp(iredu)=ij2
              porb_redu(iredu,1)=kl
            end if

          end do
        end do

       ! Write down the gradients or RHS
        open(unit=110,file="RHS_DMRG_MCLR.txt")
          ij=0
          do i=1,ncls
            do j=ncls+1,norb
              ij=ij+1
              write(110,*)ij,MPS(rlx)%A(i,j)-MPS(rlx)%A(j,i)
            end do
          end do
          do i=ncls+1,nocc
            do j=nocc+1,norb
              ij=ij+1
              write(110,*)ij,MPS(rlx)%A(i,j)-MPS(rlx)%A(j,i)
            end do
          end do
        close(110)


        ij=0
        do i=1,nij
          do j=1,nij_add/2
            if(porb_tmp(j).eq.porb(i))then
              ij=ij+1
              porb_redu(ij,2)=i
            end if
          end do
        end do

        write(1,*)""
        write(1,*)"porb : "
        write(1,*) porb
        if(nij_add.gt.0)then
          write(1,*)"porb_redu : "
          write(1,*) porb_redu(:,1)  !     redundant index
          write(1,*) porb_redu(:,2)  ! non-redundant index (counter-part in prob)
        end if
        write(1,*)" ===== "
        call flush(1)
        deallocate(porb_tmp)

! to Leon :!      only the effective part of Hessian will be constructed

        ! non-redundant orbital-orbital part
        do i=1,nij
          L_R%Y(i)= G(porb(i))
          do j=1,nij    ! needed Hessian
            L_R%HOO(i,j)=H(porb(i),porb(j))
          end do
        end do

        write(1,*)"Effective Hessian (orb-orb)"
        call print_mat1(nij,nij,L_R%HOO,120)
!        call print_mat(nij,nij,L_R%HOO,6)
!        write(6,*)"===="

! Notice Y .ne. G
        write(1,*)"RHS (residual) : "
        write(1,*) -1.0d0*L_R%Y

        offset=norb*norb
        ! non-redundant orbital-MPS part
        ici=0
        i0=0
        do iroot=0,dmrg_nstates-1
          !nC=MPS(iroot+1)%nC
          nC=L_R%nC
          do iC=1,nC
            ici=ici+1
            do i=1,nij
              L_R%HCO(ici,i)=H(offset+iC+i0,porb(i))
              L_R%HOC(i,ici)=H(porb(i),offset+iC+i0)
            end do
          end do
          i0=i0+nC
        end do

        write(1,*)" LR%HCO "
!        call print_mat(nC0,nij,L_R%HCO,1)
        call flush(1)
        write(1,*)" LR%HOC "
        call print_mat1(nij,nC0,L_R%HOC,1)
        call flush(1)

        ! non-redundant MPS-MPS part
        ici=0
        i0=0
        do iroot=0,dmrg_nstates-1
          !nC=MPS(iroot+1)%nC
          nC=L_R%nC
          do iC=1,nC
            do jC=1,nC
              L_R%HCC(iC+i0,jC+i0)=H(offset+iC+i0,offset+jC+i0)
            end do
          end do
          i0=i0+nC
        end do

        write(1,*)" LR%HCC "
        call print_mat1(nC0,nC0,L_R%HCC,1)
        call flush(1)
        !stop "fort.1 check for the Hessian (to do) "

        allocate( LP(nij+nC0));          LP=0.0d0
        allocate(LP0(nij+nC0));         LP0=0.0d0
        allocate(LPx(nij+nC0,100));     LPx=0.0d0
        allocate(  YP(nij+nC0));         YP=0.0d0
        allocate( YP0(nij+nC0));        YP0=0.0d0
        allocate(YP00(nij+nC0));       YP00=0.0d0
        allocate(HP(nij+nC0,nij+nC0));   HP=0.0d0
        allocate(HP2(nij+nC0,nij+nC0)); HP2=0.0d0

        YP(1:nij)                       = -1.0d0*L_R%Y   ! Should be minus due to -g
        HP(1:nij,1:nij)                 =        L_R%HOO
        HP(nij+1:nij+nC0,1:nij)         =  1.0d0*L_R%HCO
        HP(1:nij,nij+1:nij+nC0)         =  1.0d0*L_R%HOC
        HP(nij+1:nij+nC0,nij+1:nij+nC0) =  1.0d0*L_R%HCC

        HP2=HP  ! Save for backup and recover

!        call print_mat(nij,nij,HP2(1:nij,1:nij),666)

        write(1111,*)" The  SA Hessian : "
        call print_mat1(nij+nC0,nij+nC0,HP,1111)
        ! Change to Davidson when Hessian is too large
        call eigtql2(nij+nC0,HP,YP00)
        write(1111,*)" The  eigenvalues of SA Hessian : "
        call print_mat1(1,nij+nC0,YP00,1111)
        call flush(1111)
!        write(1111,*)" The eigenvectors of SA Hessian : "
!        call print_mat1(nij+nC0,nij+nC0,HP,1111)
!        call flush(1111)

        write(6,*)"Print the lowest 5 eigenvalues of full-SA Hessian : "
        call print_mat1(1,5,YP00(1:5),6)
        call flush(6)
        write(2,*)"Print the lowest 5 eigenvalues of full-SA Hessian : "
        call print_mat1(1,5,YP00(1:5),2)
        call flush(2)

        if(YP00(1).lt.-1.0d-5)then
          write(6,*)"Warning : The Hessian is not positive-defined "
          write(6,*)"          Print the lowest 5 eigenvalues : "
          call print_mat1(1,5,YP00(1:5),6)
          write(6,*)"Warning : The Hessian is not positive-defined "
          write(6,*)"If SA states in a same irreps, it implies one "
          write(6,*)"may need to correct the sign in Hess(orb-MPS) "
          !stop
        end if

        HP=HP2 ! first time recover back the Hessian
        YP00=0.0d0

!        stop

        if(.true.)then
          write(6,*)"Notice, the non-preconditioner &
      Hessian elements are set to 0 as inital"

          do i=1,norb
            ii=0
            do j=1,nij0/2
              if(valid(j,1).eq.i)then
                ii=ii+1
              end if
              nprei(i)=ii
            end do
            do j=1,nij_add/2
              if(valid_add(j,1).eq.i)then
                ii=ii+1
              end if
              nprei(i)=ii
            end do
          end do

          write(6,*)"nprei",nprei
          call flush(6)

          ! Diagional as preconditioner
          HP(1:nij,1:nij)=0.0d0

          ii1=1
          do i=1,norb
            ii2=ii1+nprei(i)-1

!            write(6,*)"ii1,ii2",ii1,ii2
!            call flush(6)
            if(ii1.le.ii2)then
!              write(6,*)ii1,ii2,HP2(ii1:ii2,ii1:ii2)
!              call flush(6)
              HP(ii1:ii2,ii1:ii2)=HP2(ii1:ii2,ii1:ii2)
            end if

            if(nprei(i).gt.0.and.nprei(i+1).gt.0)then
              ii1=ii2+1
            end if
            if(nprei(i).gt.0.and.nprei(i+1).eq.0)then
              ii1=ii2+1
            end if
          end do

          write(116,*)"The dimension for HP",nij
          call print_mat(nij,nij,HP(1:nij,1:nij),116)
          call flush(116)

        end if

        if(allocated(valid_add))then
          deallocate(valid_add)
        end if

! to Leon :!        Check for the RHS, i.e. orbital gradients for a specific state

!        write(1,*)" LHS "
!        call print_mat(nij+nC0,nij+nC0,HP,1)


        write(1,*)" RHS "
        call print_mat(      1,nij+nC0,YP,1)
        call flush(1)

! to Leon :!        Start the PCG solver

        allocate(SAvec(nC0));               SAvec=0.0d0
        allocate(HPM(nij+nC0,nij+nC0));     HPM=0.0d0
        isa=0
        do iroot=0,dmrg_nstates-1
          !do iC=1,MPS(iroot+1)%nC
          do iC=1,L_R%nC
            isa=isa+1
            SAvec(isa)=MPS(iroot+1)%Ci(iC)%dv
!           write(1,*)ic,MPS(iroot+1)%Ci(iC)%dv,SAvec(nij+isa)
          end do
        end do

        do i=1,nij+nC0
          !HPM(i,i)=1.0/HP(i,i)
          HPM(i,i)=HP(i,i)
        end do

        neff=nij+nC0-nij_add/2

        write(6,*)"neff",neff
        call flush(6)

        allocate(sigma(neff));sigma=0.0d0
        allocate( xvec(neff)); xvec=0.0d0
        allocate( rvec(neff)); rvec=0.0d0

        allocate(   V1(neff));  V1 =0.0d0
        allocate(   V2(neff));  V2 =0.0d0
        allocate(   V3(neff));  V3 =0.0d0

        iloop=0
        !iloopmax=max_iter
        iloopmax=1
        do

          iloop=iloop+1

          ij=0
          kl=0
!           do i=1,norb
!             do j=1,norb ! only half is needed
!               ij=ij+1
! !               G(ij)=MPS(rlx)%A(i,j)-MPS(rlx)%A(j,i)
!             end do
!           end do

!          do i=1,nij
!            L_R%Y(i)= G(porb(i))
!          end do
!          YP(1:nij)                       = -1.0d0*L_R%Y

          if(.true.)then

            write(1,*)"SAvec "
!         call print_mat(1,nij+nC0,SAvec,1)

          !call Solver_LR(nij+nC0,HP,YP,LP)
          !call Solver_CG_SA(nij+nC0,HP,YP,LP,SAvec)
          ! Solved results from Molcas

            if(.not.allocated(conf))then
              allocate(conf(dmrg_nstates))
              conf=0
            else
              conf=0
            end if

            do iroot=0,dmrg_nstates-1
              conf(iroot+1)=MPS(iroot+1)%nC
            end do

            !call Solver_LR(nij+nC0,HP,YP,LP)
            !call Solver_CG(nij+nC0,HP,YP,LP)
            !call Solver_PCG(nij+nC0,HP,YP,HPM,LP)
            !call Solver_PCG_SA(nij+nC0,HP,YP,HPM,LP,SAvec)

! Solve as a whole
            if(.true.)then

              write(6,*)"Solve as a whole in line 1348 of ConDMRGLR.f90"
              call flush(6)

              write(6,*) " == nij == ",nij
!              call print_mat(nij,nij,HP ,6)
              write(6,*)"Printed HP"
              call flush(6)
!              call print_mat(nij,nij,YP ,6)
              write(6,*)"Printed YP"
              call flush(6)
!              call print_mat(nij,nij,HPM,6)
              write(6,*)"Printed HPM"
              call flush(6)
!              call print_mat(nij,nij,LP ,6)
              write(6,*)"Printed LP "
              call flush(6)

! If couple the MPS/ci or not
!            if(iloop.eq.1)then ! no need to consider MPS/ci part
              call Solver_PCG&
              (nij,HP(1:nij,1:nij),YP(1:nij),HPM(1:nij,1:nij),LP(1:nij))
!            else               !

              LP0=LP
              write(6,*)"Finished the orbital-only PCG solver"
              call flush(6)

              ! ==========================================
              !     Write down the orbital-lag

               allocate(TM1(norb,norb)); TM1=0.0d0
! Lagrange - orbital part
               do i=1,nij
                 i1=valid(i,1)
                 j1=valid(i,2)
!! recover the anti-symmetric rotations
                 TM1(i1,j1)=  LP(i)/1.0d0
                 TM1(j1,i1)= -LP(i)/1.0d0
               end do

               open(unit=110,file="orb_lag_initial.txt")
                 ij=0
                 do i=1,ncls
                   do j=ncls+1,norb
                     ij=ij+1
                     write(110,*)ij,TM1(i,j)
                   end do
                 end do
                 do i=ncls+1,nocc
                   do j=nocc+1,norb
                     ij=ij+1
                     write(110,*)ij,TM1(i,j)
                   end do
                 end do
               close(110)

               deallocate(TM1)
              ! ==========================================


              HP=HP2 ! In case the Hessian is modofied

              write(6,*)"Finished the HP=HP2 "
              call flush(6)

              if(.false.)then

                allocate(TM1(nij+nC0,nij+nC0));     TM1=0.0d0
                call precondition(nij,nC0,HPM,TM1,dmrg_nstates,conf,SAvec)

                write(31,*)"HPM"
                call print_mat(nij+nC0,nij+nC0,HPM,31)
                write(31,*)"TM1"
                call print_mat(nij+nC0,nij+nC0,TM1,31)

                ! The modified preconditioner (TM1)
                call Solver_PCG_SA2&
                     (nij,nC0,HP,YP,LP,HPM,dmrg_nstates,conf,SAvec)

                deallocate(TM1)
              else
                write(6,*)"Start the Solver_CG_SA"
                call flush(6)
                call Solver_CG_SA&
                     (nij,nC0,HP,YP,LP,dmrg_nstates,conf,SAvec)
!                call Solver_PCG&
!                     (nij+nC0,HP,YP,HPM,LP)
                write(6,*)"Finished the Solver_CG_SA"
                call flush(6)
              end if

            else

! Solve as partions, i.e. -> Hoo-Roo -> Hoc-Rco -> Hcc-Rcc ->

              ! Total number of (P)CG iterations
              nCG=nij+nC0 !(tempoary OK, should be all dimension)
              YP00=YP

              ! initialize parameters that used to control PCG/CG process
              rAlphaORB=0.0d0
              rdeltaORB=0.0d0
              rsigmaORB=0.0d0
              resORB=0.0d0

              rAlphaMPS=0.0d0
              rdeltaMPS=0.0d0
              rsigmaMPS=0.0d0
              resMPS=0.0d0

              rAlpha = 0.0d0
              rdelta = 0.0d0
              res=0.0d0

              ! The start iteration (only orbital lags)

              do i=1,1

                write(1,*)" ======================================= "
                write(1,*)"    The residual in iteration ",i
                write(1,*)" ======================================= "
                write(1,*)" "

                 LP=0.0d0
                YP0=0.0d0

                ! Hoo-Ro=Residual -> get Ro
                call Solver_PCG&
             (nij, HP(1:nij,1:nij),YP(1:nij),HPM(1:nij,1:nij),LP(1:nij))

                write(1,*)"The Rotations (LP) orb"
                call print_mat(1,nij,LP(1:nij),1)

                goto 2222

                  ! H*R check (Ax check, automatically remove the redundant)
                 call Hall_Rall&
                     (nij,nC0,HP2,LP,YP0,nij_add/2,porb_redu,sigma,xvec)

                  call rm_redundant(nij,nC0,nij_add/2,porb_redu,YP,V1)
                  call inner_product(neff,xvec,   V1,rdelta)

                  ! alpha - parameter
                  call inner_product(neff,xvec, xvec,dv1)
                  call inner_product(neff,xvec,sigma,dv2)
                  write(1,*)"<x|x>=",dv1," <x|A|x>=",dv2

                  ralpha=rdelta/dv2
                  write(1,*)"ralpha",ralpha,"rdelta (<x|r_old>)",rdelta

                  call flush(1)

                  ! update residual
                  call RminusAx&
                 (nij,nC0,YP0,nij_add/2,porb_redu,YP,rvec,ralpha)

                  write(1,*)"The Residual YP (modified redundant)"
                  call print_mat(1,nij+nC0,YP,1)

                  ! Save R, and updated R for the next iterations
                  LPx(:,i)=ralpha*LP+LPx(:,i)
                  LP=LPx(:,i)

                  write(1,*)"The rotations LP (modified redundant)"
                  call print_mat(1,nij+nC0,LP,1)

                  call rm_redundant(nij,nC0,nij_add/2,porb_redu,LP,V1)
                  call inner_product(neff, V1,   V1,  res)

!                ! H*R check (Ax check, automatically remove the redundant)
!                call Hall_Rall&
!               (nij,nC0,HP2,LP,YP0,nij_add/2,porb_redu,sigma,xvec)

!                ! update residual
!                YP=YP00
!                call RminusAx&
!               (nij,nC0,YP0,nij_add/2,porb_redu,YP,rvec)

                  YP00=YP

!                write(1,*)"The Residual YP (modified redundant),iter-1 "
!                call print_mat(1,nij+nC0,YP,1)
!                write(1,*)"The Residual YP ( without redundant),iter-1 "
!                call print_mat(1,neff,rvec,1)

2222             continue

              end do

              write(2,*)&
              "Finish the diagional LR-orb part as the initial guess"
!              call flush(2)
!              stop

              naddCG=20

              do i=1,nCG+naddCG
              !do i=1,5

                write(1,*)" ======================================= "
                write(1,*)"    The residual (RHS) in iteration ",i
                write(1,*)" ======================================= "
                write(1,*)" "
                call print_mat(1,nij+nC0,YP,1)

                ! H*R check (Ax check, automatically remove the redundant)
                call Hall_Rall&
                    (nij,nC0,HP2,LP,YP0,nij_add/2,porb_redu,sigma,xvec)

                ! alpha/beta/delta - parameter
                call rm_redundant(nij,nC0,nij_add/2,porb_redu,YP,V1)

                write(1,*)"The Residual YP   (non-redundant)"
                call print_mat(1,neff,  V1,1)
                write(1,*)"The LR parameters (non-redundant)"
                !call print_mat(1,neff,xvec,1)
                call print_mat1(1,neff,xvec,1)

                ! orbital first
                call inner_product&
                (neff-nC0,xvec(1:neff-nC0),  V1(1:neff-nC0),rdeltaorb)
                call inner_product&
                (neff-nC0,xvec(1:neff-nC0),xvec(1:neff-nC0),resorb)
                call inner_product&
                (neff-nC0,xvec(1:neff-nC0),sigma(1:neff-nC0),dv1)

                ! MPS second
                call inner_product&
            (nC0,xvec(neff-nC0+1:neff),   V1(neff-nC0+1:neff),rdeltaMPS)
                call inner_product&
            (nC0,xvec(neff-nC0+1:neff), xvec(neff-nC0+1:neff),resMPS)
                call inner_product&
            (nC0,xvec(neff-nC0+1:neff),sigma(neff-nC0+1:neff),dv2)

                ralphaorb=rdeltaorb/dv1
                ralphamps=rdeltamps/dv2

                rdelta=rdeltaORB+rdeltaMPS
                ralpha=rdelta/(dv1+dv2)

!                rbeta=(rdeltaORB+rdeltaMPS)/rdelta
!                rdelta=rdeltaORB+rdeltaMPS

                write(1,*)"orb: <x|x>=",resorb," <x|A|x>=",dv1,&
                          "rdelta (<x|r_old>)",rdeltaorb
                write(1,*)"MPS: <x|x>=",resMPS," <x|A|x>=",dv2,&
                          "rdelta (<x|r_old>)",rdeltaMPS

                write(1,*)&
              "ralpha",ralpha,"rbeta",rbeta,"rdelta (<x|r_old>)",rdelta

                ! update residual
                call RminusAx&
               (nij,nC0,YP0,nij_add/2,porb_redu,YP,rvec,ralpha)

                write(1,*)"The Residual YP (modified redundant)"
                call print_mat1(1,nij+nC0,YP,1)

                ! Save R, and updated R for the next iterations
                LPx(:,i)=ralpha*LP+LPx(:,i)

                write(1,*)"The LR parameters (modified redundant)"
                call print_mat1(1,nij+nC0,LPx(:,i),1)

!                stop

                ! LP=0.0d0
                 LP=0.0d0
                YP0=0.0d0

                !if(i.gt.1)then

                ! ================= The MPS part ===================

                  ! Note the HPM

                  if(.false.)then
                    ! Hcc-Rc=Residual -> get Rc
                    call Solver_PCG&
                        (nC0,HP(nij+1:nij+nC0,nij+1:nij+nC0)/1.0d0,  &
                             YP(nij+1:nij+nC0),                &
                             HPM(nij+1:nij+nC0,nij+1:nij+nC0)/1.0d0,  &
                              LP(nij+1:nij+nC0))
                  end if
                  if(.false.)then

                    ! If need the preconditioned version ??
                    !call Solver_PCG_SA(ndim,LHS,RHS,PM,VEC,SAvec)

                  else

                    call Solver_CG_SA&
                         (0,nC0,HP(nij+1:nij+nC0,nij+1:nij+nC0),  &
                            YP(nij+1:nij+nC0),LP(nij+1:nij+nC0),  &
                             dmrg_nstates,conf,SAvec)
                  end if

                  write(1,*)"====================================="
                  write(1,*)"  The rotations parameters of MPS"
                  call print_mat(1,nC0,LP(nij+1:nij+nC0),1)
                  write(1,*)"====================================="
!                  write(1,*)LP(nij+1:nij+nC0)

!                  stop

                !  LPx(:,i)=LP

                ! ================= The MPS part ===================

                !end if

                !LP=0.0d0
                ! ================= The ORB part ===================

                ! Hoo-Ro=Residual -> get Ro
                call Solver_PCG&
             (nij, HP(1:nij,1:nij),YP(1:nij),HPM(1:nij,1:nij),LP(1:nij))

                write(1,*)"The Rotations (LP) orb"
                call print_mat(1,nij,LP(1:nij),1)

                ! ================= The ORB part ===================

                ! updated the X-vector (LR)
                call rm_redundant(nij,nC0,nij_add/2,porb_redu,LP,xvec)

                write(1,*)"The updated Residual   (non-redundant)"
                call print_mat(1,neff,rvec,1)
                write(1,*)&

               "The LR parameters, after preconditioner (non-redundant)"
                call print_mat(1,neff,xvec,1)

                ! orbital
                call inner_product&
                (neff-nC0,xvec(1:neff-nC0),rvec(1:neff-nC0),rdeltaorb)
                ! MPS
                call inner_product&
           (nC0,xvec(neff-nC0+1:neff),rvec(neff-nC0+1:neff),rdeltaMPS)


                rbeta=(rdeltaORB+rdeltaMPS)/rdelta
                rdelta=rdeltaORB+rdeltaMPS

                write(1,*)&
            "rdeltaORB",rdeltaORB,"rdeltaMPS",rdeltaMPS,"rdelta",rdelta
                rbeta=0.0d0
                write(1,*)"rbeta set to 0",rbeta

                LP=LP+rbeta*LP

                write(1,*)"The Rotations (LP) all and non-redundant"
                call print_mat1(1,nij+nC0,LP,1)

                YP00=YP

              end do

              LP=0.0d0
              do i=1,nCG+naddCG
                write(1,*)"The delta Lagrange parameters in (P)CG : "
                call print_mat1(1,nij+nC0,LPx(:,i),1)
                Lp=LP+LPx(:,i)
              end do

            end if

! to Leon :!         printout the Lagrange multipliers for orbital and MPS

            write(1,*)"The Lagrange parameters just after (P)CG : "
            call print_mat(1,nij+nC0,LP,1)
            call flush(1)

            ! Check the results

            allocate(TL1(nij+nC0)); TL1=0.0d0

            call MXV(nij+nC0,nij+nC0,HP,LP,TL1)

            write(1,*)"  "
            write(1,*)"The Hess*Lag"
            call print_mat1(1,nij+nC0,TL1,1)
            write(1,*)"The gradients for specific state"
            call print_mat1(1,nij+nC0, YP,1)
            write(1,*)"   "

            deallocate(TL1)

            write(1,*)"after finish LR solver"
            call flush(1)

          else

          end if

!          LRsign=.false.
!          write(1,*)"The Lagrange parameters : "
!          write(1,*)" ================================================ "
!          write(1,*)" Correct the sigh if RHS(DMRG-LR)=-RHS(Molcas)  "
!          write(1,*)"        Add automatic sign check later  "
!          write(1,*)" ================================================ "
!          write(6,*)" ================================================ "
!          write(6,*)" Correct the sigh if RHS(DMRG-LR)=-RHS(Molcas)  "
!          write(6,*)"        Add automatic sign check later  "
!          write(6,*)" ================================================ "
!          !LRsign=.true.
!          if(LRsign)then ! yma :
!            LP=-LP
!          end if
          call print_mat(1,nij+nC0,LP,1)
          call flush(1)

!          stop "Done the whole LR solving"

          if(.not.(allocated(mat2%deltR)))then
            allocate(mat2%deltR(norb,norb))
            mat2%deltR=0.0d0
          else
            mat2%deltR=0.0d0
          end if

! Lagrange - orbital part
          do i=1,nij
            i1=valid(i,1)
            j1=valid(i,2)
!! recover the anti-symmetric rotations
!          mat2%deltR(i1,j1)=  XP(i+1)/2.0d0
!          mat2%deltR(j1,i1)= -XP(i+1)/2.0d0
            mat2%deltR(i1,j1)=  LP(i)/1.0d0
            mat2%deltR(j1,i1)= -LP(i)/1.0d0
            write(1,*)i,"-th R Lagrange in",i1,j1,mat2%deltR(i1,j1)
            if(dabs(mat2%deltR(i1,j1)).gt.0.10)then
              write(2,*)"notice : large rotation between",i1,j1,&
                        " : ",mat2%deltR(i1,j1)
            end if
          end do

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
                    ! avoid possible over-shooting by 1/2 (no need)
              mat2%deltR(ii,jj)=mat2%deltR(ii,jj)/1.0d0
              mat2%deltR(kk,ll)=mat2%deltR(ii,jj)*redu%sign_mat(kk,ll)
              write(1,*)ii,jj,mat2%deltR(ii,jj),kk,ll,mat2%deltR(kk,ll)
                  end do
                end do
              end if
              js0=js0+orb%total(jsym)
            end do
            is0=is0+orb%total(isym)
          end do


          L_R%R=mat2%deltR
          write(1,*)"The orbital Lagrange : "
          call print_mat(norb,norb,L_R%R,1)
          call flush(1)

! Lagrange - MPSci part
          offset=nij
          ! non-redundant orbital-MPS part
          ici=0
          do iroot=0,dmrg_nstates-1
            !nC=MPS(iroot+1)%nC
            nC=L_R%nC
            do iC=1,nC
              ici=ici+1
              MPS(iroot+1)%Ci(iC)%lv=LP(offset+ici)
              write(1,*)"Davidson & Lagrange MPSci",&
                         MPS(iroot+1)%Ci(iC)%dv,    &
                         MPS(iroot+1)%Ci(iC)%lv
            end do
          end do
          call flush(1)

          !mat2%deltR=0.0d0
         ! 1) Orb modification
            allocate(TM1(norb,norb)); TM1=0.0d0
            allocate(TM2(norb,norb)); TM2=0.0d0
            call CAL_T(norb,mat2%deltR,1.0d0,150,TM1)
            call T_TO_U(norb,TM1,TM2)
            call MXM(norb,mat2%U,TM2,TM1)
            mat2%U=TM1
            mat2%T=mat2%U
            do i=1,norb
              do j=1,norb
                if(dabs(mat2%U(i,j)).lt.1.0d-16)then
                  mat2%U(i,j)=0.0d0
                end if
              end do
              mat2%T(i,i)=mat2%U(i,i)-1.0d0
            end do
            deallocate(TM1,TM2)

! Since the specific state is improved, the linear-equation also need to be relaxed

            write(1,*)"The 1st order T, start"
            call flush(1)

          if(iloop.ge.iloopmax)exit

        end do  !PCG iterations

        deallocate(porb_redu)
        deallocate(sigma)
        deallocate(xvec)
        deallocate(rvec)
        deallocate(V1,V2,V3)

        deallocate(Y,G)
        deallocate(SAvec)
        !write(1,*)"deallocateed Y "
        !call flush(1)
        if(allocated(valid))then
          deallocate(valid)
        end if

        deallocate(HPM)

        deallocate(H)
        deallocate(porb)
        deallocate(LP,YP,HP)

        deallocate(L_R%Y)
        deallocate(L_R%HOO)
        deallocate(L_R%HCO)
        deallocate(L_R%HOC)
        deallocate(L_R%HCC)

        deallocate(conf)

        write(1,*)"deallocateed valid  "
        call flush(1)

      End Subroutine ConDMRGLR


