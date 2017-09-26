      Subroutine PreDMRGLR(iloop)
 
        use global_control
        use matrix
 
        character    ctmp1,ctmp
        character*2  ctmp2
        character*2  ctmp3
        character*192 string0,string1,string2,string3,string4
       
        integer :: iloop

        double precision dv 
! Read-in the Hami with augmented dim       
        double precision,allocatable::Hami_all(:,:)

! For the test case; check the completeness of RDM-deri       
        character*72 oneRDMr
        character*72 twoRDMr
! Whole 1' & 2' RDMs for checking correctness
        double precision oRDMchk(nocc,nocc)
        double precision tRDMchk(nocc,nocc,nocc,nocc)
! Some temparary matrix
        double precision,allocatable::TM1(:,:),TM2(:,:)
        double precision,allocatable::TM3(:,:),TM4(:,:)
        double precision,allocatable::TM5(:,:),TM6(:,:)
    
        character*72 oneRDMfile,twoRDMfile
 
        double precision,allocatable::GM1(:,:,:,:)

        if(.not.allocated(mps))then
          allocate(mps(dmrg_nstates)) 
        end if

        ! Allocate orbital related parameters
        do i=1,dmrg_nstates

          ! Half part of the gradients
          if(.not.allocated(mps(i)%A))then
            allocate(mps(i)%A(norb,norb));   mps(i)%A=0.0d0
          else
            mps(i)%A=0.0d0 
          end if
          ! gradients for specific state 
          if(.not.allocated(mps(i)%G))then
            allocate(mps(i)%G(norb,norb));   mps(i)%G=0.0d0
          else
            mps(i)%G=0.0d0
          end if  
          ! to be solved orbital roations for specfic state
          if(.not.allocated(mps(i)%R))then
            allocate(mps(i)%R(norb,norb));   mps(i)%R=0.0d0
          else
            mps(i)%R=0.0d0
          end if
          ! RDMs
          if(.not.allocated(mps(i)%d))then
            allocate(mps(i)%d(nact,nact));   mps(i)%d=0.0d0
          else
            mps(i)%d=0.0d0
          end if
          if(.not.allocated(mps(i)%p))then
            allocate(mps(i)%p(nact,nact,nact,nact))  
            mps(i)%p=0.0d0          
          else
            mps(i)%p=0.0d0          
          end if

        end do 

        write(6,*)"allocated mps info."
        call flush(6)
          
! Here should be read in the number of davidson elements, 
       ! Should from text file
       ! List all the 1'-RDM-deri, for counting the davidson elements 
        do iroot=0,dmrg_nstates-1

          ! recover the closed shell structure          

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if       

          open(unit=200,file="rdms-deri_stat")
            string1=""
            string2=""
            string3=""
            string4=""
            string1="ls oneRDM."//trim(ctmp2)//".deriL_* > info1.tmp"
            string2="ls twoRDM."//trim(ctmp2)//".deriL_* > info2.tmp"
            string3="ls oneRDM."//trim(ctmp2)//".deriR_* > info3.tmp"
            string4="ls twoRDM."//trim(ctmp2)//".deriR_* > info4.tmp"  
            write(200,*)trim(string1)
            write(200,*)trim(string2)
            write(200,*)trim(string3)
            write(200,*)trim(string4)
            string0="wc -l MPSCi_value."//trim(ctmp2)//&
                    ".info > info_line_all.tmp"
            write(200,*)trim(string0)
          close(200) 
          call system("chmod +x rdms-deri_stat")
          call system("./rdms-deri_stat")
          call system("wc -l info1.tmp > info_line.tmp")

          open(unit=100,file="info_line_all.tmp")
            read(100,*)nC_all
          close(100)
          write(2,*)"Number of all Davidson-vector elems :"
          write(2,*) nC_all
          MPS(iroot+1)%nC_all=nC_all

          open(unit=100,file="info_line.tmp")
            read(100,*)nC
          close(100)
          write(2,*)"Number of Davidson-vector elems to be coupled :"
          write(2,*) nC
          MPS(iroot+1)%nC=nC

          ! Ci is a class for handling coupling
          if(.not.allocated(MPS(iroot+1)%ci))then
            allocate(MPS(iroot+1)%ci(nC))
          end if

          if(.true.)then
            !write(*,*)"!!! > Readin cofficient from MPSCi.x.info"
            string1=""
            string1="MPSCi."//trim(ctmp2)//".info"
            iC=0
            open(unit=100,file=trim(string1))
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
                      MPS(iroot+1)%Ci(iC)%dv=dv
                      MPS(iroot+1)%Ci(iC)%num=inum
                    end if
                  end do
                end do
              end do
!              nC_all=inum
            close(100)
          else
            string1=""
            string1="MPSCi_value."//trim(ctmp2)//".info"
            open(unit=100,file=trim(string1))
              inum=0
              do iC=1,nC_all
                read(100,*)dv,i,j,k
                if(dabs(dv).lt.1.0d-12)then
                  !Ci(iC)%dv=0.0d0
                else
                  inum=inum+1 
                  MPS(iroot+1)%Ci(inum)%dv=dv
                  MPS(iroot+1)%Ci(inum)%num=iC
                end if
              end do
            close(100)
          end if

          write(2,*)"Initial Davidson vectors"
          call print_mat(1,nc,MPS(iroot+1)%Ci%dv,2)
          write(2,*)"     mapping index  "
          write(2,*)MPS(iroot+1)%Ci%num
          call flush(2) 

        ! read in the sign info. from overlap 
          string1=""
          string1="Loverlap."//trim(ctmp2)//".txt"
          open(unit=100,file=trim(string1))
            do iC=1,nC
              read(100,*)ctmp,ctmp,ctmp,ctmp, dv
              read(100,*)
              MPS(iroot+1)%Ci(iC)%dsL = sign(1.0d0,dv)
            end do
          close(100)
          string1=""
          string1="Roverlap."//trim(ctmp2)//".txt"
          open(unit=100,file=trim(string1))
            do iC=1,nC
              read(100,*)ctmp,ctmp,ctmp,ctmp, dv
              read(100,*)
              MPS(iroot+1)%Ci(iC)%dsR = sign(1.0d0,dv)
              !write(2,*)MPS(iroot+1)%Ci(iC)%dsR
            end do
          close(100)

!          write(2,*)"Sign correct vector L"
!          call print_mat(1,nc,MPS(iroot+1)%Ci%dsL,2)
!          call flush(2) 
!          write(2,*)"Sign correct vector R"
!          call print_mat(1,nc,MPS(iroot+1)%Ci%dsR,2)
!          call flush(2) 
!          write(2,*)"Attention!! Unified the Sign L and R using L"
!          write(2,*)"Attention!! Unified the Sign L and R using L"
!          MPS(iroot+1)%Ci%dsR = MPS(iroot+1)%Ci%dsL

! For read in 1-RDM-deri (<i|o|mps> part) 
          open(unit=100,file="info1.tmp")
            do iC=1,nC
              oneRDMr=""
! symmetrized 1-RDM-deri
              allocate(MPS(iroot+1)%Ci(iC)%D(nact,nact))
              MPS(iroot+1)%Ci(iC)%D=0.0d0
! Left side part
              allocate(MPS(iroot+1)%Ci(iC)%DL(nact,nact))
              MPS(iroot+1)%Ci(iC)%DL=0.0d0
              read(100,*,iostat=ierr)oneRDMr
              if(ierr.ne.0)exit
              open(unit=101,file=oneRDMr)
                read(101,*)ij0
                do i=1,nact
                  do j=1,nact
                    read(101,*)ij,ji,dv
                    ! if NAN set to 0
                    if (isnan(dv).or.dabs(dv).lt.dthrs)then
                      dv = 0.0d0
                    end if
                    MPS(iroot+1)%Ci(iC)%DL(iorder(ij+1),iorder(ji+1))=dv
                  end do
                end do
              close(101)

!              write(2,*)"1-RDM-deri L for element",iC
!              call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%DL,2)
!              call flush(2)
 
            end do
          close(100)
          write(2,*)"Finish importing 1-RDM-deri L for state",iroot
          call flush(2)

! For read in 1-RDM-deri (<mps|o|i> part) 
          open(unit=100,file="info3.tmp")
            do iC=1,nC
              oneRDMr=""
! Right side part
              allocate(MPS(iroot+1)%Ci(iC)%DR(nact,nact))
              MPS(iroot+1)%Ci(iC)%DR=0.0d0
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
                    MPS(iroot+1)%Ci(iC)%DR(iorder(ij+1),iorder(ji+1))=dv
                  end do
                end do
              close(101)

!              write(2,*)"1-RDM-deri R for element",iC
!              call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%DR,2)
!              call flush(2) 

            end do
          close(100)
          write(2,*)"Finish importing 1-RDM-deri R for state",iroot
          call flush(2)


! For read in 2-RDM-deri (L)
          open(unit=100,file="info2.tmp")
            do iC=1,nC
              twoRDMr=""
! Avergaed 
              allocate(MPS(iroot+1)%Ci(iC)%P(nact,nact,nact,nact))
              MPS(iroot+1)%Ci(iC)%P=0.0d0
! Left side part
              allocate(MPS(iroot+1)%Ci(iC)%PL(nact,nact,nact,nact))
              MPS(iroot+1)%Ci(iC)%PL=0.0d0
              read(100,*,iostat=ierr)twoRDMr
              if(ierr.ne.0)exit
              open(unit=101,file=twoRDMr)
                read(101,*)ijkl
                ijkl=0  ! Just a counter
                do
                  read(101,*,iostat=ierr)ij,jk,kl,li,dv
                  if(ierr.ne.0)exit
                  ijkl=ijkl+1
                  MPS(iroot+1)%Ci(iC)%PL(&
                            iorder(ij+1),&
                            iorder(jk+1),&
                            iorder(kl+1),&
                            iorder(li+1))=dv
                  MPS(iroot+1)%Ci(iC)%PL(&
                            iorder(kl+1),&
                            iorder(li+1),&
                            iorder(ij+1),&
                            iorder(jk+1))=dv
                  MPS(iroot+1)%Ci(iC)%PL(&
                            iorder(jk+1),&
                            iorder(ij+1),&
                            iorder(li+1),&
                            iorder(kl+1))=dv
                  MPS(iroot+1)%Ci(iC)%PL(&
                            iorder(li+1),&
                            iorder(kl+1),&
                            iorder(jk+1),&
                            iorder(ij+1))=dv
                end do
              close(101)

              write(322,*)"2-RDM-deri (L) for iC",iC
              call print_gat &
              (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PL,322)

              MPS(iroot+1)%Ci(iC)%PL=MPS(iroot+1)%Ci(iC)%PL*2.0d0
            end do
          close(100)
          write(2,*)"Finish importing 2-RDM-deri L for state",iroot
          call flush(2)

! For read in 2-RDM-deri (R)
          open(unit=100,file="info4.tmp")
            do iC=1,nC
              twoRDMr=""
              allocate(MPS(iroot+1)%Ci(iC)%PR(nact,nact,nact,nact))
              MPS(iroot+1)%Ci(iC)%PR=0.0d0
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
                  MPS(iroot+1)%Ci(iC)%PR(&
                            iorder(ij+1),&
                            iorder(jk+1),&
                            iorder(kl+1),&
                            iorder(li+1))=dv
                  MPS(iroot+1)%Ci(iC)%PR(&
                            iorder(kl+1),&
                            iorder(li+1),&
                            iorder(ij+1),&
                            iorder(jk+1))=dv
                  MPS(iroot+1)%Ci(iC)%PR(&
                            iorder(jk+1),&
                            iorder(ij+1),&
                            iorder(li+1),&
                            iorder(kl+1))=dv
                  MPS(iroot+1)%Ci(iC)%PR(&
                            iorder(li+1),&
                            iorder(kl+1),&
                            iorder(jk+1),&
                            iorder(ij+1))=dv
                end do
              close(101)

              write(324,*)"2-RDM-deri (R) for iC",iC
              call print_gat &
              (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PR,324)

              MPS(iroot+1)%Ci(iC)%PR=MPS(iroot+1)%Ci(iC)%PR*2.0d0
            end do
          close(100)
          write(2,*)"Finish importing 2-RDM-deri R for state",iroot
          call flush(2)

! Closed shell update for all RDM-deri 

          do iC=1,nC

            call CS_dim_fro2cls()
            write(*,*)"Start the ", iC,"iterations for RDM-deri"

            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
!            call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%dL,6)
!            call print_gat(nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%p,6)

            ! L
! Fix the sign for RDM-deri (closed shell part)
              if(MPS(iroot+1)%Ci(iC)%dv.eq.&
              dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsL)then
                dv= 1.0d0*MPS(iroot+1)%Ci(iC)%dv
              else
                dv=-1.0d0*MPS(iroot+1)%Ci(iC)%dv
              end if
            call closed_shell_update&
            (MPS(iroot+1)%Ci(iC)%dL,MPS(iroot+1)%Ci(iC)%PL,TM1,GM1,dv)  ! CCCC
            deallocate(MPS(iroot+1)%Ci(iC)%dL,MPS(iroot+1)%Ci(iC)%PL)
            allocate(MPS(iroot+1)%Ci(iC)%dL(nocc,nocc))
            allocate(MPS(iroot+1)%Ci(iC)%PL(nocc,nocc,nocc,nocc))
            MPS(iroot+1)%Ci(iC)%dL=TM1
            MPS(iroot+1)%Ci(iC)%pL=GM1
            deallocate(TM1,GM1)
!           call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%dL,6)
!           call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%pL,6)

          end do

          do iC=1,nC

            call CS_dim_fro2cls()
            write(*,*)"Start the ", iC,"iterations for RDM-deri"

            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
!            call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%dR,6)
!            call print_gat(nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%p,6)

            ! R
! Fix the sign for RDM-deri (closed shell part)
              if(MPS(iroot+1)%Ci(iC)%dv.eq.&
              dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsR)then
                dv= 1.0d0*MPS(iroot+1)%Ci(iC)%dv
              else
                dv=-1.0d0*MPS(iroot+1)%Ci(iC)%dv
              end if

            call closed_shell_update&
            (MPS(iroot+1)%Ci(iC)%dR,MPS(iroot+1)%Ci(iC)%PR,TM1,GM1,dv)  ! CCCC
            deallocate(MPS(iroot+1)%Ci(iC)%dR,MPS(iroot+1)%Ci(iC)%PR)
            allocate(MPS(iroot+1)%Ci(iC)%dR(nocc,nocc))
            allocate(MPS(iroot+1)%Ci(iC)%PR(nocc,nocc,nocc,nocc))
            MPS(iroot+1)%Ci(iC)%dR=TM1
            MPS(iroot+1)%Ci(iC)%pR=GM1
            deallocate(TM1,GM1)

          end do
          call CS_dim_fro2cls()

!          stop

!          call print_mat(nact,nact,oRDMchk,2)
!          call flush(2) 

! Testing for the correctness
          oRDMchk=0.0d0
          tRDMchk=0.0d0
          write(2,*)"Before testing for correctness"
          call flush(2)
          do iC=1,nC
            write(2,*)"Davidson C_",iC," : ", MPS(iroot+1)%Ci(iC)%dv, &
                       MPS(iroot+1)%Ci(iC)%dsL, MPS(iroot+1)%Ci(iC)%dsR
!                 dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsL, &
!                 dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsR
!            call flush(2)

            if(.true.)then
              ! Fix the sign for RDM-deri
              if(MPS(iroot+1)%Ci(iC)%dv.eq.&
              dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsL)then
              else
                MPS(iroot+1)%Ci(iC)%DL=-1.0d0*MPS(iroot+1)%Ci(iC)%DL
                MPS(iroot+1)%Ci(iC)%PL=-1.0d0*MPS(iroot+1)%Ci(iC)%PL
              end if
              if(MPS(iroot+1)%Ci(iC)%dv.eq.&
              dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsR)then
              else
                MPS(iroot+1)%Ci(iC)%DR=-1.0d0*MPS(iroot+1)%Ci(iC)%DR
                MPS(iroot+1)%Ci(iC)%PR=-1.0d0*MPS(iroot+1)%Ci(iC)%PR
              end if
            else
                MPS(iroot+1)%Ci(iC)%DL=MPS(iroot+1)%Ci(iC)%dsL*MPS(iroot+1)%Ci(iC)%DL
                MPS(iroot+1)%Ci(iC)%PL=MPS(iroot+1)%Ci(iC)%dsL*MPS(iroot+1)%Ci(iC)%PL
                MPS(iroot+1)%Ci(iC)%DR=MPS(iroot+1)%Ci(iC)%dsR*MPS(iroot+1)%Ci(iC)%DR
                MPS(iroot+1)%Ci(iC)%PR=MPS(iroot+1)%Ci(iC)%dsR*MPS(iroot+1)%Ci(iC)%PR
            end if

            write(2,*)"iC DL"
            call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%DL,2)
            write(2,*)"iC DR"
            call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%DR,2)

! Sum up L and R, then average
            MPS(iroot+1)%Ci(iC)%D=MPS(iroot+1)%Ci(iC)%DL&
                                 +MPS(iroot+1)%Ci(iC)%DR
            MPS(iroot+1)%Ci(iC)%P=MPS(iroot+1)%Ci(iC)%PL&
                                 +MPS(iroot+1)%Ci(iC)%PR
            MPS(iroot+1)%Ci(iC)%D=MPS(iroot+1)%Ci(iC)%D/2.0d0
            MPS(iroot+1)%Ci(iC)%P=MPS(iroot+1)%Ci(iC)%P/2.0d0

            oRDMchk=oRDMchk+MPS(iroot+1)%Ci(iC)%D*MPS(iroot+1)%Ci(iC)%dv
            tRDMchk=tRDMchk+MPS(iroot+1)%Ci(iC)%P*MPS(iroot+1)%Ci(iC)%dv
          end do
 
          
          write(2,*)"resembled oRDM "
          call print_mat(nocc,nocc,oRDMchk,2)

          write(2,*)" !! Corrected RDM-deri (L & R) is used. !!"
          call flush(2)


!  For read-in normal RDM for each state

! 1-RDMs
          if(.not.allocated(mat2%d))then
             allocate(mat2%d(nact,nact))           ; mat2%d=0.0d0
          else
             mat2%d=0.0d0
          end if
         !allocate(TM1(nact,nact))                 ;    TM1=0.0d0

! 2-RDMs
          if(.not.allocated(mat2%p))then
            allocate(mat2%p(nact,nact,nact,nact)) ; mat2%p=0.0d0
          else
             mat2%p=0.0d0
          end if
         !allocate(GM1(nact,nact,nact,nact))       ;    GM1=0.0d0

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
          MPS(iroot+1)%D=mat2%d

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
          MPS(iroot+1)%P=mat2%p*2.0d0

          
            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
            !call print_mat(nact,nact,MPS(iroot+1)%d,6)
            !call print_gat(nact,nact,nact,nact,MPS(iroot+1)%p,6)
            call closed_shell_update&
            (MPS(iroot+1)%d,MPS(iroot+1)%P,TM1,GM1,1.0d0)  ! CCCC
            deallocate(MPS(iroot+1)%d,MPS(iroot+1)%P)
            allocate(MPS(iroot+1)%d(nocc,nocc))
            allocate(MPS(iroot+1)%P(nocc,nocc,nocc,nocc))
            MPS(iroot+1)%d=TM1
            MPS(iroot+1)%p=GM1
            !call print_mat(nocc,nocc,MPS(iroot+1)%d,6)
            !call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%p,6)
            deallocate(TM1,GM1)

! Check the energy for each state
          call energy_gen(orb%nsub,orb%act,orb%total,nact,norb, &
                             T,U,MPS(iroot+1)%D,MPS(iroot+1)%P, &
                                 MPS(iroot+1)%RDM_energy)
 
          write(2,*)"DMRG(RDM) energy for state",iroot
          write(2,*) MPS(iroot+1)%RDM_energy
          write(2,*) MPS(iroot+1)%RDM_energy+FRONRE0


! Check the MPSci gradients for each Davidson vectors
  
          do iC=1,nC

!            call CS_dim_fro2cls()
            write(*,*)"Start the ", iC,"iterations" 

            ! Closed shell update
!            allocate(TM1(nocc,nocc))
!            allocate(GM1(nocc,nocc,nocc,nocc))
!            call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%d,6)
!            call print_gat(nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%p,6)

!            dv=MPS(iroot+1)%Ci(iC)%dv
!            call closed_shell_update&
!            (MPS(iroot+1)%Ci(iC)%d,MPS(iroot+1)%Ci(iC)%P,TM1,GM1,dv)  ! CCCC
!            deallocate(MPS(iroot+1)%Ci(iC)%d,MPS(iroot+1)%Ci(iC)%P)
!            allocate(MPS(iroot+1)%Ci(iC)%d(nocc,nocc))
!            allocate(MPS(iroot+1)%Ci(iC)%P(nocc,nocc,nocc,nocc))
!            MPS(iroot+1)%Ci(iC)%d=TM1
!            MPS(iroot+1)%Ci(iC)%p=GM1
!            call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%d,6)
!            call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%p,6)
!            deallocate(TM1,GM1)

          ! Gradient that calculated by energy_gen
          call energy_gen(orb%nsub,orb%act,orb%total,nact,norb,T,U,  &
                       MPS(iroot+1)%Ci(iC)%D,MPS(iroot+1)%Ci(iC)%P,  &
                       MPS(iroot+1)%Ci(iC)%grad)

            write(2,*)&
          "Energy component(mpsci, ref) & gradient for Davidson vectors",iC
            write(2,*) MPS(iroot+1)%Ci(iC)%grad, &
                       MPS(iroot+1)%RDM_energy*MPS(iroot+1)%Ci(iC)%dv, &
                       MPS(iroot+1)%Ci(iC)%grad &
                      -MPS(iroot+1)%RDM_energy*MPS(iroot+1)%Ci(iC)%dv

          end do

          string1=""
          string1="oneRDM."//trim(ctmp2)//".reformed"
          open(unit=100,file=trim(string1))
            do i=1,nocc
              do j=1,nocc
                write(100,*)i-1,j-1,oRDMchk(i,j)
              end do
            end do
          close(100)
          MPS(iroot+1)%D=oRDMchk

          write(2,*)" Refactored one-RDM "        
          call print_mat(nocc,nocc,MPS(iroot+1)%D,2)
          call flush(2)            

          string1=""
          string1="twoRDM."//trim(ctmp2)//".reformed"
          open(unit=100,file=trim(string1))
            do i=1,nocc
              do j=1,nocc
                do k=1,nocc
                  do l=1,nocc
                    write(100,*)i-1,j-1,k-1,l-1,tRDMchk(i,j,k,l)/2.0d0
                  end do
                end do
              end do
            end do
          close(100)
          MPS(iroot+1)%P=tRDMchk

          write(6,*)"The reformed two-RDMs for checking, fort.18"
!           
!          call print_mat(nocc,nocc,MPS(iroot+1)%d,6)
          call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%p,18)

!          stop

! Allocate the orb-CI and CI-CI Hessian
          allocate(MPS(iroot+1)%HCR(nC,norb,norb)); MPS(iroot+1)%HCR=0.0d0
          allocate(MPS(iroot+1)%KCR(nC,norb,norb)); MPS(iroot+1)%KCR=0.0d0
          allocate(MPS(iroot+1)%HCC(nC,nC));        MPS(iroot+1)%HCC=0.0d0
! Allocate for the full Hamiltonian
          allocate(MPS(iroot+1)%Hami(nC,nC)) ;      MPS(iroot+1)%Hami=0.0d0
          allocate(MPS(iroot+1)%Hvec(nC,nC)) ;      MPS(iroot+1)%Hvec=0.0d0
          allocate(MPS(iroot+1)%Hval(nC))    ;      MPS(iroot+1)%Hval=0.0d0

! Here should be the read in the full Hamiltonian
          ! read in code 
          allocate(Hami_all(nC_all,nC_all))
          Hami_all=0.0d0
          string1=""
          string1="Hami."//trim(ctmp2)//".txt"
          open(unit=100,file=trim(string1))
            do iC=1,nC_all
              read(100,*)(Hami_all(iC,j),j=1,nC_all)
            end do
          close(100)
          do iC=1,nC
            do jC=1,nC
                       MPS(iroot+1)%Hami(iC,jC)= &
              Hami_all(MPS(iroot+1)%Ci(iC)%num,MPS(iroot+1)%Ci(jC)%num)
            end do             
          end do

          write(13,*)"Hamiltonian Matrix all "        
          call print_mat(nC_all,nC_all,Hami_all,13)
          call flush(13)   

          write(13,*)"Hamiltonian Matrix "        
          call print_mat(nC,nC,MPS(iroot+1)%Hami,13)
          call flush(13)   
          deallocate(Hami_all)

! Set very small value to 0
          do i=1,nC
            do j=1,nC
              if(dabs(MPS(iroot+1)%Hami(i,j)).lt.1.0d-12)then
                MPS(iroot+1)%Hami(i,j)=0.0d0
              end if
            end do
          end do
          write(2,*)"Hamiltonian Matrix "        
          call print_mat(nC,nC,MPS(iroot+1)%Hami,2)
          call flush(2) 

! Re-calculate the RDM_energies base on the refactored RDMs

          call energy_gen(orb%nsub,orb%act,orb%total,nact,norb,T,U,  &
                                     MPS(iroot+1)%D,MPS(iroot+1)%P,  &
                                     MPS(iroot+1)%RDM_energy)

          write(2,*)"The REM-refactored E (no-NRE and with-NRE) is", &
                    MPS(iroot+1)%RDM_energy, MPS(iroot+1)%RDM_energy+FRONRE0

            do iC=1,nC
              do jC=1,nC
                MPS(iroot+1)%HCC(iC,jC)=MPS(iroot+1)%Hami(iC,jC)
                if(.false.)then

                  write(3,*)"The SA coupling in c-c is also considered" 

             write(1,*)"SP_E, SA_E",MPS(iroot+1)%RDM_energy,RDM_energy

                  dv=MPS(iroot+1)%Ci(iC)%dv*MPS(iroot+1)%Ci(jC)%dv &
                    *(MPS(iroot+1)%RDM_energy-RDM_energy)          &
                     *dmrg_weight(iroot+1)*dmrg_weight(iroot+1) 

             write(1,*)"newly term", dv

                  MPS(iroot+1)%HCC(iC,jC)=MPS(iroot+1)%HCC(iC,jC)+dv

                end if
              end do

              write(6,*)" FRONRE,FRONRE0,FRONREA",FRONRE,FRONRE0,FRONREA

              write(6,*)" = = = = = = = = = = "
              write(6,*)" iroot, iC, jC ",iroot,iC,jC
              write(6,*)MPS(iroot+1)%HCC(iC,iC),MPS(iroot+1)%RDM_energy

              MPS(iroot+1)%HCC(iC,iC)=MPS(iroot+1)%HCC(iC,iC)&
                                     +FRONREA                &
                                     -MPS(iroot+1)%RDM_energy&
                                     -FRONRE0 

              write(6,*)"  updated HCC "

              write(6,*) MPS(iroot+1)%HCC(iC,iC)
              write(6,*)" = = = = = = = = = = "

            end do

          MPS(iroot+1)%HCC=MPS(iroot+1)%HCC*2.0d0          
          MPS(iroot+1)%HCC=MPS(iroot+1)%HCC*dmrg_weight(iroot+1)

       ! Construct the A(Fock) matrix for each state
        call fock_gen(orb%nsub,orb%act,orb%total,     &
                      nact,norb,T,U,MPS(iroot+1)%D,   &
                      MPS(iroot+1)%P,MPS(iroot+1)%A)

          write(2,*)"The A-matrix for state",iroot
          call print_mat(norb,norb,MPS(iroot+1)%A,2)

          do i=1,norb
            do j=1,norb
              MPS(iroot+1)%G(i,j) = &
              MPS(iroot+1)%A(i,j)-MPS(iroot+1)%A(j,i)
            end do
          end do

          ! Construct the A^I matrix for every Davidson elemments for each state
 
          do iC=1,nC
            allocate(MPS(iroot+1)%ci(iC)%A(norb,norb))
            allocate(MPS(iroot+1)%ci(iC)%AL(norb,norb))
            allocate(MPS(iroot+1)%ci(iC)%AR(norb,norb))
            MPS(iroot+1)%ci(iC)%A=0.0d0
            MPS(iroot+1)%ci(iC)%AL=0.0d0
            MPS(iroot+1)%ci(iC)%AR=0.0d0

            call fock_gen(orb%nsub,orb%act,orb%total,     &
                   nact,norb,T,U,MPS(iroot+1)%Ci(iC)%D,   &
                   MPS(iroot+1)%Ci(iC)%P,MPS(iroot+1)%Ci(iC)%A)

            write(2,*)"The A^I-matrix for Davidson elem",iC
            call print_mat(norb,norb,MPS(iroot+1)%Ci(iC)%A,2)

          end do  ! loop the MPS_ci elements

! Construct the Coupling part of the coupled Hessian (only for single state)
!          do iC=1,nC
!            do i=1,norb
!              do j=1,norb
!                MPS(iroot+1)%KCR(iC,i,j)=MPS(iroot+1)%Ci(iC)%A(i,j)&
!                           -MPS(iroot+1)%Ci(iC)%A(j,i) ! For property check

!                MPS(iroot+1)%HCR(iC,i,j)=MPS(iroot+1)%Ci(iC)%A(i,j)&
!                                        -MPS(iroot+1)%Ci(iC)%A(j,i)&
!                                -MPS(iroot+1)%Ci(iC)%dv*mat2%A(i,j)&
!                                +MPS(iroot+1)%Ci(iC)%dv*mat2%A(j,i) ! <-- second line ()

!               write(*,*)" more than two states also need TDM"

!              end do
!            end do
!          end do
!          HCR=1.0d0*HCR*2.0d0

          ! Recover the old closed shell structure
          call CS_dim_fro2cls()

        end do  ! Finish read-in all properties for every states

        write(6,*)" Finish read-in all properties for every states"
        call flush(6)

        ! Calculate the preconditioner for MPS part
        call SA_prec_mps() 

! Read-in state-transfer properties
        allocate(MPSTDMs(dmrg_nstates,dmrg_nstates))
     
        do iroot=0,dmrg_nstates-1
          ir=iroot+1
          do jroot=iroot+1,dmrg_nstates-1
            jr=jroot+1

            ctmp2=""
            if(iroot.ge.10)then
              write(ctmp2,"(I2)")iroot
            else
              write(ctmp1,"(I1)")iroot
              ctmp2(1:1)=ctmp1
            end if

            ctmp3=""
            if(jroot.ge.10)then
              write(ctmp3,"(I2)")jroot
            else
              write(ctmp1,"(I1)")jroot
              ctmp3(1:1)=ctmp1
            end if

            ! one-body TDMs (pair)
            allocate(MPSTDMs(ir,jr)%D(nact,nact)) 
            MPSTDMs(ir,jr)%D=0.0d0
            allocate(MPSTDMs(jr,ir)%D(nact,nact))
            MPSTDMs(jr,ir)%D=0.0d0
            ! two-body TDMs (pair)
            allocate(MPSTDMs(ir,jr)%P(nact,nact,nact,nact)) 
            MPSTDMs(ir,jr)%P=0.0d0
            allocate(MPSTDMs(jr,ir)%P(nact,nact,nact,nact)) 
            MPSTDMs(jr,ir)%P=0.0d0
           
            ! oneRDM.iroot.jroot 
            string1=""
            string1="oneRDM."//trim(ctmp2)//"."//trim(ctmp3)
            open(unit=101,file=trim(string1))
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  MPSTDMs(ir,jr)%D(i,j)=dv
                end do
              end do
            close(101)
            write(2,*)"oneTDM ",iroot,jroot
            call print_mat(nact,nact,MPSTDMs(ir,jr)%D,2)
            call flush(2) 

            ! oneRDM.jroot.iroot 
            string1=""
            string1="oneRDM."//trim(ctmp3)//"."//trim(ctmp2)
            open(unit=101,file=trim(string1))
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  MPSTDMs(jr,ir)%D(i,j)=dv
                end do
              end do
            close(101)             
            write(2,*)"oneTDM ",jroot,iroot
            call print_mat(nact,nact,MPSTDMs(jr,ir)%D,2)
            call flush(2) 

            ! twoRDM.iroot.jroot 
            string1=""
            string1="twoRDM."//trim(ctmp2)//"."//trim(ctmp3)
            open(unit=101,file=trim(string1))
              read(101,*)ijkl
              ijkl=0  ! Just a counter
              do
                read(101,*,iostat=ierr)ij,jk,kl,li,dv
                if(ierr.ne.0)exit
                ijkl=ijkl+1
                MPSTDMs(ir,jr)%P(      &
                          iorder(ij+1),&
                          iorder(jk+1),&
                          iorder(kl+1),&
                          iorder(li+1))=dv
                MPSTDMs(ir,jr)%P(      &
                          iorder(kl+1),&
                          iorder(li+1),&
                          iorder(ij+1),&
                          iorder(jk+1))=dv
                MPSTDMs(ir,jr)%P(      &
                          iorder(jk+1),&
                          iorder(ij+1),&
                          iorder(li+1),&
                          iorder(kl+1))=dv
                MPSTDMs(ir,jr)%P(      &
                          iorder(li+1),&
                          iorder(kl+1),&
                          iorder(jk+1),&
                          iorder(ij+1))=dv
              end do
            close(101)
            MPSTDMs(ir,jr)%P=MPSTDMs(ir,jr)%P*2.0d0

            ! twoRDM.jroot.iroot 
            string1=""
            string1="twoRDM."//trim(ctmp3)//"."//trim(ctmp2)
            open(unit=101,file=trim(string1))
              read(101,*)ijkl
              ijkl=0  ! Just a counter
              do
                read(101,*,iostat=ierr)ij,jk,kl,li,dv
                if(ierr.ne.0)exit
                ijkl=ijkl+1
                MPSTDMs(jr,ir)%P(      &
                          iorder(ij+1),&
                          iorder(jk+1),&
                          iorder(kl+1),&
                          iorder(li+1))=dv
                MPSTDMs(jr,ir)%P(      &
                          iorder(kl+1),&
                          iorder(li+1),&
                          iorder(ij+1),&
                          iorder(jk+1))=dv
                MPSTDMs(jr,ir)%P(      &
                          iorder(jk+1),&
                          iorder(ij+1),&
                          iorder(li+1),&
                          iorder(kl+1))=dv
                MPSTDMs(jr,ir)%P(      &
                          iorder(li+1),&
                          iorder(kl+1),&
                          iorder(jk+1),&
                          iorder(ij+1))=dv
              end do
            close(101)
            MPSTDMs(jr,ir)%P=MPSTDMs(jr,ir)%P*2.0d0

            write(2,*)"Finished all read-in TDMs "
            call flush(2)

            ! Symmetrilize TDMs (not 100% sure if need, PS: in derivative need)
            MPSTDMs(ir,jr)%D=(MPSTDMs(ir,jr)%D+MPSTDMs(jr,ir)%D)/2.0d0
            MPSTDMs(jr,ir)%D= MPSTDMs(ir,jr)%D
            
            MPSTDMs(ir,jr)%P=(MPSTDMs(ir,jr)%P+MPSTDMs(jr,ir)%P)/2.0d0
            MPSTDMs(jr,ir)%P= MPSTDMs(ir,jr)%P

            ! A-matrix can be constructed using a specific subroutine, then it can reduce many lines
            ! Cnstruct the A (gradients contribute) for transition
            ! one-body part
            allocate(MPSTDMs(ir,jr)%A(norb,norb))
            allocate(MPSTDMs(jr,ir)%A(norb,norb))
            MPSTDMs(ir,jr)%A=0.0d0
            MPSTDMs(jr,ir)%A=0.0d0

            call fock_gen(orb%nsub,orb%act,orb%total,    &
                      nact,norb,T,U,MPSTDMs(ir,jr)%D,    &
                      MPSTDMs(ir,jr)%P,MPSTDMs(ir,jr)%A)

            MPSTDMs(jr,ir)%A=MPSTDMs(ir,jr)%A

!            write(2,*)"Gradient contribution for states transfering"
!            call print_mat(norb,norb,MPSTDMs(ir,jr)%A,2)
!            call flush(2) 

          end do
          MPSTDMs(ir,ir)%D=MPS(iroot+1)%D
          MPSTDMs(ir,ir)%P=MPS(iroot+1)%P
          MPSTDMs(ir,ir)%A=MPS(iroot+1)%A

!          write(2,*)"Finish the PreDMRGLR part"
!          call flush(2)        

        end do

        write(2,*)"Finish the PreDMRGLR part"
        call flush(2)        

!        stop  
       
      End Subroutine PreDMRGLR


