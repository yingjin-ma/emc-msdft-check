      Subroutine PreDMRGLR(iloop)

        use global_control
        use matrix

        character    ctmp1,ctmp
        character*2  ctmp2
        character*2  ctmp3
        character*192 string0,string1,string2,string3,string4

        integer :: iloop

        double precision dv,dv1,dv2,Esa,dvPhiPhj
! Read-in the Hami with augmented dim
        double precision,allocatable::Hami_all(:,:)

! For the test case; check the completeness of RDM-deri
        character*72 oneRDMr
        character*72 twoRDMr
! Whole 1' & 2' RDMs for checking correctness
        double precision,allocatable::oRDMchk(:,:)
        double precision,allocatable::tRDMchk(:,:,:,:)
! Some temparary matrix
        double precision,allocatable::TM1(:,:),TM2(:,:)
        double precision,allocatable::TM3(:,:),TM4(:,:)
        double precision,allocatable::TM5(:,:),TM6(:,:)

        character*72 oneRDMfile,twoRDMfile

        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
!        double precision thresMPSHami 
!        double precision signTDM 

        if(.not.allocated(mps))then
          allocate(mps(dmrg_nstates))
        end if

        if(.not.allocated(mpsall))then
          allocate(mpsall(dmrg_nstates))
        end if

        allocate(oRDMchk(nocc,nocc))
        oRDMchk=0.0d0
        allocate(tRDMchk(nocc,nocc,nocc,nocc))
        tRDMchk=0.0d0

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

! ===================================================================
! to Leon :
!          1) read in the coefficient for the MPS (Davidson vector)
!          2) read in the RDM-deris/TDMs for each Davidson vector
!          3) correct the signs for RDM-deris
! ===================================================================


! Here should be read in the number of Davidson elements,
       ! Should from text file
       ! List all the 1'-RDM-deri, for counting the davidson elements
        do iroot=0,dmrg_nstates-1

          ! recover the closed shell structure
!          write(6,*)"orb%nsub,orb%act,orb%total,nact,norb"
!          write(6,*)orb%nsub,orb%act,orb%total,nact,norb  

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

          call system("cat info_line_all.tmp")
          open(unit=100,file="info_line_all.tmp")
            read(100,*)nC_all
          close(100)
          write(2,*)"Number of Davidson-vector elems for state",iroot
          write(2,*) nC_all
          MPS(iroot+1)%nC_all=nC_all

          if(.false.)then
            open(unit=100,file="info_line.tmp")
              read(100,*)nC
            close(100)
            write(2,*)"Number of Davidson-vector elems to be coupled :"
            write(2,*) nC
            MPS(iroot+1)%nC=nC
          else  
            ! record the effective MPS-elem. position
            call system("python pickup.py info1.tmp")
            call system("mv pickup.vec pickup."//trim(ctmp2))
            open(unit=100,file="pickup."//trim(ctmp2)) 
            nC=0
            do 
              read(100,*,iostat=ierr)
              if(ierr.ne.0)exit
              nC=nC+1 
            end do
            close(100) 
            MPS(iroot+1)%nC_list=nC
            write(6,*)"nC",MPS(iroot+1)%nC_list          
!            if(.not.allocated(MPS(iroot+1)%position))then
            allocate(MPS(iroot+1)%position(nC)) 
            allocate(MPS(iroot+1)%pos(nC,5)) 
            MPS(iroot+1)%pos=-100
            MPS(iroot+1)%position(:)="" 
!            end if
            open(unit=100,file="pickup."//trim(ctmp2)) 
            iC=1
            do 
              read(100,*,iostat=ierr)MPS(iroot+1)%position(iC)
              if(ierr.ne.0)exit
              iC=iC+1       
            end do
            close(100)            
            call system("python pickup_eff.py pickup."//trim(ctmp2))
            open(unit=100,file='pickup.pos')
            do i=1,MPS(iroot+1)%nC_list
              read(100,*)MPS(iroot+1)%pos(i,1),MPS(iroot+1)%pos(i,2),&
                         MPS(iroot+1)%pos(i,3)
              MPS(iroot+1)%pos(i,4)=MPS(iroot+1)%pos(i,1)*1000000+&
                                    MPS(iroot+1)%pos(i,2)*1000+&
                                    MPS(iroot+1)%pos(i,3) 
            end do
            close(100)
          end if

        end do

        ieff=0 
        do iroot=0,dmrg_nstates-2
          do jroot=iroot+1,dmrg_nstates-1
            do iC=1,MPS(iroot+1)%nC_list
              string1=""
              string1=trim(MPS(iroot+1)%position(iC))  
              string2=""
              string2=trim(MPS(jroot+1)%position(iC)) 
              if(trim(string1).eq.trim(string2))then
              else
                ieff=-1
                write(6,*)"At least 1 different MPS-elem. is found" 
                write(6,*)trim(string1)," state-",iroot+1,".",iC
                write(6,*)trim(string2)," state-",jroot+1,".",iC
                exit 
              end if
!               write(6,*)trim(MPS(iroot+1)%position(iC))
!               call flush(6) 
            end do
          end do   
        end do

! ===================================================================
!             Record all the effective MPSci vector  
! ===================================================================

        allocate(mps_head)
        mps_tail => mps_head
        mps_tail%p => null()
        mps_tail%string=trim(MPS(1)%position(1))

        ! unified MPS vector : record the first one
        do iC=2,MPS(1)%nC_list
!          write(6,*)iC,trim(MPS(1)%position(iC))          
          allocate(mps_tail%p)
          mps_tail => mps_tail%p
          mps_tail%p => null()
          mps_tail%string = trim(MPS(1)%position(iC))
        end do
!        mps_tail=> mps_ptr

!        mps_ptr%p => null()

        mps_ptr => mps_head
        do iC=1,MPS(1)%nC_list
!          write(6,*)"MPS base ",mps_ptr%string
          mps_ptr => mps_ptr%p
        end do

        ! not unified MPS vector : additonal recording
        if(ieff.ne.0)then
!          write(6,*)"additonal recording" 
          do iroot=1,dmrg_nstates-1  
            do iC=1,MPS(iroot+1)%nC_list
              imatch=0
              !mps_ptr%p => null()
!              mps_tail => mps_head
              mps_ptr => mps_head
!              write(6,*)" == ",MPS(iroot+1)%position(iC) 
              do 
                if(.not.associated(mps_ptr)) exit
                if(trim(mps_ptr%string).eq.&
                   trim(MPS(iroot+1)%position(iC)))then
                  imatch=1     
                end if 
!                write(6,*)"MPS listing : ", mps_ptr%string 
                mps_ptr => mps_ptr%p
              end do  
              if(imatch.eq.0)then
!                write(6,*)"additonal recording -- invoke start" 
                !mps_ptr => mps_ptr%p
                !mps_head%p => null()
                !mps_ptr => mps_tail
                allocate(mps_tail%p)
                mps_tail => mps_tail%p
                mps_tail%p => null()
                mps_tail%string = trim(MPS(iroot+1)%position(iC))
!                write(6,*)" : ", trim(MPS(iroot+1)%position(iC)) 
!                write(6,*)"additonal recording -- invoke done" 
              end if        
            end do   
          end do   
        end if     

        ! Print to check the all effecitve MPS elements from all state
        !mps_ptr%p => null()
        open(unit=100,file="pickup.all")
        mps_ptr => mps_head
        nmps_list=0
        do
          if(.not.associated(mps_ptr)) exit
          nmps_list=nmps_list+1
          write(100,*)mps_ptr%string
          mps_ptr => mps_ptr%p
        end do
        close(100) 
        write(*,*)"Done the effecitve MPS elements part ", nmps_list
!        stop
!        deallocate(mps_ptr%p => null())
!        mps_head  => null()
!        mps_tail  => null()
!        mps_ptr%p => null()

        call system("python pickup_eff.py pickup.all")
        open(unit=100,file='pickup.pos')
        allocate(mpspos(nmps_list))
        do i=1,nmps_list
          mpspos(i)%pos(:)=0
          read(100,*)mpspos(i)%pos(1),mpspos(i)%pos(2),mpspos(i)%pos(3)
          mpspos(i)%ipos=&
              mpspos(i)%pos(1)*1000000+mpspos(i)%pos(2)*1000+&
              mpspos(i)%pos(3)
          !write(6,*)i,"-pos :", mpspos(i)%pos(1),mpspos(i)%pos(2),mpspos(i)%pos(3)  
        end do
        close(100)
        do i=1,nmps_list-1
          do j=i+1,nmps_list
            if(mpspos(i)%ipos.gt.mpspos(j)%ipos)then
              itmp=mpspos(j)%ipos
              mpspos(j)%ipos=mpspos(i)%ipos
              mpspos(i)%ipos=itmp

              ip1=mpspos(j)%pos(1)
              ip2=mpspos(j)%pos(2)
              ip3=mpspos(j)%pos(3)
              mpspos(j)%pos(1)=mpspos(i)%pos(1)
              mpspos(j)%pos(2)=mpspos(i)%pos(2)
              mpspos(j)%pos(3)=mpspos(i)%pos(3)
              mpspos(i)%pos(1)=ip1
              mpspos(i)%pos(2)=ip2
              mpspos(i)%pos(3)=ip3
           
            end if
          end do
        end do
!        do i=1,nlr_list
!          write(6,*)i,"-pos :", mpspos(i)%pos(1),mpspos(i)%pos(2),mpspos(i)%pos(3)  
!          call flush(6) 
!        end do
 
        do iroot=0,dmrg_nstates-1

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
          close(200)
          call system("chmod +x rdms-deri_stat")
          call system("./rdms-deri_stat")

          ! keep info. for all vec.
          allocate(mpsall(iroot+1)%ipos(nC_all))
          mpsall(iroot+1)%ipos=-100

          ! Important :
          ! This coefficient must be the same as that in abs(vjk)>1e-12 of ts_optimize.hpp
          ! If SA-DMRG, this part need no worries.
          ! Currently it may need to be considered,
          ! since the number of MPSci elements may be different for separated states
          if(.true.)then
            !write(*,*)"!!! > Readin coefficient from MPSCi.x.info"
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
                    ! in order to avoid the redundant
                    inum=inum+1
                    read(100,*)dv,ii,jj,kk
                    if(dabs(dv).gt.thresMPSHami)then
!                      write(6,*)thresMPSHami, dv, ii,jj,kk
                      iC=iC+1
                    end if
                    ! in order to avoid the redundant
                    mpsall(iroot+1)%ipos(inum)=&
                           ii*1000000+jj*1000+kk
                  end do
                end do
              end do
            close(100)

            MPS(iroot+1)%nC=iC
            nC=iC
            write(6,*)"nC, iC", nc, iC  
            ! Ci is a class for handling coupling
            ! May need to be extended
            if(.not.allocated(MPS(iroot+1)%ci))then
              allocate(MPS(iroot+1)%ci(nC)) 
            end if
            MPS(iroot+1)%ci%n=0 
            MPS(iroot+1)%ci%num=0 
            MPS(iroot+1)%ci%ipos=-1 
            MPS(iroot+1)%ci%list=0 

            iC=0
            open(unit=100,file=trim(string1))
              read(100,*)ntype
              inum=0
              do i=1,ntype
                read(100,*)j0,k0
                do j=1,j0
                  do k=1,k0
                    ! in order to avoid the redundant
                    inum=inum+1
                    read(100,*)dv,ii,jj,kk
                    if(dabs(dv).gt.thresMPSHami)then
                      ieff=0
                      ipos=-1
                      do l=1,nmps_list
                        if(mpspos(l)%pos(1).eq.ii)then
                          if(mpspos(l)%pos(2).eq.jj)then
                            if(mpspos(l)%pos(3).eq.kk)then
                              ieff=l
                              ipos=ii*1000000+jj*1000+kk
                            end if
                          end if
                        end if
                      end do
                      iC=iC+1
                      MPS(iroot+1)%Ci(iC)%dv   =   dv
                      MPS(iroot+1)%Ci(iC)%num  = inum
                      MPS(iroot+1)%Ci(iC)%list = ieff
                      MPS(iroot+1)%Ci(iC)%ipos = ipos
                      MPS(iroot+1)%Ci(iC)%n    = iC
                      write(6,*)iC,dv,ii,jj,kk,inum,ieff,ipos
                    end if
                    ! in order to avoid the redundant
                  end do
                end do
              end do
            close(100) 
          else
            string1=""
            string1="MPSCi_value."//trim(ctmp2)//".info"
            open(unit=100,file=trim(string1))
              inum=0
              do iC=1,nC_all
                read(100,*)dv,i,j,k
                if(dabs(dv).lt.thresMPSHami)then
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
          write(2,*)"     mapping index(num)  "
          write(2,*)MPS(iroot+1)%Ci%num
          write(2,*)"     mapping index(n)  "
          write(2,*)MPS(iroot+1)%Ci%n
          write(2,*)"     mapping index(list)  "
          write(2,*)MPS(iroot+1)%Ci%list
          write(2,*)"     mapping index(ipos)  "
          write(2,*)MPS(iroot+1)%Ci%ipos
          call flush(2)

        ! read in the sign info. from overlap
!           string1=""
!           string1="Loverlap."//trim(ctmp2)//".txt"
!           open(unit=100,file=trim(string1))
            do iC=1,nC
              !Leon: modified this to parse the MPS overlap output with a newline before the number
!               read(100,*)
!               read(100,*)
!               read(100,*)
!               read(100,*) dv
              MPS(iroot+1)%Ci(iC)%dsL = 1.0d0
            end do
!           close(100)
!           string1=""
!           string1="Roverlap."//trim(ctmp2)//".txt"
!           open(unit=100,file=trim(string1))
            do iC=1,nC
              !Leon: modified this to parse the MPS overlap output with a newline before the number
!               read(100,*)
!               read(100,*)
!               read(100,*)
!               read(100,*) dv
              MPS(iroot+1)%Ci(iC)%dsR = 1.0d0
              !write(2,*)MPS(iroot+1)%Ci(iC)%dsR
            end do
!           close(100)

!          write(2,*)"Sign correct vector L"
!          call print_mat(1,nc,MPS(iroot+1)%Ci%dsL,2)
!          call flush(2)
!          write(2,*)"Sign correct vector R"
!          call print_mat(1,nc,MPS(iroot+1)%Ci%dsR,2)
!          call flush(2)
!          write(2,*)"Attention!! Unified the Sign L and R using L"
!          write(2,*)"Attention!! Unified the Sign L and R using L"
!          MPS(iroot+1)%Ci%dsR = MPS(iroot+1)%Ci%dsL


          call CS_dim_fro2cls()

! For read in 1-RDM-deri (<i|o|mps> part)
          iC=0
          open(unit=100,file="info1.tmp")
            do iC0=1,MPS(iroot+1)%nC_list
              oneRDMr=""

              read(100,*,iostat=ierr)oneRDMr
              if(ierr.ne.0)exit

              ieff=0         
              ipos=MPS(iroot+1)%pos(iC0,4) 
              do ii=1,MPS(iroot+1)%nC
!                write(6,*)"ipos",ipos,"ii-pos",MPS(iroot+1)%Ci(ii)%ipos
                if(ipos.eq.MPS(iroot+1)%Ci(ii)%ipos)then
                  iC=iC+1
                  ieff=iC
                end if  
              end do 
         
              write(6,*)"iC0, iC : ",iC0, iC 
              
              if(ieff.ne.0)then
! symmetrized 1-RDM-deri
                allocate(MPS(iroot+1)%Ci(iC)%D(nact,nact))
                MPS(iroot+1)%Ci(iC)%D=0.0d0
!              write(6,*)MPS(iroot+1)%Ci(iC)%D
!              call flush(6)
!                write(6,*)"iC in read, nact",iC,nact                
! Left side part
                allocate(MPS(iroot+1)%Ci(iC)%DL(nact,nact))
                MPS(iroot+1)%Ci(iC)%DL=0.0d0
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
!                write(2,*)"1-RDM-deri L for element",iC
!                call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%DL,2)
!                call flush(2)
              end if

            end do
          close(100)
          write(2,*)"Finish importing 1-RDM-deri L for state",iroot
          call flush(2)

! For read in 1-RDM-deri (<mps|o|i> part)
          iC=0
          open(unit=100,file="info3.tmp")
            do iC0=1,MPS(iroot+1)%nC_all
              oneRDMr=""

              read(100,*,iostat=ierr)oneRDMr
              if(ierr.ne.0)exit

              ieff=0
              ipos=MPS(iroot+1)%pos(iC0,4)
              do ii=1,MPS(iroot+1)%nC
!                write(6,*)"ipos",ipos,"ii-pos",MPS(iroot+1)%Ci(ii)%ipos
                if(ipos.eq.MPS(iroot+1)%Ci(ii)%ipos)then
                  iC=iC+1
                  ieff=iC
                end if
              end do

              if(ieff.ne.0)then
! Right side part
                allocate(MPS(iroot+1)%Ci(iC)%DR(nact,nact))
                MPS(iroot+1)%Ci(iC)%DR=0.0d0
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
!                write(2,*)"1-RDM-deri R for element",iC
!                call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%DR,2)
!                call flush(2)
              end if
            end do
          close(100)
          write(2,*)"Finish importing 1-RDM-deri R for state",iroot
          call flush(2)


! For read in 2-RDM-deri (L)
          iC=0  
          open(unit=100,file="info2.tmp")
            do iC0=1,MPS(iroot+1)%nC_all
              twoRDMr=""

              read(100,*,iostat=ierr)twoRDMr
              if(ierr.ne.0)exit

              ieff=0
              ipos=MPS(iroot+1)%pos(iC0,4)
              do ii=1,MPS(iroot+1)%nC
!                write(6,*)"ipos",ipos,"ii-pos",MPS(iroot+1)%Ci(ii)%ipos
                if(ipos.eq.MPS(iroot+1)%Ci(ii)%ipos)then
                  iC=iC+1
                  ieff=iC
                end if
              end do              

              if(ieff.ne.0)then
! Left side part
                allocate(MPS(iroot+1)%Ci(iC)%PL(nact,nact,nact,nact))
                MPS(iroot+1)%Ci(iC)%PL=0.0d0

!                write(322,*)"2-RDM-deri (L)(all 0) for iC",iC
!                call print_gat &
!                (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PL,322)

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
!                  MPS(iroot+1)%Ci(iC)%PL(&
!                            iorder(kl+1),&
!                            iorder(li+1),&
!                            iorder(ij+1),&
!                            iorder(jk+1))=dv
!                  MPS(iroot+1)%Ci(iC)%PL(&
!                            iorder(jk+1),&
!                            iorder(ij+1),&
!                            iorder(li+1),&
!                            iorder(kl+1))=dv
!                  MPS(iroot+1)%Ci(iC)%PL(&
!                            iorder(li+1),&
!                            iorder(kl+1),&
!                            iorder(jk+1),&
!                            iorder(ij+1))=dv
                    write(11210,*)ic,"iC",ij,jk,kl,li,gv
                  end do
                close(101)

!                write(322,*)"2-RDM-deri (L) for iC",iC
!              call print_gat &
!              (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PL,322)

                MPS(iroot+1)%Ci(iC)%PL=MPS(iroot+1)%Ci(iC)%PL*2.0d0
              end if
            end do
          close(100)
          write(2,*)"Finish importing 2-RDM-deri L for state",iroot
          call flush(2)

! For read in 2-RDM-deri (R)
          iC=0
          open(unit=100,file="info4.tmp")
            do iC0=1,MPS(iroot+1)%nC_all
              twoRDMr=""

              read(100,*,iostat=ierr)twoRDMr
              if(ierr.ne.0)exit

              ieff=0
              ipos=MPS(iroot+1)%pos(iC0,4)
              do ii=1,MPS(iroot+1)%nC
!                write(6,*)"ipos",ipos,"ii-pos",MPS(iroot+1)%Ci(ii)%ipos
                if(ipos.eq.MPS(iroot+1)%Ci(ii)%ipos)then
                  iC=iC+1
                  ieff=iC
                end if
              end do              

              if(ieff.ne.0)then
                allocate(MPS(iroot+1)%Ci(iC)%PR(nact,nact,nact,nact))
                MPS(iroot+1)%Ci(iC)%PR=0.0d0
!                write(324,*)"2-RDM-deri (R) (all 0) for iC",iC
!                call print_gat &
!                (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PR,324)
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
!                  MPS(iroot+1)%Ci(iC)%PR(&
!                            iorder(kl+1),&
!                            iorder(li+1),&
!                            iorder(ij+1),&
!                            iorder(jk+1))=dv
!                  MPS(iroot+1)%Ci(iC)%PR(&
!                            iorder(jk+1),&
!                            iorder(ij+1),&
!                            iorder(li+1),&
!                            iorder(kl+1))=dv
!                  MPS(iroot+1)%Ci(iC)%PR(&
!                            iorder(li+1),&
!                            iorder(kl+1),&
!                            iorder(jk+1),&
!                            iorder(ij+1))=dv
!                    write(11220,*)ic,"iC",ij,jk,kl,li,dv
                  end do
                close(101)

!                write(324,*)"2-RDM-deri (R) for iC",iC
!                call print_gat &
!                (nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%PR,324)

                MPS(iroot+1)%Ci(iC)%PR=MPS(iroot+1)%Ci(iC)%PR*2.0d0
              end if
            end do
          close(100)
          write(2,*)"Finish importing 2-RDM-deri R for state",iroot
          call flush(2)

! Closed shell update for all RDM-deri

          do iC=1,MPS(iroot+1)%nC

            call CS_dim_fro2cls()
            write(*,*)"Start the ", iC,"iterations for RDM-deri"

            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
            write(6,*)"iroot, nact, iC", iroot, nact, iC
!            do i=1,nact
!              do j=1,nact
!                write(6,*)MPS(iroot+1)%Ci(iC)%dL(i,j)
!                call flush(6) 
!              end do
!            end do
!            call flush(6)
!            call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%dL,6)
!            call flush(6)
!            call print_gat(nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%p,6)

            ! L
! Fix the sign for RDM-deri (closed shell part)

            dv = dabs(MPS(iroot+1)%Ci(iC)%dv) ! Leon to Yingjin: why is the closed shell update without sign?
            ! Leon: Maybe it makes sense to merge the closed shell updates or am I missing something?
            call closed_shell_update&
            (MPS(iroot+1)%Ci(iC)%dL,MPS(iroot+1)%Ci(iC)%PL,TM1,GM1,dv)  ! CCCC
            deallocate(MPS(iroot+1)%Ci(iC)%dL,MPS(iroot+1)%Ci(iC)%PL)
            allocate(MPS(iroot+1)%Ci(iC)%dL(nocc,nocc))
            allocate(MPS(iroot+1)%Ci(iC)%PL(nocc,nocc,nocc,nocc))
            MPS(iroot+1)%Ci(iC)%dL=TM1
            MPS(iroot+1)%Ci(iC)%pL=GM1
            deallocate(TM1,GM1)
!            write(6,*)iroot,"MPS(iroot+1)%Ci(iC)%dL"
!            call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%dL,6)
!            call flush(6)
!           call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%pL,6)
!            write(2,*)iC,"-th sign",MPS(iroot+1)%Ci(iC)%sign

          end do

          do iC=1,MPS(iroot+1)%nC

            call CS_dim_fro2cls()
!            write(*,*)"Start the ", iC,"iterations for RDM-deri"

            ! Closed shell update
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
!            call print_mat(nact,nact,MPS(iroot+1)%Ci(iC)%dR,6)
!            call print_gat(nact,nact,nact,nact,MPS(iroot+1)%Ci(iC)%p,6)

            ! R
! Fix the sign for RDM-deri (closed shell part)
            dv = dabs(MPS(iroot+1)%Ci(iC)%dv) ! see above

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
          do iC=1,MPS(iroot+1)%nC
! Avergaed
            allocate(MPS(iroot+1)%Ci(iC)%P(nact,nact,nact,nact))
            MPS(iroot+1)%Ci(iC)%P=0.0d0

            write(2,*)"Davidson C_",iC,"(",MPS(iroot+1)%Ci(iC)%n,") : ", &
                       MPS(iroot+1)%Ci(iC)%dv, &
                       MPS(iroot+1)%Ci(iC)%dsL, MPS(iroot+1)%Ci(iC)%dsR
!                 dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsL, &
!                 dabs(MPS(iroot+1)%Ci(iC)%dv)*MPS(iroot+1)%Ci(iC)%dsR
!            call flush(2)

              ! Fix the sign for RDM-deri
              if(MPS(iroot+1)%Ci(iC)%dv.lt.0.0d0)then
                MPS(iroot+1)%Ci(iC)%DL=-MPS(iroot+1)%Ci(iC)%DL
                MPS(iroot+1)%Ci(iC)%PL=-MPS(iroot+1)%Ci(iC)%PL
              end if
              if(MPS(iroot+1)%Ci(iC)%dv.lt.0.0d0)then
                MPS(iroot+1)%Ci(iC)%DR=-MPS(iroot+1)%Ci(iC)%DR
                MPS(iroot+1)%Ci(iC)%PR=-MPS(iroot+1)%Ci(iC)%PR
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

            write(1121,*)iroot+1,"state",iC,"-th tRDM derivative (L) "
            call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%PL,1121)
            write(1122,*)iroot+1,"state",iC,"-th tRDM derivative (R) "
            call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%PR,1122)

          end do

          write(2,*)"resembled oRDM "
          call print_mat(nocc,nocc,oRDMchk,2)

!          write(12,*)"resembled tRDM "
!          call print_gat(nocc,nocc,nocc,nocc,tRDMchk,12)

          write(2,*)" !! Corrected RDM-deri (L & R) is used. !!"
          call flush(2)

! to Leon :
!       read-in normal RDM for each state

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
            write(15,*)"the closed-shell_updated 1-RDMs",iroot
            call print_mat1(nocc,nocc,MPS(iroot+1)%d,15)
            write(16,*)"the closed-shell_updated 2-RDMs",iroot
            call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%p,16)
            deallocate(TM1,GM1)

!            call CS_dim_fro2cls()

! to Leon :
!         Check the energy for each state
          write(6,*)"orb%nsub,orb%act,orb%total,nact,norb"
          write(6,*)orb%nsub,orb%act,orb%total,nact,norb  
          call energy_gen(orb%nsub,orb%act,orb%total,nact,norb, &
                             T,U,MPS(iroot+1)%D,MPS(iroot+1)%P, &
                                 MPS(iroot+1)%RDM_energy)

          write(2,*)"DMRG(RDM) energy for state",iroot
          write(2,*) MPS(iroot+1)%RDM_energy
          write(2,*) MPS(iroot+1)%RDM_energy+FRONRE0

! to Leon :
!          All MPSci gradients should exactly be 0 since the converged SA-DMRG-SCF

! Check the MPSci gradients for each Davidson vectors

          do iC=1,MPS(iroot+1)%nC

!            call CS_dim_fro2cls()
!            write(*,*)"Start the ", iC,"iterations"

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

             write(111,*)iroot+1,"state",iC,"-th oRDM derivative (sym) "
             call print_mat(nocc,nocc,MPS(iroot+1)%Ci(iC)%D,111)

!             write(112,*)iroot+1,"state",iC,"-th tRDM derivative (sym) "
!             call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%Ci(iC)%P,112)

             write(111,*)" 1e integrals T"
             call print_mat(norb,norb,T,111)

!             write(112,*)" 2e integrals U"
!             call print_gat(norb,norb,norb,norb,U,112)

          ! Gradient that calculated by energy_gen
            MPS(iroot+1)%Ci(iC)%grad=0.0d0
            call energy_gen(orb%nsub,orb%act,orb%total,nact,norb,T,U,  &
                         MPS(iroot+1)%Ci(iC)%D,MPS(iroot+1)%Ci(iC)%P,  &
                         MPS(iroot+1)%Ci(iC)%grad)

            write(2,*)&
          "Energy component(mpsci, ref) & gradient for Davidson vectors",iC
            write(2,*) MPS(iroot+1)%Ci(iC)%grad, &
                       MPS(iroot+1)%RDM_energy*MPS(iroot+1)%Ci(iC)%dv, &
                       MPS(iroot+1)%Ci(iC)%grad &
                      -MPS(iroot+1)%RDM_energy*MPS(iroot+1)%Ci(iC)%dv

            dv=MPS(iroot+1)%Ci(iC)%grad &
              -MPS(iroot+1)%RDM_energy*MPS(iroot+1)%Ci(iC)%dv

            if(dabs(dv).gt.1.0e-8)then
              !write(6,*)orb%nsub,nact,norb
              !write(6,*)orb%act
              !write(6,*)orb%total
         write(6,*)"==================================================="
         write(6,*)"Derivation found in the RDM-derivatives calculation"
         write(6,*)"State",iroot+1,"the",iC,"-th MPSci element"
         write(6,*)"     ---- !!!!! warning warning !!!!! ---- "
         write(6,*)"The mismatch is",dv
         write(6,*)"==================================================="
              call flush(6)
!              stop
            end if

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
          call print_mat1(nocc,nocc,MPS(iroot+1)%D,2)
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

!          write(6,*)"The reformed two-RDMs for checking, fort.18"
!          call print_mat(nocc,nocc,MPS(iroot+1)%d,6)
!          call print_gat(nocc,nocc,nocc,nocc,MPS(iroot+1)%p,18)
!          stop
!          call CS_dim_fro2cls()

        end do

!        stop "All info read in for reduced mpsci" 

        ! unified LR MPS vector : record in ascending order list
        do iroot=0,dmrg_nstates-1
          do iC=1,MPS(iroot+1)%nC

            allocate(lr_ptr)
            lr_ptr%ipos=MPS(iroot+1)%Ci(iC)%ipos
!            write(6,*)iC," - ", lr_ptr%ipos

            if(.not.associated(lr_head))then
              lr_head  => lr_ptr
              lr_tail  => lr_head
              nullify(lr_ptr%p)
            else  
              if(lr_ptr%ipos.lt.lr_head%ipos)then
                lr_ptr%p  => lr_head
                lr_head   => lr_ptr
              else if(lr_ptr%ipos.gt.lr_tail%ipos)then
                lr_tail%p => lr_ptr
                lr_tail   => lr_ptr
                nullify(lr_tail%p)
              else if(lr_ptr%ipos.eq.lr_tail%ipos)then
                lr_tail   => lr_ptr
              else  
                lr_ptr1   => lr_head 
                lr_ptr2   => lr_ptr1%p
                do 
                  if((lr_ptr%ipos.ge.lr_ptr1%ipos).and.&
                     (lr_ptr%ipos.lt.lr_ptr2%ipos))then
                    if(lr_ptr%ipos.eq.lr_ptr1%ipos)then
                      exit
                    else        
                      lr_ptr%p  => lr_ptr2 
                      lr_ptr1%p => lr_ptr
                      exit
                    end if  
                  end if 
                  lr_ptr1 => lr_ptr2       
                  lr_ptr2 => lr_ptr2%p
                end do 
              end if         
            end if

          end do
        end do  

        ! Print to check the all effecitve LR mps elements from all state
        !lr_ptr%p => null()
        lr_ptr => lr_head
        nlr_list=0
        do
          if(.not.associated(lr_ptr)) exit
          nlr_list=nlr_list+1
          write(6,*)nlr_list, " lr mps vec. ", lr_ptr%ipos
          lr_ptr => lr_ptr%p
        end do
        L_R%nC=nlr_list 
!        stop "Done the effecitve LR mps elements part "

        do iroot=0,dmrg_nstates-1 
          do iC=1,MPS(iroot+1)%nC_list
            ipos=MPS(iroot+1)%pos(iC,4)
            imatch=0

            lr_ptr => lr_head
            ii=0
            do
              if(.not.associated(lr_ptr)) exit
              ii=ii+1
              if(lr_ptr%ipos.eq.ipos)then
                imatch=ii
              end if  
              lr_ptr => lr_ptr%p
            end do
            if(imatch.ne.0)then
              MPS(iroot+1)%pos(iC,5)=iC
            end if
          end do  
!          write(6,*)"pos 5 :",MPS(iroot+1)%pos(:,5)
        end do

        write(6,*)"MPS all",mpsall(1)%ipos

        allocate(L_R%idx(L_R%nC))
        allocate(L_R%pos(L_R%nC))
        L_R%idx=0
        L_R%pos=-100
        lr_ptr => lr_head
        ii=0
        do
          if(.not.associated(lr_ptr)) exit
          ii=ii+1
          write(6,*)ii,lr_ptr%ipos
          imatch=0
          do iC=1,nC_all
            ipos=mpsall(iroor+1)%ipos(iC)
            if(lr_ptr%ipos.eq.ipos)then
              imatch=iC
              exit
            end if
          end do
          if(imatch.ne.0)then
            L_R%idx(ii)=iC
            L_R%pos(ii)=ipos
          else
            stop "MPS vec missing or wrong"
          end if
          lr_ptr => lr_ptr%p
        end do
        write(6,*)"L_R%idx : ",L_R%idx(:)
        write(6,*)"L_R%idx : ",L_R%pos(:)
!        stop 

        do iroot=0,dmrg_nstates-1

          write(6,*)"I'm not sure if this is needed, open eyes" 
          call CS_dim_cls2fro()          

! Allocate the orb-CI and CI-CI Hessian
          allocate(MPS(iroot+1)%HCR(nlr_list,norb,norb)); MPS(iroot+1)%HCR=0.0d0
          allocate(MPS(iroot+1)%KCR(nlr_list,norb,norb)); MPS(iroot+1)%KCR=0.0d0
          allocate(MPS(iroot+1)%HCC(nlr_list,nlr_list)) ; MPS(iroot+1)%HCC=0.0d0
! Allocate for the full Hamiltonian
          allocate(MPS(iroot+1)%Hami(nlr_list,nlr_list)) ; MPS(iroot+1)%Hami=0.0d0
          allocate(MPS(iroot+1)%Hvec(nlr_list,nlr_list)) ; MPS(iroot+1)%Hvec=0.0d0
          allocate(MPS(iroot+1)%Hval(nlr_list))          ; MPS(iroot+1)%Hval=0.0d0

! to Leon :
!         read in the full local Hamiltonian

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
          do iC=1,nlr_list
            do jC=1,nlr_list
              MPS(iroot+1)%Hami(iC,jC)= &
               Hami_all(L_R%idx(iC),L_R%idx(jC))
            end do
          end do

          write(13,*)"Hamiltonian Matrix all "
          call print_mat(nC_all,nC_all,Hami_all,13)
          call flush(13)

          write(13,*)"Hamiltonian Matrix "
          call print_mat(nlr_list,nlr_list,MPS(iroot+1)%Hami,13)
          call flush(13)
          deallocate(Hami_all)

! Set very small value to 0
          do i=1,nlr_list
            do j=1,nlr_list
              if(dabs(MPS(iroot+1)%Hami(i,j)).lt.thresMPSHami)then
                MPS(iroot+1)%Hami(i,j)=0.0d0
              end if
            end do
          end do
          write(2,*)"Hamiltonian Matrix "
          call print_mat(nlr_list,nlr_list,MPS(iroot+1)%Hami,2)
          call flush(2)

! Re-calculate the RDM_energies base on the refactored RDMs

          call energy_gen(orb%nsub,orb%act,orb%total,nact,norb,T,U,  &
                                     MPS(iroot+1)%D,MPS(iroot+1)%P,  &
                                     MPS(iroot+1)%RDM_energy)

          write(2,*)"The REM-refactored E (no-NRE and with-NRE) is", &
                    MPS(iroot+1)%RDM_energy, MPS(iroot+1)%RDM_energy+FRONRE0

! to Leon :
!          Construct the MPS-MPS Hessian in LR equations

            do iC=1,nlr_list
              do jC=1,nlr_list
                MPS(iroot+1)%HCC(iC,jC)=MPS(iroot+1)%Hami(iC,jC)
              end do

!              write(6,*)" FRONRE,FRONRE0,FRONREA",FRONRE,FRONRE0,FRONREA
!
!              write(6,*)" = = = = = = = = = = "
!              write(6,*)" iroot, iC, jC ",iroot,iC,jC
!              write(6,*)MPS(iroot+1)%HCC(iC,iC),MPS(iroot+1)%RDM_energy

              MPS(iroot+1)%HCC(iC,iC)=MPS(iroot+1)%HCC(iC,iC)&
                                     +FRONREA                &
                                     -MPS(iroot+1)%RDM_energy&
                                     -FRONRE0

!              write(6,*)"  updated HCC (no state-transfer part)"

!              write(6,*) MPS(iroot+1)%HCC(iC,iC)
!              write(6,*)" = = = = = = = = = = "

            end do

          MPS(iroot+1)%HCC=MPS(iroot+1)%HCC*dmrg_weight(iroot+1)

       ! Construct the A(Fock) matrix for each state
        call fock_gen(orb%nsub,orb%act,orb%total,     &
                      nact,norb,T,U,MPS(iroot+1)%D,   &
                      MPS(iroot+1)%P,MPS(iroot+1)%A)

          write(2,*)"The A-matrix for state",iroot
          call print_mat(norb,norb,MPS(iroot+1)%A,2)
          write(21,*)"The A-matrix for state",iroot
          call print_mat(norb,norb,MPS(iroot+1)%A,21)

! to Leon :
!         orbital gradients for each specific state

          do i=1,norb
            do j=1,norb
              MPS(iroot+1)%G(i,j) = &
              MPS(iroot+1)%A(i,j)-MPS(iroot+1)%A(j,i)
            end do
          end do

          ! Construct the A^I matrix for every Davidson elemments for each state
          allocate(TM1(norb,norb));TM1=0.0d0
          allocate(TM2(norb,norb));TM2=0.0d0
          do iC=1,MPS(iroot+1)%nC
            allocate(MPS(iroot+1)%ci(iC)%A(norb,norb))
            allocate(MPS(iroot+1)%ci(iC)%AL(norb,norb))
            allocate(MPS(iroot+1)%ci(iC)%AR(norb,norb))
            MPS(iroot+1)%ci(iC)%A=0.0d0
            MPS(iroot+1)%ci(iC)%AL=0.0d0
            MPS(iroot+1)%ci(iC)%AR=0.0d0

            call fock_gen(orb%nsub,orb%act,orb%total,     &
                   nact,norb,T,U,MPS(iroot+1)%Ci(iC)%D,   &
                   MPS(iroot+1)%Ci(iC)%P,MPS(iroot+1)%Ci(iC)%A)

            write(21,*)"The A^I-matrix for Davidson elem",iC
            call print_mat(norb,norb,MPS(iroot+1)%Ci(iC)%A,21)
            TM1=TM1+MPS(iroot+1)%Ci(iC)%A
            TM2=TM2+MPS(iroot+1)%Ci(iC)%A*MPS(iroot+1)%Ci(iC)%dv

          end do  ! loop the MPS_ci elements
          write(21,*)"The sum A-matrix for all Davidson vectors"
          call print_mat(norb,norb,TM1,21)
          write(21,*)"The sum A-matrix for all Davidson vectors with c"
          call print_mat(norb,norb,TM2,21)
          deallocate(TM1)
          deallocate(TM2)

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

! to Leon : closed shell approximation related subroutines
          ! Recover the old closed shell structure
          call CS_dim_fro2cls()

        end do  ! Finish read-in all properties for every states

! Some deallocate
        deallocate(mpspos)
        
! ===============================================
!    check the closed-shell updated RDMs : start
! ===============================================

        allocate(TM1(nocc,nocc))
        allocate(GM1(nocc,nocc,nocc,nocc))
        TM1=0.0d0
        GM1=0.0d0
        do iroot=0,dmrg_nstates-1
          ir=iroot+1
          write(6,*)" the ir ",ir
          call flush(6)
          TM1=TM1+MPS(ir)%d*dmrg_weight(ir)
          GM1=GM1+MPS(ir)%P*dmrg_weight(ir)
        end do

        write(15,*)"the closed-shell_updated averaged 1-RDMs"
        call print_mat1(nocc,nocc,TM1,15)
        write(16,*)"the closed-shell_updated averaged 2-RDMs"
        call print_gat(nocc,nocc,nocc,nocc,GM1,16)
        deallocate(TM1,GM1)

! ===============================================
!    check the closed-shell updated RDMs : done
! ===============================================

        write(6,*)" Finish read-in all properties for every states"
        call flush(6)

! ===============================================
!        The SA state-transition portion
! ===============================================
        SA_RDMenergy=0.0d0
        do iroot=0,dmrg_nstates-1
          SA_RDMenergy=SA_RDMenergy&
            +MPS(iroot+1)%RDM_energy*dmrg_weight(iroot+1)
        end do

!        MPS(iroot+1)%Ci(iC)%num

        do iroot=0,dmrg_nstates-1

!          write(6,*)"==== Hcc (1st term in Eq.30) ====",iroot,"state"
!          call flush(6)
          irt=iroot+1
          nCi=MPS(iroot+1)%nC
          write(2,*)"==== Hcc (1st term in Eq.30) ====",iroot,"state"
          call print_mat1(nlr_list,nlr_list,MPS(iroot+1)%HCC,2)
          call flush(2)

          allocate(TM1(nlr_list,nlr_list))
          TM1=0.0d0

          do jroot=0,dmrg_nstates-1
            jrt=jroot+1
            nCj=MPS(jrt)%nC

            !if(nCi.ne.nCj)then
            !  write(6,*)"The local MPS basis for state-", iroot,&
            !  " and state-",jroot,"is different"
            !  stop 
            !end if         

            call trace(nlr_list,MPS(irt)%Hami,dv1)
            call trace(nlr_list,MPS(jrt)%Hami,dv2)

            if(dabs(dv1-dv2).lt.1.0e-8)then
              ! If the same local Hamitlonian, e.g. same irreps
              write(6,*)"State",iroot," and ",jroot,&
              "should in a same irreps"
              do iC=1,nCj

                ! check if this MPS element exist
!                iexist=0
!                do iCx=1,nCi
!                  if(MPS(jrt)%Ci(iC)%num.eq.MPS(irt)%Ci(iCx)%num)then
!                    iexist=1
!                  end if
!                end do
!                if(iexist.eq.1)then

                  do jC=1,nCj

                    ! check if this MPS element exist
!                    jexist=0
!                    do jCx=1,nCi
!                      if(MPS(jrt)%Ci(jC)%num.eq.MPS(irt)%Ci(jCx)%num)then
!                        jexist=1
!                      end if
!                    end do
!                    if(jexist.eq.1)then

!                    if(MPS(jrt)%Ci(jC)%num)
                      dv =0.0d0
                      dv1=0.0d0
                      dv2=0.0d0
                      dv1=dv1+MPS(jrt)%Ci(iC)%dv!*dmrg_weight(jrt)
                      dv2=dv2+MPS(jrt)%Ci(jC)%dv!*dmrg_weight(jrt)
                      dv=MPS(iroot+1)%RDM_energy-MPS(jroot+1)%RDM_energy

                     ! This is assume MPSx have the same Local Hamiltonian, careful !! need further check
                     ! This is assume MPSx have the same Local Hamiltonian, careful !! need further check
                     ! This is assume MPSx have the same Local Hamiltonian, careful !! need further check
                      iipos=MPS(jrt)%Ci(iC)%ipos
                      jjpos=MPS(jrt)%Ci(jC)%ipos

                      call &
                      get_position(iipos,L_R%nC,L_R%idx,L_R%pos,iiC,0)
                      call &
                      get_position(jjpos,L_R%nC,L_R%idx,L_R%pos,jjC,0)
!                      write(6,*)iipos,jjpos,"iiC,jjC",iiC,jjC
!                      call flush(6)
!                      stop

                      TM1(iiC,jjC)=TM1(iiC,jjC)+dv1*dv2*dv*dmrg_weight(jrt)
!                    end if

                  end do
!                end if

              end do  ! for the iC-th jroot
            else
              ! If the different local Hamitlonian, e.g. different irreps
              write(6,*)"State",iroot," and ",jroot,&
              "should in two different irreps"
              write(2,*)"State",iroot," and ",jroot,&
              "should in two different irreps"
            end if
            write(2,*)"==== TM1 (2nd term in Eq.30) ====",iroot,jroot
            call print_mat1(nlr_list,nlr_list,TM1,2)
            write(2,*)"=================================",iroot,jroot
            call flush(2)
          end do  ! for iroot

          call print_mat1(nlr_list,nlr_list,TM1,2)
          call flush(2)
          call print_mat1(nlr_list,nlr_list,MPS(iroot+1)%HCC,2)
          call flush(2)

!              write(2,*)iC,jC," dv1, dv2 : ",dv1,dv2
!              write(2,*)"E(i),E^sa",MPS(iroot+1)%RDM_energy,SA_RDMenergy
          MPS(iroot+1)%HCC=MPS(iroot+1)%HCC+TM1
          deallocate(TM1)
          write(2,*)"finish loop and deallocated TM1"
          call flush(2)
        end do
!        stop "HCC almost done"

! Final step for Hcc part; make it separately
        do iroot=0,dmrg_nstates-1
          MPS(iroot+1)%HCC=MPS(iroot+1)%HCC*2.0d0
          !nc=MPS(iroot+1)%nC
          write(2,*)"==== the Hcc part ====",iroot,"state"
          call print_mat1(nlr_list,nlr_list,MPS(iroot+1)%HCC,2)
          call flush(2)
        end do

        ! Calculate the preconditioner for MPS part (not used here)
        call SA_prec_mps()

! ================================================
! Read-in state-transfer properties
! ================================================
!       re-coding this part later
! ================================================
! Sign problem after SA-solver?  2018.9.21     
!   (Force the sign of TDMs foe testing)        
! ================================================

!        signTDM=-1.0d0

        allocate(MPSTDMs(dmrg_nstates,dmrg_nstates))

! use a simple CI routine to check the TDMs

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
            dv1=0.0d0  ! trace of TDMs
            open(unit=101,file=trim(string1))
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  MPSTDMs(ir,jr)%D(i,j)=dv
                end do
                dv1=dv1+dabs(MPSTDMs(ir,jr)%D(i,i)) 
              end do
              ! Hotfix sign if using sa-solver  
              dv1=sign(dv1,MPSTDMs(ir,jr)%D(1,1)) ! check the sign 
            close(101)
            write(2,*)"oneTDM ",iroot,jroot
!            if(dv1.lt.0.0d0)then
!              MPSTDMs(ir,jr)%D=-1.0d0*MPSTDMs(ir,jr)%D
!              write(2,*)"==> The TDMs(D) sign is changed <==" 
!            end if   
            MPSTDMs(ir,jr)%D=signTDM*MPSTDMs(ir,jr)%D
            call print_mat1(nact,nact,MPSTDMs(ir,jr)%D,2)
            call flush(2)

            ! oneRDM.jroot.iroot
            string1=""
            string1="oneRDM."//trim(ctmp3)//"."//trim(ctmp2)
            dv2=0.0d0  ! trace of TDMs
            open(unit=101,file=trim(string1))
              read(101,*)ij0
              do i=1,nact
                do j=1,nact
                  read(101,*)ij,ji,dv
                  MPSTDMs(jr,ir)%D(i,j)=dv
                end do
                dv2=dv2+dabs(MPSTDMs(jr,ir)%D(i,i))
              end do
              dv2=sign(dv2,MPSTDMs(jr,ir)%D(1,1)) ! check the sign 
            close(101)
            write(2,*)"oneTDM ",jroot,iroot
            if(dv1+dv2.lt.1.0d-6)then
              MPSTDMs(jr,ir)%D=-1.0d0*MPSTDMs(jr,ir)%D
              write(2,*)"==> The TDMs(D) sign is changed <==" 
            end if       
            MPSTDMs(jr,ir)%D=signTDM*MPSTDMs(jr,ir)%D       
            call print_mat1(nact,nact,MPSTDMs(jr,ir)%D,2)
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
!                MPSTDMs(ir,jr)%P(      &
!                          iorder(kl+1),&
!                          iorder(li+1),&
!                          iorder(ij+1),&
!                          iorder(jk+1))=dv
!                MPSTDMs(ir,jr)%P(      &
!                          iorder(jk+1),&
!                          iorder(ij+1),&
!                          iorder(li+1),&
!                          iorder(kl+1))=dv
!                MPSTDMs(ir,jr)%P(      &
!                          iorder(li+1),&
!                          iorder(kl+1),&
!                          iorder(jk+1),&
!                          iorder(ij+1))=dv
              end do
            close(101)
            MPSTDMs(ir,jr)%P=MPSTDMs(ir,jr)%P*2.0d0
            MPSTDMs(ir,jr)%P=signTDM*MPSTDMs(ir,jr)%P
!            if(dv1.lt.0.0d0)then
!              MPSTDMs(ir,jr)%P=-1.0d0*MPSTDMs(ir,jr)%P
!              write(2,*)"==> The TDMs(P) sign is changed <==" 
!            end if         

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
!                MPSTDMs(jr,ir)%P(      &
!                          iorder(kl+1),&
!                          iorder(li+1),&
!                          iorder(ij+1),&
!                          iorder(jk+1))=dv
!                MPSTDMs(jr,ir)%P(      &
!                          iorder(jk+1),&
!                          iorder(ij+1),&
!                          iorder(li+1),&
!                          iorder(kl+1))=dv
!                MPSTDMs(jr,ir)%P(      &
!                          iorder(li+1),&
!                          iorder(kl+1),&
!                          iorder(jk+1),&
!                          iorder(ij+1))=dv
              end do
            close(101)
            MPSTDMs(jr,ir)%P=MPSTDMs(jr,ir)%P*2.0d0
            !if(dv2.lt.0.0d0)then
            if(dv1+dv2.lt.1.0d-6)then
              MPSTDMs(jr,ir)%P=-1.0d0*MPSTDMs(jr,ir)%P
              write(2,*)"==> The TDMs(P) sign is changed <==" 
            end if         
            MPSTDMs(jr,ir)%P=signTDM*MPSTDMs(jr,ir)%P

            write(2,*)"Finished all read-in TDMs "
            call flush(2)

            ! Symmetrilize TDMs (Need)  ! option-1
            !MPSTDMs(ir,jr)%D=(MPSTDMs(ir,jr)%D+MPSTDMs(jr,ir)%D)/2.0d0
            !MPSTDMs(jr,ir)%D= MPSTDMs(ir,jr)%D

            !MPSTDMs(ir,jr)%P=(MPSTDMs(ir,jr)%P+MPSTDMs(jr,ir)%P)/2.0d0
            !MPSTDMs(jr,ir)%P= MPSTDMs(ir,jr)%P

            ! ====  Closed shell update ====
            allocate(TM1(nocc,nocc))
            allocate(GM1(nocc,nocc,nocc,nocc))
            allocate(TM2(nocc,nocc))
            allocate(GM2(nocc,nocc,nocc,nocc))
            TM1=0.0d0
            TM2=0.0d0
            GM1=0.0d0
            GM2=0.0d0
            !write(6,*)"The symmetrical TDMs (1-body)"
            !call print_mat(nact,nact,MPSTDMs(ir,jr)%d,6)
            !write(6,*)"The symmetrical TDMs (2-body)"
            !call print_gat(nact,nact,nact,nact,MPSTDMs(ir,jr)%p,6)

            ! the overlap between state-i and state-j
            ! Assuming it is 0
            dvphiphj=0.0d0

            call closed_shell_update&
            (MPSTDMs(ir,jr)%d,MPSTDMs(ir,jr)%P,TM1,GM1,dvphiphj)  ! CCCC
            call CS_dim_fro2cls()
            call closed_shell_update&
            (MPSTDMs(jr,ir)%d,MPSTDMs(jr,ir)%P,TM2,GM2,dvphiphj)  ! CCCC
            deallocate(MPSTDMs(ir,jr)%d,MPSTDMs(ir,jr)%P)
            deallocate(MPSTDMs(jr,ir)%d,MPSTDMs(jr,ir)%P)
            allocate(MPSTDMs(ir,jr)%d(nocc,nocc))
            allocate(MPSTDMs(ir,jr)%P(nocc,nocc,nocc,nocc))
            allocate(MPSTDMs(jr,ir)%d(nocc,nocc))
            allocate(MPSTDMs(jr,ir)%P(nocc,nocc,nocc,nocc))
            MPSTDMs(ir,jr)%d=0.0d0
            MPSTDMs(ir,jr)%p=0.0d0
            MPSTDMs(jr,ir)%d=0.0d0
            MPSTDMs(jr,ir)%p=0.0d0
            MPSTDMs(ir,jr)%d=TM1
            MPSTDMs(ir,jr)%p=GM1
            MPSTDMs(jr,ir)%d=TM2
            MPSTDMs(jr,ir)%p=GM2
            write(17,*)"The closed shell updated 1-TDMs"
!            call print_mat1(nocc,nocc,MPSTDMs(ir,jr)%d,171)
!            call print_mat1(nocc,nocc,MPSTDMs(jr,ir)%d,172)
            call print_mat1(nocc,nocc,(MPSTDMs(ir,jr)%d+MPSTDMs(jr,ir)%d)/2.0d0,17)
            write(18,*)"The closed shell updated 2-TDMs"
!            call print_gat(nocc,nocc,nocc,nocc,MPSTDMs(ir,jr)%p,181)
!            call print_gat(nocc,nocc,nocc,nocc,MPSTDMs(jr,ir)%p,182)
!            call print_gat(nocc,nocc,nocc,nocc,(MPSTDMs(ir,jr)%p+MPSTDMs(jr,ir)%p)/2.0d0,18)
!            MPSTDMs(ir,jr)%d=(MPSTDMs(ir,jr)%d+MPS(iroot+1)%D)/2.0d0
!            MPSTDMs(jr,ir)%d=(MPSTDMs(jr,ir)%d+MPS(iroot+1)%D)/2.0d0
!            MPSTDMs(ir,jr)%P=(MPSTDMs(ir,jr)%P+MPS(iroot+1)%P)/2.0d0
!            MPSTDMs(jr,ir)%P=(MPSTDMs(jr,ir)%P+MPS(iroot+1)%P)/2.0d0
!            call print_mat(nocc,nocc,MPSTDMs(ir,jr)%d,2)
            deallocate(TM1,GM1)
            ! ====  Closed shell update ====
!            call flush(2)

! symmatrical TDMs  ! option-2 (both OK, option-1 is better; shift back later)
            MPSTDMs(ir,jr)%d=(MPSTDMs(ir,jr)%d+MPSTDMs(jr,ir)%d)/2.0d0
            MPSTDMs(jr,ir)%d=MPSTDMs(ir,jr)%d
            MPSTDMs(ir,jr)%p=(MPSTDMs(ir,jr)%p+MPSTDMs(jr,ir)%p)/2.0d0
            MPSTDMs(jr,ir)%p=MPSTDMs(ir,jr)%p
!            stop

            ! A-matrix can be constructed using a specific subroutine, then it can reduce many lines
            ! Cnstruct the A (gradients contribute) for transition
            ! one-body part
            allocate(MPSTDMs(ir,jr)%A(norb,norb))
            allocate(MPSTDMs(jr,ir)%A(norb,norb))
            MPSTDMs(ir,jr)%A=0.0d0
            MPSTDMs(jr,ir)%A=0.0d0

            ! A-matrix for state-transfer
            call fock_gen(orb%nsub,orb%act,orb%total,    &
                      nact,norb,T,U,MPSTDMs(ir,jr)%D,    &
                      MPSTDMs(ir,jr)%P,MPSTDMs(ir,jr)%A)

            MPSTDMs(jr,ir)%A=MPSTDMs(ir,jr)%A

            write(2,*)"(omit print)"
            write(2,*)"Gradient contribution for states transfering"
            write(2,*)"state",ir,jr
            call print_mat(norb,norb,MPSTDMs(ir,jr)%A,2)
            write(2,*)"state",jr,ir
            call print_mat(norb,norb,MPSTDMs(jr,ir)%A,2)
            call flush(2)

          end do
          MPSTDMs(ir,ir)%D=MPS(iroot+1)%D
          MPSTDMs(ir,ir)%P=MPS(iroot+1)%P
          MPSTDMs(ir,ir)%A=MPS(iroot+1)%A

!          write(2,*)"Finish the PreDMRGLR part"
!          call flush(2)

        end do

        write(2,*)"Finish the PreDMRGLR part"
        call flush(2)
        call CS_dim_fro2cls()

        deallocate(oRDMchk)
        deallocate(tRDMchk)
!        stop

      End Subroutine PreDMRGLR


