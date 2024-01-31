      Subroutine input()
      
        use global_control

! Ax: temporary char array
        character    A0
        character*72 A1
        character*10 A2
        character*10 A3 
        integer IL(8)
        double precision D1,D2
        logical        :: is_available
        character(200) :: line
        character(len=255)::ctemp

!        write(*,*)"start infile"

        ! initialized values
        act_omit=.true.
        old_code=.false.

        ! If (default) further SA-DMRG-SCF is needed for 
        !     Coupled-Perturbation/Lagrange calculations
        DMRG_SCF= .true.
        ! If (default) Lagrange multipliers for specific states
        DMRG_LAG=.false.
        ! DIIS acceleration
        !LDIIS   = .true.
        LDIIS   = .false.
        LDIIS_updated =.false.
        ! Interval(m) for DIIS
        ITERVAL_DIIS = 3
        ! MCSCF-type rotation
        LMCSCF  =.false.
        ! Number of TDMs that need to be solved
        ntdms=0
        ! By default assuming the same irreps for the states in DMRG-MCLR
        LRirreps = .true.

        Lcorrected_phase=.false. 

        ctemp=""
        Cworkdir=""! 
        ctemp=trim(Cworkdir)//"infile.inp"
        write(6,*)"ctemp",ctemp

       inquire(file=ctemp, exist=is_available)
       if(.not.is_available)then
         write(6,*) ' Input file infile.inp not available, ERROR'
         stop -1
       else
         write(6,'(/a)') ' Content of the input file'
         write(6,*) '-------------------------'
         open(100,                  &
              file   = 'infile.inp',&
              status = 'unknown',   &
              form   = 'formatted', &
              access = 'sequential')
         rewind(100)
         do while (.true.)
           read(100, '(a200)', end=1) line 
           write(6, '(a)') line(1:132)
         end do
 1       continue
         close(100, status = 'keep')
         write(6,'(/a//)') '-------------------------'
       end if

        open(unit=100,file=ctemp)
          do while (.true.)
            read(100,1000,iostat=io1)A1
            !write(*,*)A1(1:15),io1
            if(io1.ne.0) exit
            if(A1(1:7).eq."$ORBINT")then
              do while(.true.) 
                read(100,*)A2,A3
                !write(*,*)A2,A3
                if(A2(1:3).eq."sym")then
                  orb%sym=A3
                  if(orb%sym.eq."d2h")then 
                    orb%nsub=8 
                  else if(orb%sym.eq."c2v")then
                    orb%nsub=4
                  else if(orb%sym.eq."c1")then
                    orb%nsub=1             
                  end if 

                  allocate(orb%frozen(orb%nsub)); orb%frozen=0 !define
                  allocate(orb%closed(orb%nsub)); orb%closed=0
                  allocate(orb%   occ(orb%nsub)); orb%   occ=0
                  allocate(orb%   act(orb%nsub)); orb%   act=0
                  allocate(orb%  act2(orb%nsub)); orb%  act2=0
                  allocate(orb% extnl(orb%nsub)); orb% extnl=0
                  allocate(orb% total(orb%nsub)); orb% total=0
                  allocate(orb%   ele(orb%nsub)); orb%   ele=0

                  allocate(orbbk%frozen(orb%nsub)); orbbk%frozen=0 !define
                  allocate(orbbk%closed(orb%nsub)); orbbk%closed=0
                  allocate(orbbk%   occ(orb%nsub)); orbbk%   occ=0
                  allocate(orbbk%   act(orb%nsub)); orbbk%   act=0
                  allocate(orbbk%  act2(orb%nsub)); orbbk%  act2=0
                  allocate(orbbk% extnl(orb%nsub)); orbbk% extnl=0
                  allocate(orbbk% total(orb%nsub)); orbbk% total=0
                  allocate(orbbk%   ele(orb%nsub)); orbbk%   ele=0

                  allocate(orb%grouptable(orb%nsub,orb%nsub))
                  orb%grouptable=0
                  allocate(sym(orb%nsub))
                  sym%norb=0 

                  IL=0
                end if
 
                if(A2(1:8).eq."orbitals")FCIDUMP=A3
                if(A2(1:8).eq."electron")then
                  backspace(100)
                  read(100,*)A2,(orb%ele(i),i=1,orb%nsub)
                end if 
                if(A2.eq."$END")then
                  backspace(100)
                  exit 
                end if
              end do
              call group_table() 

! active space orbitals
            else if(A1(1:7).eq."$ORBSPC")then
              do while(.true.)
                read(100,*,iostat=io2)A2,(IL(i),i=1,orb%nsub)
                !open(unit=401,status="scratch")
                !  write()
                ! write(*,*)"A2",A2,IL(1),IL(2),IL(3),IL(4)IL(5)IL(6)IL(7)IL(8)
                if(io2.ne.0.or.A2(1:4).eq."$END") then
                  backspace(100)
                  exit
                end if 
                if(A2.eq."frozen")then
                  do i=1,orb%nsub
                    orb%frozen(i)=IL(i)  
                  end do
                else if(A2.eq."closed")then
                  !write(*,*)"closed"
                  do i=1,orb%nsub
                    orb%closed(i)=IL(i)  
                  end do
                else if(A2.eq."occ")then
                  do i=1,orb%nsub
                    orb%occ(i)=IL(i)
                  end do  
                end if
                A2=''
              end do

              orb%act=orb%occ-orb%frozen-orb%closed
              call FCIDUMP_READ()
!              write(*,*)"fcidump_after"
              !stop 
            else if(A1(1:11).eq."$DMRG_input")then

              dmrg_nstates=1
              dmrg_inputs=""
              dmrg_weight=0.0d0
              dmrg_rlxstate=1

              do 
                i=i+1
                read(100,"(A5)",iostat=ist)A2(1:5)
                if(ist.ne.0)exit 
                if(A2(1:5).eq."binar")then
!                 write(*,*)"Block of the input"
                  backspace(100) 
                  read(100,*)A2,dmrg_binary_name
                else if(A2(1:5).eq."input")then
                  backspace(100) 
!                 write(*,*)"Block of the input 1"
                  if(dmrg_nstates.eq.1)then
                    read(100,*)A2,dmrg_input_name                    
                    dmrg_inputs(1)=dmrg_input_name
                    dmrg_weight(1)=1.0d0
                  else
                    read(100,*)A2,(dmrg_inputs(i),i=1,dmrg_nstates)
                  end if
 
                  do i=1,dmrg_nstates
                    if(index(trim(dmrg_inputs(i)),trim(Cworkdir))&
                      .eq.0)then
                      dmrg_inputs(i)=trim(Cworkdir)&
                                   //trim(dmrg_inputs(i))
                    else
                      dmrg_inputs(i)=trim(dmrg_inputs(i))               
                    end if
                  end do
                  ! write(6,*) dmrg_nstates
                  ! write(6,*)"ctemp**********",dmrg_inputs(1)

                else if(A2(1:5).eq."tdm_i")then
                  backspace(100)                  
                  if(dmrg_nstates.eq.1)then
                  else 
                    do j=1,dmrg_nstates-1
                      ntdms=ntdms+j
                    end do
                    ntdms=ntdms*2  
                    read(100,*)A2,(dmrg_inputs_tdm(i),i=1,ntdms)
                    write(*,*)&
                    "Attention : "  
                    write(*,*)&
                    "Current input for TDMs is order sensitive"  

                    do i=1,ntdms
                      if(index(trim(dmrg_inputs_tdm(i)),&
                               trim(Cworkdir)).eq.0)then
                         dmrg_inputs_tdm(i)=&
                         trim(Cworkdir)//trim(dmrg_inputs_tdm(i))
                      else
                        dmrg_inputs_tdm(i)=trim(dmrg_inputs_tdm(i))
                      end if
                    end do

                  end if                   
!              write(*,*)"Block of the input 2"
                else if(A2(1:5).eq."reord")then
                  backspace(100) 
                  read(100,*)A2,dmrg_reorder_name
!              write(*,*)"Block of the input 3"
                else if(A2(1:5).eq."outpu")then
                  backspace(100) 
                  read(100,*)A2,dmrg_output_name
!              write(*,*)"Block of the input 4"
                else if(A2(1:5).eq."mpiru")then
                  backspace(100) 
                  read(100,*)A2,A3
                  if(A3.eq."yes".or.A3.eq."YES")then
                    ifmpirun=.true.
                  end if
                  read(100,*)A2,nproc_in_mpirun
                else if(A2(1:5).eq."nstat")then
                  backspace(100)
                  ! Total states and the state for relaxing
                  read(100,*)A2,dmrg_nstates
                else if(A2(1:5).eq."weigh")then
                  backspace(100)
                  read(100,*)A2,(dmrg_weight(i),i=1,dmrg_nstates)     
                end if
                if(A2(1:4).eq."$END")exit
              end do 

              d1=0.0d0 
              do i=1,dmrg_nstates
                d1=d1+dmrg_weight(i)
              end do
              dmrg_weight=dmrg_weight/d1

!              call  MODIFIED_GRAM_SCHMIDT(1,dmrg_nstates,dmrg_weight(1:dmrg_nstates))

!              write(*,*)"Block",dmrg_binary_name,dmrg_input_name,dmrg_rdm1_name,dmrg_rdm2_name
              !write(*,*)"Block",dmrg_output_name,ifmpirun,nproc_in_mpirun
            else if(A1(1:11).eq."$DMRG-Lagra")then

              CKP_LR=.true.
              dmrg_rlxstate=1
              signtdm=1.0d0
              thresMPSHami=1.0d-9 

              do 
                read(100,"(A7)",iostat=ist)A2(1:7)
                if(ist.ne.0)exit 

                if(A2(1:7).eq."DMRG-SC")then
                  backspace(100) 
                  read(100,*)A2,DMRG_SCF 
                else if(A2(1:7).eq."DMRG-LR")then 
                  backspace(100) 
                  read(100,*)A2,DMRG_LAG
                else if(A2(1:7).eq."ckp_re-")then 
                  backspace(100) 
                  read(100,*)A2,CKP_LR
                else if(A2(1:7).eq."rlxstat")then 
                  backspace(100) 
                  read(100,*)A2,dmrg_rlxstate
                else if(A2(1:6).eq."irreps")then 
                  backspace(100) 
                  read(100,*)A2,LRirreps
                else if(A2(1:7).eq."signtdm")then 
                  backspace(100) 
                  read(100,*)A2,signtdm
                else if(A2(1:7).eq."thresMP")then 
                  backspace(100) 
                  read(100,*)A2,thresMPSHami
                end if

                if(A2(1:4).eq."$END")exit
              end do 
             
              !write(6,*)DMRG_SCF,DMRG_LAG 
              
            else if(A1(1:12).eq."$DMRG-CASSCF")then 
              !write(*,*)"DMRG-CASSCF read in"
              ! The default orb-opt strategy

              ! initialize parameters 
              thrs%r=1.0d-8   ! rotation
              thrs%s=1.0d-8   ! symmetry
              thrs%e=1.0d-8   ! energy
              thrs%d=1.0d-12  ! davidson
              thrs%c=1.0d-4   ! coupling
              thrs%g=1.0d-2   ! gradients
              Nci_update=1000
              NH2_update=0
              step_damping=0.0d0
              restart_maquis=.false. 

              do i=1,99
                method(i)="direct4   "
              end do 
              i=0
              do 
                i=i+1
                read(100,"(A6)",iostat=ist)A2(1:6)
                if(ist.ne.0)exit 
                if(A2(1:6).eq."method")then
                  backspace(100) 
                  read(100,*)A2,method(i),nstep
                  do j=1,nstep-1
                    method(i+j)=method(i)
                  end do
                  i=i+nstep-1
                  !write(*,*)"000"
                else if(A2(1:6).eq."maxcyc")then 
                  backspace(100) 
                  read(100,*)A2,ncycle
                  method(i)=method(i+1)
                else if(A2(1:6).eq."CP_int")then 
                  CP_integrals=.true.
                  backspace(100) 
                  read(100,*)A2,ith_INTE                   
                else if(A2(1:6).eq."CP_Hes")then 
                  backspace(100) 
                  read(100,*)A2,ith_HESS                   
                else if(A2(1:6).eq."MPSci_")then 
                  backspace(100) 
                  read(100,*)A2,Nci_update
                else if(A2(1:6).eq."Step_d")then 
                  backspace(100) 
                  read(100,*)A2,step_damping
                else if(A2(1:6).eq."H2_upd")then 
                  backspace(100) 
                  read(100,*)A2,NH2_update
                else if(A2(1:6).eq."MAX_it")then 
                  backspace(100) 
                  read(100,*)A2,Max_iter
                else if(A2(1:6).eq."THRS_d")then 
                  backspace(100) 
                  read(100,*)A2,thrs%d
                else if(A2(1:6).eq."THRS_s")then 
                  backspace(100) 
                  read(100,*)A2,thrs%s
                else if(A2(1:6).eq."THRS_r")then 
                  backspace(100) 
                  read(100,*)A2,thrs%r
                else if(A2(1:6).eq."THRS_e")then 
                  backspace(100) 
                  read(100,*)A2,thrs%e
                else if(A2(1:6).eq."THRS_c")then 
                  backspace(100) 
                  read(100,*)A2,thrs%c
                else if(A2(1:6).eq."THRS_g")then 
                  backspace(100) 
                  read(100,*)A2,thrs%g
                else if(A2(1:6).eq."omit_a")then 
                      act_omit=.true.
                else if(A2(1:6).eq."restar")then 
                      restart_maquis=.true.
                else if(A2(1:5).eq."mcscf")then 
                        LMCSCF=.true.
                else if(A2(1:6).eq."old_co")then
                      old_code=.true.
                end if
                ! write(6,*)"ctemp ",method
              end do
              !write(*,*)"DMRG-CASSCF",ncycle
              ! write(6,*)"ctemp ",i
            else if(A1(1:12).eq."$SCF-CONTROL")then

            end if 
          end do
        close(100) 
1000    format(A72)


        !write(*,*)act_omit
        !stop

        !write(*,*)"method",method(1:ncycle)
!        write(*,*)"DMRG-CASSCF ncycle ",ncycle
!        call flush(6)

!        do i=1,norb
!          write(*,*)U(1,1,i,1)
!        end do 
!        stop
        !int
        !integer

        ncore=0
        do i=1,orb%nsub
          ncore=ncore+orb%closed(i)
          ncore=ncore+orb%frozen(i)
        end do    
        if(ncore.gt.0)then
          ifcore=.true.
        else
          ifcore=.false.
        end if   
        ifcore2=ifcore

!        do i=1,orb%nsub
          !write(*,*)orb%extnl(i)
!        end do

        call integrals_active("FCIDUMP_ACTIVE    ",.false.)
!        write(*,*)"After FCIDUMP_ACTIVE,FRONREA",FRONREA

      end Subroutine input

! ---------------------------------------------------------------------
! Read in the active RDMs, not including the closed shell part
      Subroutine RDMs_READ()

        use global_control
        use matrix

        character    ctmp1
        character*2  ctmp2
        character*72 oneRDMfile,twoRDMfile

        double precision dv 
        double precision,allocatable::TM1(:,:)
        double precision,allocatable::GM1(:,:,:,:)

        ! 1-RDMs
        if(.not.allocated(mat2%d))then
          allocate(mat2%d(nact,nact))           ; mat2%d=0.0d0
        else
          deallocate(mat2%d)
          allocate(mat2%d(nact,nact))  
          mat2%d=0.0d0
        end if
        allocate(TM1(nact,nact))              ;    TM1=0.0d0
        ! 2-RDMs
        if(.not.allocated(mat2%p))then
          allocate(mat2%p(nact,nact,nact,nact)) ; mat2%p=0.0d0
        else
          deallocate(mat2%p)
          allocate(mat2%p(nact,nact,nact,nact))
          mat2%p=0.0d0
        end if

        allocate(GM1(nact,nact,nact,nact))    ;    GM1=0.0d0

        !        write(6,*)" iorder ",iorder; call flush(6)

        do iroot=0,dmrg_nstates-1

          !          write(*,*)" iroot ",iroot,"RDMs importing"; call flush(6)

                    ctmp2=""
                    if(iroot.ge.10)then
                      write(ctmp2,"(I2)")iroot
                    else
                      write(ctmp1,"(I1)")iroot
                      ctmp2(1:1)=ctmp1
                    end if

                    if(dmrg_binary_name.eq."block")then
                      oneRDMfile="spatial_onepdm."//trim(ctmp2)//"."//trim(ctmp2)
                      twoRDMfile="spatial_twopdm."//trim(ctmp2)//"."//trim(ctmp2)
                    else if(dmrg_binary_name.eq."mrci")then
                      oneRDMfile="1rdm.out"
                      twoRDMfile="2rdm.out"
                    else
                      oneRDMfile="oneRDM."//trim(ctmp2)//"."//trim(ctmp2)
                      twoRDMfile="twoRDM."//trim(ctmp2)//"."//trim(ctmp2)
                    end if

          !          oneRDMfile=trim(Cworkdir)//trim(oneRDMfile)
          !          twoRDMfile=trim(Cworkdir)//trim(twoRDMfile)

        !          write(6,*)"iorder",iorder

                  mat2%d=0.0d0
                  open(unit=110,file=oneRDMfile)
        !           read(110,*)ij
                    do i=1,nact
                      do j=1,nact
                        !                write(6,*)"Import 1-RDMs"
                        read(110,*)ij,ji,mat2%d(iorder(i),iorder(j))
                        !                write(6,*)"Imported 1-RDMs", mat2%d(iorder(i),iorder(j))
                        !                call flush(6)
                        end do
                      end do
                    close(110)
                    TM1=TM1+mat2%d*dmrg_weight(iroot+1)

          !          write(6,*)"Imported 1-RDMs"; call flush(6)

          mat2%p=0.0d0
          open(unit=110,file=twoRDMfile)
          !          READ(110,*)ijkl
          if(dmrg_binary_name.eq."maquis")then
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
        !              write(*,*)ij+1,jk+1,kl+1,li+1,rdm2P(ij+1,jk+1,kl+1,li+1)
            end do
          else if (dmrg_binary_name.eq."mrci")then 
            ijkl=0  ! Just a counter
            do
              read(110,*,iostat=ierr)ij,jk,kl,li,dv
              if(ierr.ne.0)exit
              ijkl=ijkl+1
              mat2%P(iorder(ij),&
                     iorder(jk),&
                     iorder(kl),&
                     iorder(li))=dv
              mat2%P(iorder(kl),&
                     iorder(li),&
                     iorder(ij),&
                     iorder(jk))=dv
              mat2%P(iorder(jk),&
                     iorder(ij),&
                     iorder(li),&
                     iorder(kl))=dv
              mat2%P(iorder(li),&
                     iorder(kl),&
                     iorder(jk),&
                     iorder(ij))=dv
          !              write(*,*)ij,jk,kl,li,rdm2P(ij+1,jk+1,kl+1,li+1)
            end do                  
          else
            do i=1,nact
              do j=1,nact
                do k=1,nact
                  do l=1,nact
                    READ(110,*)ij,jk,kl,li,&
                    mat2%p(iorder(i),iorder(j),iorder(k),iorder(l))
                    !mat2%p(iorder(i),iorder(j),iorder(k),iorder(l))
                  end do
                end do
              end do
            end do
          end if
          close(110)
          GM1=GM1+mat2%p*dmrg_weight(iroot+1)

          write(31,*)"  2-RDMs for state", iroot
        !          write(*,*)"   2-RDMs for state"; call flush(6)
          call print_gat(nact,nact,nact,nact,mat2%P,31)

        end do

        ! return the weighted RDMs back 
        mat2%d=TM1
        mat2%p=GM1
        deallocate(TM1,GM1)

        ! In order to match the energy expression containing the "1/2"
        mat2%p=mat2%p*2.0d0

        ! Not needed
        !        call symmetrize()

      end Subroutine RDMs_READ
        ! ---------------------------------------------------------------------

      Subroutine FCIDUMP_READ_MATRIX()

        use global_control

        character*10 A2,A3
        integer,allocatable::IL(:)
        double precision D1,D2

        !null(head%dv)
        !null(ptr)
        !null(tail)
        nullify(head,ptr,tail)
        if(index(trim(FCIDUMP),trim(Cworkdir)).eq.0)then
          FCIDUMP=trim(Cworkdir)//trim(FCIDUMP)
        else
          FCIDUMP=trim(FCIDUMP)                
        end if
! total orbitals
              open(unit=101,file=FCIDUMP)
                read(101,*)A2,A3,norb 
                allocate(IL(norb))                         !define  
                read(101,"(A9)",advance='no')A2
                read(101,*)(IL(i),i=1,norb)
                do i=1,orb%nsub
                  sym(i)%norb=0
                end do
                do i=1,norb
                  do j=1,orb%nsub
                    if(IL(i).eq.j)then
                      sym(j)%norb=sym(j)%norb+1
                    end if
                  end do
                end do
                do i=1,orb%nsub
                  orb%total(i)=sym(i)%norb              
                end do
                orb%extnl=orb%total-orb%occ
                orb%act=orb%occ-orb%closed  ! Since it is directly from FCIDUMP file, 
                                            !       the frozens are always 0
                nact=0
                nocc=0
                do i=1,orb%nsub
                  nact=nact+orb%act(i)
                  nocc=nocc+orb%occ(i)
                end do 
                deallocate(IL)                             !dedefine                 

! Read the integrals & distribute integrals
              ! working integrals
              allocate(T(norb,norb)) ; T=0.0d0            !define
              allocate(U(norb,norb,norb,norb)) ; U=0.0d0
              ! origionl in each macro-iter
              allocate(T_origin(norb,norb)) ; T_origin=0.0d0            !define
              allocate(U_origin(norb,norb,norb,norb)) ; U_origin=0.0d0
              ! 2-index transform with full dim              
              allocate(T2H(norb,norb)); T2H=0.0d0         !define
              allocate(U2H(norb,norb,norb,norb)); U2H=0.0d0

! The same time, allocate the 1-index transformed Integrals
              allocate(T1dx(norb,norb)) ; T1dx=0.0d0      
              allocate(J1dx(norb,norb,norb,norb)) ; J1dx=0.0d0
              allocate(K1dx(norb,norb,norb,norb)) ; K1dx=0.0d0
! The same time, allocate the 2-index transformed Integrals
              allocate(T2nd(nocc,nocc)) ; T2nd=0.0d0      
              allocate(J2dx(norb,norb,norb,norb)) ; J2dx=0.0d0
              allocate(K2dx(norb,norb,norb,norb)) ; K2dx=0.0d0
              allocate(U2nd(nocc,nocc,nocc,nocc)) ; U2nd=0.0d0

                n1e=0
                n2e=0
                do while(.true.) 
                  read(101,*,iostat=io3)A2
                  if(io3.ne.0) exit
                  if(A2.eq."&END")then
                    ival=0
                    do while(.true.)
                      read(101,*,iostat=io4)d1,i1,i2,i3,i4
                      if(io4.ne.0) exit
                      ival=ival+1
                      if(i3.eq.0)then
                        if(i1.ne.0)then
                          n1e=n1e+1
                          T(i1,i2)=d1
                          T(i2,i1)=d1
                        else
                          FRONRE0=d1
                        end if 
                      else
                        n2e=n2e+1
                        U(i1,i2,i3,i4)=d1
                        U(i1,i2,i4,i3)=d1
                        U(i2,i1,i3,i4)=d1
                        U(i2,i1,i4,i3)=d1
                        U(i3,i4,i1,i2)=d1
                        U(i3,i4,i2,i1)=d1
                        U(i4,i3,i1,i2)=d1
                        U(i4,i3,i2,i1)=d1
                      end if
                    end do   
                  end if
                end do
              close(101)

              T_origin=T
              U_origin=U
              T2H=T
              U2H=U

              FRONRE=FRONRE0

      end Subroutine FCIDUMP_READ_MATRIX

      Subroutine FCIDUMP_READ()

        use global_control

        character*10 A2,A3
        integer,allocatable::IL(:)
        double precision D1,D2

        !null(head%dv)
        !null(ptr)
        !null(tail)
        nullify(head,ptr,tail)
        if(index(trim(FCIDUMP),trim(Cworkdir)).eq.0)then
          FCIDUMP=trim(Cworkdir)//trim(FCIDUMP) 
        else
          FCIDUMP=trim(FCIDUMP) 
        end if
! total orbitals
              open(unit=101,file=FCIDUMP)
              !open(unit=121,file='dump.tmp')
                read(101,*)A2,A3,norb 
                ! write(121,*)A2,A3,norb 
                allocate(IL(norb))                         !define  
                ! write(*,*)"fsdgdsgvcsdgchjsdcjhshjc",norb
                read(101,"(A9)",advance='no')A2
                !write(121,"(A9)",advance='no')A2
                read(101,*)(IL(i),i=1,norb)
                write(6,*)"UUUUUUUUUUUUUUUUUUUUUu",norb
                ! do i=1,norb
                !  write(121,"(I2A1)",advance='no')IL(i),","
                !  write(*,*)"fsdgdsgvcsdgchjsdcjhshjc",IL(i)
                ! end do
                !write(121,*)
                !write(121,*)"&END"
                do i=1,orb%nsub
                  sym(i)%norb=0
                end do
                do i=1,norb
                  do j=1,orb%nsub
                    if(IL(i).eq.j)then
                      sym(j)%norb=sym(j)%norb+1
                    end if
                  end do
                end do
                ! write(*,*)sym%norb
                do i=1,orb%nsub
                  orb%total(i)=sym(i)%norb              
                end do
                orb%extnl=orb%total-orb%occ
                orb%act=orb%occ-orb%closed  ! Since it is directly from FCIDUMP file, 
                                            !       the frozens are always 0
                nact=0
                nocc=0
                do i=1,orb%nsub
                  nact=nact+orb%act(i)
                  nocc=nocc+orb%occ(i)
                end do 
                deallocate(IL)                             !dedefine 

! Read the integrals
                do while(.true.) 
                  read(101,*,iostat=io3)A2
                  if(io3.ne.0) exit
                  if(A2.eq."&END")then
                    !write(*,*)"&end"
                    ival=0
                    do while(.true.)
                      read(101,*,iostat=io4)d1,i1,i2,i3,i4
                      !write(*,*)d1,i1,i2,i3,i4
                      !stop
                      if(io4.ne.0) exit
                      ival=ival+1
                      if(.not.associated(head))then
                        allocate(head,STAT=istat)
                        tail => head
                        nullify(tail%p)
                        tail%i=i1
                        tail%j=i2
                        tail%k=i3
                        tail%l=i4
                        tail%dv=d1
                        !deallocate(head)
                      else
                        allocate(tail%p,STAT=istat)
                        tail => tail%p
                        nullify(tail%p)
                        tail%i=i1
                        tail%j=i2
                        tail%k=i3
                        tail%l=i4
                        tail%dv=d1
                        !deallocate(tail%p)
                      end if
                      !deallocate(head)
                      !deallocate(tail)
                    end do   
                  end if
                end do
              close(101)

! distribute integrals
              ! working integrals
              allocate(T(norb,norb)) ; T=0.0d0            !define
              allocate(U(norb,norb,norb,norb)) ; U=0.0d0
              write(6,*)"UUUUUUUUUUUUUUUUUUUUUu",norb
              U(7,7,1,1)=1.87898
              write(6,*)"UUUUUUUUUUUUUUUUUUUUUu",U(7,7,1,1)
              ! origionl in each macto-iter
              allocate(T_origin(norb,norb)) ; T_origin=0.0d0     !define
              allocate(U_origin(norb,norb,norb,norb)) ; U_origin=0.0d0
              ! integrals full dim
              allocate(T2H(norb,norb)); T2H=0.0d0         !define
              allocate(U2H(norb,norb,norb,norb)); U2H=0.0d0

! The same time, allocate the 1-index transformed Integrals
              allocate(T1dx(norb,norb)) ; T1dx=0.0d0      
              allocate(J1dx(norb,norb,norb,norb)) ; J1dx=0.0d0
              allocate(K1dx(norb,norb,norb,norb)) ; K1dx=0.0d0
! The same time, allocate the 2-index transformed Integrals
              allocate(T2nd(nocc,nocc)) ; T2nd=0.0d0      
              allocate(J2dx(norb,norb,norb,norb)) ; J2dx=0.0d0
              allocate(K2dx(norb,norb,norb,norb)) ; K2dx=0.0d0
              allocate(U2nd(nocc,nocc,nocc,nocc)) ; U2nd=0.0d0

              n1e=0
              n2e=0
              ptr => head
              !nullify(ptr%p)
              do
                if(.not.associated(ptr)) exit
                if(ptr%k.eq.0)then
                  if(ptr%i.ne.0)then
                    n1e=n1e+1
                    T(ptr%i,ptr%j)=ptr%dv
                    T(ptr%j,ptr%i)=ptr%dv
                  else
                    FRONRE0=ptr%dv
                  end if 
                else
                  n2e=n2e+1
! 8-folder symmetry
!       (i,j|k,l),(i,j|l,k),(j,i|k,l),(j,i|l,k),
!       (k,l|i,j),(k,l|j,i),(l,k|i,j),(l,k|j,i) 
                  !D2=ptr%dv
                  !if(dabs(D2).lt.1.0e-14)D2=0.0d0
                  U(ptr%i,ptr%j,ptr%k,ptr%l)=ptr%dv
                  U(ptr%i,ptr%j,ptr%l,ptr%k)=ptr%dv
                  U(ptr%j,ptr%i,ptr%k,ptr%l)=ptr%dv
                  U(ptr%j,ptr%i,ptr%l,ptr%k)=ptr%dv
                  U(ptr%k,ptr%l,ptr%i,ptr%j)=ptr%dv
                  U(ptr%k,ptr%l,ptr%j,ptr%i)=ptr%dv
                  U(ptr%l,ptr%k,ptr%i,ptr%j)=ptr%dv
                  U(ptr%l,ptr%k,ptr%j,ptr%i)=ptr%dv
                end if
                !write(*,*)ptr%dv
                ptr => ptr%p
                !nullify(ptr%p)
              end do

              T_origin=T
              U_origin=U
              T2H=T
              U2H=U
 
              FRONRE=FRONRE0

!              write(*,*)"In FCIDUMP read - test for molcas"
!              ijkl=0
!              do i=1,4
!                do j=1,i
!                  do k=1,4
!                    do l=1,7
!                      ijkl=ijkl+1
!                      write(432,*)i,j,k,l,U(l,k,j,i),ijkl  
!                    end do
!                  end do
!                end do
!              end do

!              stop 
        nullify(head,ptr,tail)

      end subroutine FCIDUMP_READ


      Subroutine FCIDUMP_READ_LR()

        use global_control

        character*10 A2,A3
        integer,allocatable::IL(:)
        double precision D1,D2

        nullify(head,ptr,tail)
        if(index(trim(FCIDUMP),trim(Cworkdir)).eq.0)then
          FCIDUMP=trim(Cworkdir)//trim(FCIDUMP)
        else
          FCIDUMP=trim(FCIDUMP)                
        end if
! total orbitals
              open(unit=101,file=FCIDUMP)
                read(101,*)A2,A3,norb 
                allocate(IL(norb))                         !define  
                read(101,"(A9)",advance='no')A2
                read(101,*)(IL(i),i=1,norb)
                do i=1,orb%nsub
                  sym(i)%norb=0
                end do
                do i=1,norb
                  do j=1,orb%nsub
                    if(IL(i).eq.j)then
                      sym(j)%norb=sym(j)%norb+1
                    end if
                  end do
                end do
                do i=1,orb%nsub
                  orb%total(i)=sym(i)%norb              
                end do
                orb%extnl=orb%total-orb%occ
                orb%act=orb%occ-orb%closed  ! Since it is directly from FCIDUMP file, 
                                            !       the frozens are always 0
                nact=0
                nocc=0
                do i=1,orb%nsub
                  nact=nact+orb%act(i)
                  nocc=nocc+orb%occ(i)
                end do 
                deallocate(IL)                             !dedefine 

! Read the integrals
                do while(.true.) 
                  read(101,*,iostat=io3)A2
                  if(io3.ne.0) exit
                  if(A2.eq."&END")then
                    ival=0
                    do while(.true.)
                      read(101,*,iostat=io4)d1,i1,i2,i3,i4
                      if(io4.ne.0) exit
                      ival=ival+1
                      if(.not.associated(head))then
                        allocate(head,STAT=istat)
                        tail => head
                        nullify(tail%p)
                        tail%i=i1
                        tail%j=i2
                        tail%k=i3
                        tail%l=i4
                        tail%dv=d1
                      else
                        allocate(tail%p,STAT=istat)
                        tail => tail%p
                        nullify(tail%p)
                        tail%i=i1
                        tail%j=i2
                        tail%k=i3
                        tail%l=i4
                        tail%dv=d1
                      end if
                    end do   
                  end if
                end do
              close(101)

! distribute integrals
              ! working integrals
              if(.not.allocated(T))then
                allocate(T(norb,norb)) ; T=0.0d0   !define
              else
                T=0.0d0
              end if
              if(.not.allocated(U))then
                allocate(U(norb,norb,norb,norb)) ; U=0.0d0
              else
                U=0.0d0
              end if 
              ! origionl in each macto-iter (not change in LR)
              if(.not.allocated(T_origin))then
                allocate(T_origin(norb,norb)) ; T_origin=0.0d0   !define
              else
                T_origin=0.0d0
              end if 
              if(.not.allocated(U_origin))then 
                allocate(U_origin(norb,norb,norb,norb)) ; U_origin=0.0d0
              else
                U_origin=0.0d0
              end if             

              n1e=0
              n2e=0
              ptr => head
              !nullify(ptr%p)
              do
                if(.not.associated(ptr)) exit
                if(ptr%k.eq.0)then
                  if(ptr%i.ne.0)then
                    n1e=n1e+1
                    T(ptr%i,ptr%j)=ptr%dv
                    T(ptr%j,ptr%i)=ptr%dv
                  else
                    FRONRE0=ptr%dv
                  end if 
                else
                  n2e=n2e+1
! 8-folder symmetry
!       (i,j|k,l),(i,j|l,k),(j,i|k,l),(j,i|l,k),
!       (k,l|i,j),(k,l|j,i),(l,k|i,j),(l,k|j,i) 
                  !D2=ptr%dv
                  !if(dabs(D2).lt.1.0e-14)D2=0.0d0
                  U(ptr%i,ptr%j,ptr%k,ptr%l)=ptr%dv
                  U(ptr%i,ptr%j,ptr%l,ptr%k)=ptr%dv
                  U(ptr%j,ptr%i,ptr%k,ptr%l)=ptr%dv
                  U(ptr%j,ptr%i,ptr%l,ptr%k)=ptr%dv
                  U(ptr%k,ptr%l,ptr%i,ptr%j)=ptr%dv
                  U(ptr%k,ptr%l,ptr%j,ptr%i)=ptr%dv
                  U(ptr%l,ptr%k,ptr%i,ptr%j)=ptr%dv
                  U(ptr%l,ptr%k,ptr%j,ptr%i)=ptr%dv
                end if
                ptr => ptr%p
              end do

              T_origin=T
              U_origin=U
              FRONRE=FRONRE0

              nullify(head,ptr,tail) ! 

      end subroutine FCIDUMP_READ_LR

