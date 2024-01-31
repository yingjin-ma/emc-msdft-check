      Subroutine LagrangeRDMs()

        use global_control
        use matrix

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
        double precision,allocatable::TM4(:,:),TM5(:,:),TM6(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)
        double precision,allocatable::GM3(:,:,:,:),GM4(:,:,:,:)

        double precision,allocatable::tmpP(:),tmpP2(:)
        double precision,allocatable::tmpT(:,:,:,:)

        integer o 

        character    ctmp1
        character*2  ctmp2
        character*192 string0,string1,string2,string3,string4

        double precision,allocatable::MPSci_lag(:)
        double precision Dscale,signMPS,signORB 

        itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

        ! treat closed-shell part as active
        call CS_dim_cls2fro()
        write(6,*)"Notice the closed-shell is treated as active"

        signMPS=1.0d0
        signORB=1.0d0

! to Leon :! update the contribution from multipliers into RDMs for a specific state

! 1) Orbital Lagrange contributations
!       (this part was omited since we can use Molcas to do the same thing)

        write(1,*)"The orbital Lagrange : "
        call print_mat(norb,norb,L_R%R,1)
        call flush(1)   

        do iroot=0,dmrg_nstates-1
!          write(6,*)ir,norb
!          call flush(6)
          ir=iroot+1
          ! allocate D,P (Lag.) for each states
          allocate(MPS(ir)%Dorb(norb,norb))           
          allocate(MPS(ir)%Porb(norb,norb,norb,norb))
          MPS(ir)%Dorb=0.0d0
          MPS(ir)%Porb=0.0d0

!      One-RDMs part

          allocate(TM1(norb,norb));  TM1=0.0d0 
     
          i0=0
          ia=0 
          do i=1,orb%nsub
            itot=orb%total(i)
            iact=orb%act(i)
            !write(6,*)itot,iact
            call flush(6) 
                    TM1(i0+1:i0+itot,i0+1:i0+itot) & 
            = MPS(ir)%D(ia+1:ia+iact,ia+1:ia+iact) 
            i0=i0+orb%total(i)
            ia=ia+orb%act(i)
          end do

!          write(1,*)"iroot",iroot,"DM"
!          call print_mat(norb,norb,TM1,1)
 
          ! First using the easy way
          if(.true.)then 
            do i=1,norb
              do j=1,norb
                do k=1,norb
                  MPS(ir)%Dorb(i,j)= &
                  MPS(ir)%Dorb(i,j)+TM1(k,j)*L_R%R(k,i)
                  MPS(ir)%Dorb(i,j)= &
                  MPS(ir)%Dorb(i,j)-TM1(i,k)*L_R%R(j,k)
                end do 
              end do 
            end do
          else
          ! Matrix muitiply & sub-irreps version later 
          end if

!          write(1,*)"iroot",iroot,"DM"
!          call print_mat(norb,norb,MPS(ir)%Dorb,1) 
 
          deallocate(TM1) 

!       Two RDM part

          allocate(GM1(norb,norb,norb,norb)) 
          GM1=0.0d0

          i0=0
          ia=0
          DO i=1,orb%nsub
            j0=0
            ja=0
            DO j=1,orb%nsub
              k0=0
              ka=0
              DO k=1,orb%nsub
                l0=0
                la=0
                DO l=1,orb%nsub
                  if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
                      GM1(i0+1:i0+orb%total(i),j0+1:j0+orb%total(j), &
                          k0+1:k0+orb%total(k),l0+1:l0+orb%total(l)) &
                    =MPS(ir)%P(ia+1:ia+orb%act(i),ja+1:ja+orb%act(j),&
                               ka+1:ka+orb%act(k),la+1:la+orb%act(l))
                  end if
                  l0=l0+orb%total(l)
                  la=la+orb%act(l)
                END DO
                k0=k0+orb%total(k)
                ka=ka+orb%act(k)
              END DO
              j0=j0+orb%total(j)
              ja=ja+orb%act(j)
            END DO
            i0=i0+orb%total(i)
            ia=ia+orb%act(i)
          END DO

!          call print_gat(nact,nact,nact,nact,MPS(ir)%P,222)
!          call print_gat(norb,norb,norb,norb,GM1,224)

          ! First using the easy way
          if(.true.)then
            do i=1,norb
              do j=1,norb
                do k=1,norb 
                  do l=1,norb 
                    do o=1,norb 
                        MPS(ir)%Porb(i,j,k,l)=MPS(ir)%Porb(i,j,k,l)+  &
                      GM1(o,j,k,l)*L_R%R(o,i)-GM1(i,o,k,l)*L_R%R(j,o) &
                     +GM1(i,j,o,l)*L_R%R(o,k)-GM1(i,j,k,o)*L_R%R(l,o) 
                    end do
                  end do
                end do
              end do
            end do
          end if

!          call print_gat(norb,norb,norb,norb,MPS(ir)%Porb,226) 
          deallocate(GM1) 

        end do        


! 2) MPSci Lagrange contributations
          
        do iroot=0,dmrg_nstates-1
 
          ir=iroot+1

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if

              nC=MPS(iroot+1)%nC
          nC_all=MPS(iroot+1)%nC_all

!          write(1,*)nc,nC_all,MPS(iroot+1)%nC,MPS(iroot+1)%nC_all
!          call flush(1)

! Allocate the space for the symmetric RDMs modification from Lagrange
          allocate(MPS(iroot+1)%Dmps(nact,nact))
          allocate(MPS(iroot+1)%pmps(nact,nact,nact,nact))
          MPS(iroot+1)%Dmps=0.0d0
          MPS(iroot+1)%Pmps=0.0d0
! Notice, may be strange that the RDMs-LR-L and RDMs-LR-R are initial symmetric
! Keep eyes on this point  -- YMA

          allocate(MPSci_lag(nc_all)) 
          MPSci_lag=0.0d0

!          write(1,*)MPSci_lag
!          call flush(1)

          Dscale=0.0d0
          do iC=1,nC
            !MPS(iroot+1)%Ci(iC)%lv
            write(1,*)"MPSci dv/Lagrange/index", &
                       MPS(iroot+1)%Ci(iC)%dv,   & 
                       MPS(iroot+1)%Ci(iC)%lv,   &
                       MPS(iroot+1)%Ci(iC)%num!  
            call flush(1)  
            MPSci_lag(MPS(ir)%Ci(iC)%num)=MPS(ir)%Ci(iC)%lv
            Dscale=Dscale+(MPS(ir)%Ci(iC)%lv)**2
          end do          

          write(1,*)"Origonal Dscale",Dscale
! Dscale is considered as re-scale the MPS (avoid bug if 2 sites) 

          call CS_dim_fro2cls() ! only to check the real active
          if(nact.le.2)then
            Dscale=dsqrt(Dscale)
          else  
            Dscale=1.0d0
          end if
          call CS_dim_cls2fro() ! Have to back to its old form

          if(LRsign) Dscale=-1.0d0*Dscale 

          write(1,*)"MPSci_lag",MPSci_lag
          write(1,*)"   Dscale",Dscale
          call flush(1)

! to Leon :   Yingjin-devel QCmaquis need to be used together with this part

          !string1="MPSci_lagrange"//trim(ctmp2)//".txt"
          !string1="MPSci_lagrange"//".txt"
          string1="updated_davidson.vec"

          write(1,*)"updated MPSci vectors for state ",iroot
          open(unit=100,file=trim(string1))
            do i=1,nC_all
              write(100,*)MPSci_lag(i)
              write(1,*)i,MPSci_lag(i)
            end do 
!             ==  Testing  ==
!              write(100,*) 0.01900536997
!              write(100,*) 0
!              write(100,*) 0 
!              write(100,*) 0.0690967511
!             ==  Testing  == 
          close(100)

          call MPSci_lagrange_update(iroot,"L")
          call MPSci_lagrange_update(iroot,"R")

        ! treat closed-shell-active back to closed-shell and active
          call CS_dim_fro2cls()
         
          call Symmetric_RDMs_lagrange(iroot,Dscale)

          ! 1'-RDMs part
          ! 2'-RDMs part         

          deallocate(MPSci_lag)

        end do        

!      ! treat "closed-shell and active" to closed-shell-active 
!        no need, since only active-space is needed
!        call CS_dim_cls2fro()
        call CS_dim_fro2cls()

        nDLMO=nact*(nact+1)/2
        nPLMO=nDLMO*(nDLMO+1)/2

        allocate(L_R%Dmps(nact,nact));           L_R%Dmps=0.0d0
        allocate(L_R%Pmps(nact,nact,nact,nact)); L_R%Pmps=0.0d0
        do iroot=0,dmrg_nstates-1         
          L_R%Dmps=L_R%Dmps+dmrg_weight(iroot+1)*MPS(iroot+1)%Dmps
          L_R%Pmps=L_R%Pmps+dmrg_weight(iroot+1)*MPS(iroot+1)%Pmps
        end do

        write(1,*)"Printed Dmps"
        call print_mat(nact,nact,L_R%Dmps,1) 
!        call print_gat(nact,nact,nact,nact,L_R%Pmps,11) 

        signMPS=1.0d0 
        if(signMPS.lt.0)then 
          Write(1,*)"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
          Write(1,*)" The sign of printed MPS-LAG(DM) will be changed"
          Write(6,*)"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
          Write(6,*)" The sign of printed MPS-LAG(DM) will be changed"
        end if
  
! Print RDMs-LAG in a compact form (only for C1 symmetry)
        open(unit=110,file="1RDM_lag.txt") 
          ij=0
          do i=1,nact
            do j=1,i
              ij=ij+1
              write(110,*)ij,L_R%Dmps(i,j)*signMPS 
            end do
          end do
        close(110)
        write(1,*)"Printed 1-mps-Lag" 
        call flush(1)

        allocate(tmpT(nact,nact,nact,nact)); tmpT=0.0d0
        do i=1,nact
          do j=1,nact
            do k=1,j
              do l=1,nact
                if(j.eq.k)then
                  tmpT(i,j,k,l)=L_R%Pmps(i,j,k,l)
                else
                  tmpT(i,j,k,l)=L_R%Pmps(i,j,k,l)+L_R%Pmps(i,k,j,l)
                end if
              end do
            end do
          end do
        end do
        write(1,*)"tmpT finish"
        call flush(1) 

        allocate(tmpP(nact**2*(nact**2+1)/2)); tmpP=0.0d0
        ij=0
        ijkl=1
        do i=1,nact
          do j=1,i
            ij=ij+1
            kl=0
            do k=1,nact
              do l=1,k
                kl=kl+1
                if(ij.ge.kl)then
                  write(1,*)ijkl,tmpT(i,k,l,j)
                  call flush(1) 
                  tmpP(ijkl)=tmpT(i,k,l,j)
                  ijkl=ijkl+1
                end if
              end do
            end do
          end do
        end do      
 

! to Leon :  All Molcas needed for the grad/geom-opt

        open(unit=110,file="2RDM_lag.txt") 
          do i=1,nPLMO
            write(110,*)i,tmpP(i)*signMPS
          end do
        close(110)

        allocate(tmpP2((nact**2+1)*nact**2/2)); tmpP2=0.0d0
        ijkl=0
        Do i=1,nAct
          Do j=1,nAct
            Do k=1,nAct
              Do l=1,nAct
                ij1=nAct*(i-1)+j
                ij2=nact*(j-1)+i
                kl1=nAct*(k-1)+l
                kl2=nact*(l-1)+k
                if(ij1.ge.kl1)then
                  ijkl=ijkl+1
                  tmpP2(itri(ij1,kl1))=L_R%Pmps(i,k,l,j)
!              write(1211,*)ijkl,itri(ij1,kl1),L_R%Pmps(i,k,l,j),i,k,l,j
                end if  
              End Do
            End Do
          End Do
        End Do

        open(unit=110,file="2RDM_lag_full.txt")
          ijkl=0
          do i=1,(nact**2+1)*nact**2/2
            ijkl=ijkl+1
            write(110,*)ijkl,2.0d0*tmpP2(i)*signMPS
          end do
        close(110)

        deallocate(tmpP,tmpP2)   
        deallocate(tmpT)

!        call CS_dim_fro2cls()

        ncls=nocc-nact
        nvir=norb-nocc

        write(6,*)" == The orb-lag is written into orb_lag.txt == "
        write(6,*)" ncls",ncls," nocc",nocc," norb",norb

        if(signMPS.lt.0)then 
          Write(1,*)"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
          Write(1,*)"   The sign of printed orb-LAG will be changed"
          Write(6,*)"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
          Write(6,*)"   The sign of printed orb-LAG will be changed"
        end if
 
        open(unit=110,file="orb_lag.txt") 
        open(unit=111,file="orb_lag_detail.txt") 
          ij=0
          do i=1,ncls
            do j=ncls+1,norb
              ij=ij+1
              write(111,*)ij,"i-p",i,j,L_R%R(i,j)
              write(110,*)ij,L_R%R(i,j)*signORB
            end do
          end do
          do i=ncls+1,nocc
            do j=nocc+1,norb
              ij=ij+1
              write(111,*)ij,"o-v",i,j,L_R%R(i,j)
              write(110,*)ij,L_R%R(i,j)*signORB
            end do
          end do
        close(111)
        close(110)

!          Do i=1,nact
!            Do j=1,i
!              ij=itri(i,j)
!              ij2=i+(j-1)*nact
!              ji2=j+(i-1)*nact
!              Do k=1,i
!                Do l=1,k
!                  kl=itri(k,l)
!                  kl2=k+(l-1)*nact
!                  lk2=l+(k-1)*nact
                  ! For symmetrizing if only <n|Op|0> is used
!                  ijkl=itri(ij2,kl2)
!                  jikl=itri(ji2,kl2)
!                  ijlk=itri(ij2,lk2)
!                  jilk=itri(ji2,lk2)
!                  write(110,*)ijkl,tmpP(ijkl)
!                End Do
!              End Do
!            End Do
!          End Do
!        close(110)


      End Subroutine LagrangeRDMs


      Subroutine Symmetric_RDMs_lagrange(ir,Dscale)

        use global_control
        use matrix

        integer :: ir
        double precision :: dscale
 
        double precision,allocatable::TM1(:,:),TM2(:,:)
        double precision,allocatable::GM1(:,:,:,:),GM2(:,:,:,:)

        character    ctmp1
        character*2  ctmp2
        character*192 string0,string1,string2,string3,string4

        ctmp2=""
        if(ir.ge.10)then
          write(ctmp2,"(I2)")ir
        else
          write(ctmp1,"(I1)")ir
          ctmp2(1:1)=ctmp1
        end if
 
        string1="oneRDM."//trim(ctmp2)//".LR-L"
        string2="oneRDM."//trim(ctmp2)//".LR-R"

        allocate(TM1(nact,nact)); TM1=0.0d0
        allocate(TM2(nact,nact)); TM2=0.0d0

        open(unit=110,file=trim(string1))
          read(110,*)ij
          do i=1,nact
            do j=1,nact
              read(110,*)ij,ji,TM1(i,j)
            end do
          end do
        close(110) 

        open(unit=110,file=trim(string2))
          read(110,*)ij
          do i=1,nact
            do j=1,nact
              read(110,*)ij,ji,TM2(i,j)
            end do
          end do
        close(110) 

!   symmetric and cumulatived
        MPS(ir+1)%Dmps=TM1+TM2
        write(1,*)"The symmetric Lag-RDMs-L : "
        call print_mat(nact,nact,TM1/1.0d0,1) 
        write(1,*)"The symmetric Lag-RDMs-R : "
        call print_mat(nact,nact,TM2/1.0d0,1) 
        MPS(ir+1)%Dmps=MPS(ir+1)%Dmps*Dscale  

        write(1,*)"Dmps and Dsacle",Dscale
        call print_mat(nact,nact,MPS(ir+1)%Dmps,1)

        string1="twoRDM."//trim(ctmp2)//".LR-L"
        string2="twoRDM."//trim(ctmp2)//".LR-R"

        allocate(GM1(nact,nact,nact,nact)); GM1=0.0d0
        allocate(GM2(nact,nact,nact,nact)); GM2=0.0d0

        open(unit=110,file=trim(string1))
           READ(110,*)ijkl
           ijkl=0  ! Just a counter
           do
             read(110,*,iostat=ierr)ij,jk,kl,li,dv
             if(ierr.ne.0)exit
             ijkl=ijkl+1
             GM1(ij+1,jk+1,kl+1,li+1)=dv
             GM1(kl+1,li+1,ij+1,jk+1)=dv
             GM1(jk+1,ij+1,li+1,kl+1)=dv
             GM1(li+1,kl+1,jk+1,ij+1)=dv
           end do       
        close(110)
 
        open(unit=110,file=trim(string2))
           READ(110,*)ijkl
           ijkl=0  ! Just a counter
           do
             read(110,*,iostat=ierr)ij,jk,kl,li,dv
             if(ierr.ne.0)exit
             ijkl=ijkl+1
             GM2(ij+1,jk+1,kl+1,li+1)=dv
             GM2(kl+1,li+1,ij+1,jk+1)=dv
             GM2(jk+1,ij+1,li+1,kl+1)=dv
             GM2(li+1,kl+1,jk+1,ij+1)=dv
           end do
        close(110)

!   symmetric and cumulatived
        MPS(ir+1)%Pmps=GM1+GM2
        MPS(ir+1)%Pmps=MPS(ir+1)%Pmps*Dscale 

        write(11,*)"Pmps"
        call print_gat(nact,nact,nact,nact,MPS(ir+1)%Pmps,11)

      End Subroutine Symmetric_RDMs_lagrange



