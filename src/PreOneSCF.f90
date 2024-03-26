      Subroutine DMRG_CASSCF()
  
        use global_control
        use matrix

! Read in the orbital informations
        ! write(6,*)nact
! order in FCIDUMP
        allocate(iorder(nact)) ; iorder=0
        open(unit=110,file=dmrg_reorder_name)
          read(110,*)(iorder(i),i=1,nact)
!          write(*,*)(iorder(i),i=1,nact)
        close(110)
        
!write(*,*)dmrg_reorder_name
!stop
        write(*,*)""
        write(*,*)"DMRG-SCF Part (orbital optimization):"
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

      Subroutine PreOneSCF(icycle)
  
        use global_control
        use matrix
        use date_time

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

        !write(*,*)"In Preone SCF ",iorder,nact
!        stop

!  ==========================================================================
!  ==============     Active the closed shell approximatation    ============
!  ==============        Modifying the RDMs (can be omited)      ============
!  ==========================================================================
goto 2222

        if(ifcore) then
          allocate(TM1(nact,nact))               ;  TM1=mat2%d
          allocate(GM1(nact,nact,nact,nact))     ;  GM1=mat2%p
          deallocate(mat2%d)          
          deallocate(mat2%p)         
          allocate(mat2%d(nocc,nocc))            ; mat2%d=0.0d0
          allocate(mat2%p(nocc,nocc,nocc,nocc))  ; mat2%p=0.0d0
          ioffset=0
          ioffset1=0
          do i=1,orb%nsub
            do j=1,orb%closed(i)
              mat2%d(j+ioffset,j+ioffset)=2.0d0
            end do
            j0=orb%closed(i)+ioffset   ! offset for occupy orbitals in each irreps
            j1=ioffset1                ! offset for active orbitals in each irreps
            mat2%d(j0+1:j0+orb%act(i),j0+1:j0+orb%act(i))=&
               TM1(j1+1:j1+orb%act(i),j1+1:j1+orb%act(i))
            ioffset=ioffset+orb%occ(i)
            ioffset1=ioffset1+orb%act(i)
          end do
                 

          !          stop
          allocate(TM2(nocc,nocc))  ; TM2=mat2%d         

          i0=0 ! for the occ
          i1=0 ! for the act ! to be continue
          do i=1,orb%nsub 
           j0=0
           j1=0 
           do j=1,orb%nsub 
            k0=0
            k1=0 
            do k=1,orb%nsub 
             l0=0
             l1=0 
             do l=1,orb%nsub 
              if(orb%grouptable(i,orb%grouptable(j,orb%grouptable(k,l))).eq.1)then
          ! This is the closed shell part in orbital induce i 
                do ii=1,orb%closed(i)
                 do jj=1,orb%closed(j)
                  do kk=1,orb%occ(k)
                   do ll=1,orb%occ(l)
                     io=ii+i0 
                     jo=jj+j0 
                     ko=kk+k0 
                     lo=ll+l0                   
                     call delt(io,jo,dij) 
                     call delt(io,ko,dik) 
                     call delt(io,lo,dil) 
                     call delt(jo,lo,djl) 
                     mat2%p(io,ko,lo,jo)=dij*TM2(ko,lo)&
                                        -dil*TM2(jo,ko)
                     mat2%p(io,ko,lo,jo)=mat2%p(io,ko,lo,jo)*2.0d0

                     mat2%p(ko,io,jo,lo)=dij*TM2(ko,lo)&
                                        -dil*TM2(jo,ko)
                     mat2%p(ko,io,jo,lo)=mat2%p(ko,io,jo,lo)*2.0d0
                   end do
                  end do
                 end do                  
                end do
                do ii=1,orb%closed(i)
                 do jj=1,orb%closed(j)
                  do kk=1,orb%occ(k)
                   do ll=1,orb%occ(l)
                     io=ii+i0 
                     jo=jj+j0 
                     ko=kk+k0 
                     lo=ll+l0                   
                     call delt(io,jo,dij) 
                     call delt(io,ko,dik) 
                     call delt(io,lo,dil) 
                     call delt(jo,lo,djl)               

                     mat2%p(lo,io,ko,jo)=2.0d0*dik*djl&
                                            -0.5d0*dij*TM2(ko,lo)
                     mat2%p(lo,io,ko,jo)=mat2%p(lo,io,ko,jo)*2.0d0

                     mat2%p(io,lo,jo,ko)=2.0d0*dik*djl&
                                            -0.5d0*dij*TM2(ko,lo)
                     mat2%p(io,lo,jo,ko)=mat2%p(io,lo,jo,ko)*2.0d0

          !                     mat2%p(io,lo,jo,ko)=TM2(io,jo)*TM2(ko,lo)&
          !             -0.25d0*TM2(io,ko)*TM2(jo,lo)-0.25d0*TM2(io,lo)*TM2(jo,ko)
          !                     write(*,*)"io,jo,ko,lo",io,jo,ko,lo,mat2%p(io,jo,ko,lo) 
          !                     write(*,*)"io,jo,ko,lo",io,jo,ko,lo,mat2%p(io,jo,ko,lo) 
                   end do
                  end do
                 end do                  
                end do
          ! This is the closed shell part
                do ii=1,orb%act(i)
                 do jj=1,orb%act(j)
                  do kk=1,orb%act(k)
                   do ll=1,orb%act(l)
                     ix=ii+i1 
                     jx=jj+j1
                     kx=kk+k1
                     lx=ll+l1
                     io=ii+i0+orb%closed(i) 
                     jo=jj+j0+orb%closed(j) 
                     ko=kk+k0+orb%closed(k) 
                     lo=ll+l0+orb%closed(l)      
          !                     write(*,*)"io,jo,ko,lo",io,jo,ko,ko   
                     mat2%p(io,jo,ko,lo)=GM1(ix,jx,kx,lx)  !!! This is correct!!
                   end do
                  end do
                 end do                  
                end do
              end if 
              l0=l0+orb%occ(l)
              l1=l1+orb%act(l)
             end do
             k0=k0+orb%occ(k)
             k1=k1+orb%act(k)
            end do
            j0=j0+orb%occ(j)
            j1=j1+orb%act(j)
           end do
           i0=i0+orb%occ(i)
           i1=i1+orb%act(i)
          end do

          !          call print_GAT(nocc,nocc,nocc,nocc,mat2%p,11)
          !          stop 

          deallocate(TM1)
          deallocate(TM2)
          deallocate(GM1)
          !  The simply way to active closed shell approximation
          ! Commented the following lines on Apr. 12, 2015 
          !          nact2=nact
          !          nact=nocc  
          !          orb%act2=orb%act 
          !          orb%act=orb%occ 
          !          ifcore=.false.
          !          ifcore2=.true.
        end if

2222    continue

        allocate(mat1(orb%nsub)) 
        !mat1=0.0d0;
        j0=0; k0=0
        do i=1,orb%nsub
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
        write(*,*)"**************************"
        write(*,*)mat2%d(1,:)
        write(*,*) size(mat2%d),size(mat2%d, dim=1),size(mat2%d, dim=2)
        write(*,*) orb%act(1)

        write(*,*) size(mat2%p),size(mat2%p, dim=1),size(mat2%P, dim=2)
        write(*,*) nocc
 
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
          !          allocate(mat1(i)%JP_KL(orb%total(i),orb%total(i)))
          mat1(i)%A=0.0d0
          mat1(i)%B=0.0d0
          mat1(i)%hD=0.0d0
          !          mat1(i)%JP=0.0d0
        end do

!        write(6,*)" === "
!        call flush(6)
!        stop

        RDM_ENERGY=0.0d0      ! initial version
!        RDM_ENERGY_MAT=0.0d0  ! MKL (dgemm) version
        if(.not. ifcore) then
          if(old_code)then ! === old code ===
            ioffset=0
            do i=1,orb%nsub
              do j=1,orb%occ(i)
                koffset=0
                do k=1,orb%nsub
                  if(i.eq.k)then
                    do l=1,orb%occ(k) 
                      RDM_ENERGY=RDM_ENERGY+T(ioffset+j,koffset+l)*mat1(i)%D(j,l)
                    ! write(*,*)j,l,T(ioffset+j,koffset+l),mat1(i)%D(j,l)
                    end do
                  end if
                  koffset=koffset+orb%total(k)
                end do
              end do
              ioffset=ioffset+orb%total(i)
            end do
          else  ! === new code ===
            ioffset=0
            do i=1,orb%nsub
              N=orb%occ(i)
              if(N.ne.0)then
                M=ioffset 
                call MXM(N,T(M+1:M+N,M+1:M+N),mat1(i)%D,mat1(i)%hD(1:N,1:N))
                call trace(N,mat1(i)%hD(1:N,1:N),dtmp)
                RDM_ENERGY=RDM_ENERGY+dtmp
              end if 
            ioffset=ioffset+orb%total(i)
            end do
          end if ! === 
          if(old_code)then !== old code ===
            ioffset=0
            ioffset1=0
            !          RDM_IJKL=0.0d0
                        do i=1,orb%nsub; do i1=1,orb%occ(i)
                        !dtemp=0.0d0
                        !write(*,*)"RDM2-I",RDM_ENERGY+FRONRE
                        !write(2,*)"RDM2-I",RDM_ENERGY+FRONRE
                        !write(3,*)"RDM2-I",RDM_ENERGY+FRONRE
                        joffset=0
                        joffset1=0
                        do j=1,orb%nsub              
            !              RDM_IJ=0.0d0  ! testing
                          do j1=1,orb%occ(j)
                          !write(*,*)"RDM2-J",RDM_ENERGY
                          RDM_KL=0.0d0
                          koffset=0
                          koffset1=0
                          do k=1,orb%nsub; do k1=1,orb%occ(k)
                            !write(*,*)"RDM2-K",RDM_ENERGY
                            loffset=0
                            loffset1=0
                            do l=1,orb%nsub; do l1=1,orb%occ(l)
                              !write(*,*)"RDM2-L-1",RDM_ENERGY
                              !write(*,*)i,i1,j,j1,k,k1,l,l1
                              !write(*,*)ioffset+i1,joffset+j1,koffset+k1,loffset+l1
                              !write(*,*)U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)
                              !write(*,*)mat2%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
                              RDM_KL=RDM_KL+&
                                0.5d0*U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)&
                            *mat2%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
                      !          write(2,*)"================================="
            !                  if(j.eq.0)then
            !                    if(U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1).ne.0)then
            !                      write(12,2211)"U(ijkl)",i,j,k,l,U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1) 
            !                      write(12,2211)"P(iklj)",i,k,l,j,mat2%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
            !                      RDM_J=RDM_J+U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)&
            !                      *mat2%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
            !                    end if
            !                  end if
                        !      dtemp=dtemp+0.5d0*U(ioffset+i1,joffset+j1,koffset+k1,loffset+l1)&
                        !                          *mat2%P(ii,kk,ll,jj)
                                !write(4,*)i,i1,j,j1,k,k1,l,l1
                  end do
                  loffset=loffset+orb%total(l) 
                  loffset1=loffset1+orb%occ(l) 
                  end do
            !                write(*,*)"RDM_KL",RDM_KL
            !                RDM_IJ=RDM_IJ+RDM_KL  
                end do
            !              write(*,*)i,j,"RDM_IJ",RDM_IJ
            !              if(j.eq.3)stop
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%occ(k)
                end do
            ! point group, considered the irreps and group product table 
            !              write(*,*)"RDM_KL is ", RDM_KL
                RDM_ENERGY=RDM_ENERGY+RDM_KL
            !              RDM_IJKL=RDM_IJKL+RDM_KL 
                          !write(*,*)"RDM_KL in trace is (total) ",dtmp2
            !              write(*,*)i,j," RDM_KL ",i1,j1,RDM_KL
                          !stop
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%occ(j)
              end do
              !write(*,*)"2-RDM energy in this orbital : ",dtemp*2.0d0
              !stop 
            end do 
            ioffset=ioffset+orb%total(i)
            ioffset1=ioffset1+orb%occ(i)
            end do
  
          else ! === new code === 
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
          2211      format(A7,I3,I3,I3,I3,G21.12)          

          ! The hD part
            if(old_code)then ! === old code ===
              noffset=0
              noffset1=0
              do n=1,orb%nsub; do n1=1,orb%total(n)
                moffset=0
                moffset1=0
                do m=1,orb%nsub; do m1=1,orb%occ(m)
                  joffset=0
                  joffset1=0   
                  do j=1,orb%nsub; do j1=1,orb%occ(j)
                    mat2%A(n1+noffset,m1+moffset)=mat2%A(n1+noffset,m1+moffset)&
                  +2.0d0*T(n1+noffset,j1+joffset)*mat2%D(m1+moffset1,j1+joffset1)
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
            else ! === new code ===

              ! The hD part -- MKL form
              write(*,*) "dfsghjkdgfsdhjgfhjsdghjdfgjsdhgfhsdgfjsdgfjhgdsjhfgjshdgfjhdgfjhsgdfjhd"
              write(*,*),orb%total
              write(*,*)
              write(*,*),orb%occ
              write(*,*)
              write(*,*),norb
              write(*,*),nocc
              write(*,*),nact
              
              ! noffset=0
              ! do n=1,orb%nsub
              !   n1=orb%total(n)
              !   n2=orb%occ(n)
              !   if(n2.ne.0)then
              !     n0=noffset 
              !     allocate(TM1(N1,N2));TM1=0.0d0
              !     allocate(TM2(N2,N2));TM2=0.0d0
              !     allocate(TM3(N1,N2));TM3=0.0d0
              !     mat1(n)%hD=0.0d0
              !     TM1=T(n0+1:n0+n1,n0+1:n0+n2)
              !     TM2=mat1(n)%D
              !     call MXMG(n1,n2,n2,TM1,TM2,TM3,"NN")
              !     TM3=TM3*2.0d0
              !     mat2%A(n0+1:n0+n1,n0+1:n0+n2)=TM3
              !     deallocate(TM1)
              !     deallocate(TM2)
              !     deallocate(TM3)
              !   end if
              !   noffset=noffset+orb%total(n)
              ! end do

              ! noffset=0
              ! noffset1=0
              ! do n=1,orb%nsub
              !   moffset=0
              !   moffset1=0
              !   do m=1,orb%nsub
              !     joffset=0
              !     joffset1=0
              !     do j=1,orb%nsub
              !       do n1=1,orb%total(n)
              !         do m1=1,orb%occ(m)
              !           do j1=1,orb%occ(j)
              !             mat2%A(n1+noffset,m1+moffset)=mat2%A(n1+noffset,m1+moffset)&
              !             +2.0d0*T(n1+noffset,j1+joffset)*mat2%D(m1+moffset1,j1+joffset1)
              !           end do
              !         end do
              !       end do
              !       joffset=joffset+orb%total(j)
              !       joffset1=joffset1+orb%occ(j)
              !     end do
              !     moffset=moffset+orb%total(m)
              !     moffset1=moffset1+orb%occ(m)
              !   end do
              !   noffset=noffset+orb%total(n)
              !   noffset1=noffset1+orb%occ(n)
              ! end do

              call operator(orb%nsub,mat2%A,T,mat2%D,norb,nocc,orb%total,orb%occ,mat2%p,U,mat2%G)
              
            end if ! ===


              ! The sum_kl[J^(kl)P^(lk)] part         
          
            if(old_code)then ! === old code === 
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
                        mat2%A(n1+noffset,m1+moffset)=mat2%A(n1+noffset,m1+moffset)&
                        +2.0d0*mat2%p(m1+moffset1,k1+koffset1,l1+loffset1,j1+joffset1)&
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

              !          call print_MAT(norb,norb,mat2%A)

              ! The sum_kl[J^(kl)P^(lk)] part --- mkl part,done, will be used later
            else ! === new code ===+

              ! call dshghdg(orb%nsub,mat2%A,T,mat2%D,norb,nocc,orb%total,orb%occ,mat2%p,U)
            
              ! allocate(TM4(norb,norb));TM4=0.0d0
              ! koffset=0
              ! koffset1=0
              ! do k=1,orb%nsub; do k1=1,orb%occ(k)
              !   loffset=0
              !   loffset1=0
              !   do l=1,orb%nsub; do l1=1,orb%occ(l)
              !     noffset=0 
              !     do n=1,orb%nsub
              !       n0=noffset
              !       n1=orb%total(n)
              !       n2=orb%occ(n)
              !       moffset=0
              !       moffset1=0
              !       do m=1,orb%nsub
              !         m0=moffset
              !         m1=orb%total(m)
              !         m2=orb%occ(m)
              !         mx=moffset1
              !         joffset=0
              !         joffset1=0
              !         do j=1,orb%nsub
              !           if(orb%occ(j).ne.0.and.m2.ne.0)then
              !             if(orb%grouptable(n,j).eq.orb%grouptable(k,l))then
              !               j0=joffset
              !               j1=orb%total(j) 
              !               j2=orb%occ(j) 
              !               jx=joffset1
              !               allocate(TM1(N1,J2));TM1=0.0d0 
              !               allocate(TM2(J2,M2));TM2=0.0d0
              !               allocate(TM3(N1,M2));TM3=0.0d0               
              !               TM1=U(n0+1:n0+n1,j0+1:j0+j2,k1+koffset,l1+loffset) 
              !               TM2=mat2%p(mx+1:mx+m2,k1+koffset1,l1+loffset1,jx+1:jx+j2)
              !               call MXMG(n1,m2,j2,TM1,TM2,TM3,"NT") 
              !               TM3=TM3*2.0d0
              !               TM4(N0+1:N0+N1,M0+1:M0+M2)=TM4(N0+1:N0+N1,M0+1:M0+M2)+TM3
              !               deallocate(TM1)
              !               deallocate(TM2)
              !               deallocate(TM3)                   
              !             end if
              !           end if
              !           joffset=joffset+orb%total(j)
              !           joffset1=joffset1+orb%occ(j)
              !         end do 
              !         moffset=moffset+orb%total(m)
              !         moffset1=moffset1+orb%occ(m)
              !       end do
              !       noffset=noffset+orb%total(n)
              !     end do
              !   end do
              !   loffset=loffset+orb%total(l)
              !   loffset1=loffset1+orb%occ(l)
              !   end do
              ! end do 
              ! koffset=koffset+orb%total(k)
              ! koffset1=koffset1+orb%occ(k)
              ! end do
              ! mat2%A=mat2%A+TM4
              ! deallocate(TM4)

              ! noffset=0
              ! noffset1=0
              ! do n=1,orb%nsub
              !   moffset=0
              !   moffset1=0
              !   do m=1,orb%nsub
              !     joffset=0
              !     joffset1=0
              !     do j=1,orb%nsub
              !       koffset=0
              !       koffset1=0
              !       do k=1,orb%nsub
              !         loffset=0
              !         loffset1=0
              !         do l=1,orb%nsub
              !           do n1=1,orb%total(n)
              !             do m1=1,orb%occ(m)
              !               do j1=1,orb%occ(j)
              !                 do k1=1,orb%occ(k)
              !                   do l1=1,orb%occ(l)
              !                     mat2%A(n1+noffset,m1+moffset)=mat2%A(n1+noffset,m1+moffset)&
              !                     +2.0d0*mat2%p(m1+moffset1,k1+koffset1,l1+loffset1,j1+joffset1)&
              !                     *U(n1+noffset,j1+joffset,k1+koffset,l1+loffset)
              !                   end do
              !                 end do
              !               end do
              !             end do
              !           end do
              !           loffset=loffset+orb%total(l)
              !           loffset1=loffset1+orb%occ(l)
              !         end do
              !         koffset=koffset+orb%total(k)
              !         koffset1=koffset1+orb%occ(k)
              !       end do
              !       joffset=joffset+orb%total(j)
              !       joffset1=joffset1+orb%occ(j)
              !     end do
              !     moffset=moffset+orb%total(m)
              !     moffset1=moffset1+orb%occ(m)
              !   end do
              !   noffset=noffset+orb%total(n)
              !   noffset1=noffset1+orb%occ(n)
              ! end do

            end if ! ===


          ! ioffset=0 
          ! do i=1,orb%nsub
          !   do k=1,orb%total(i)
          !     do l=1,orb%total(i)
          !       if(dabs(mat2%A(k+ioffset,l+ioffset)).lt.1.0e-9)then
          !          mat2%A(k+ioffset,l+ioffset)=0.0d0
          !       end if
          !       mat1(i)%A(k,l)=mat2%A(k+ioffset,l+ioffset)
          !     end do
          !   end do
          !   ioffset=ioffset+orb%total(i)
          ! end do

          ! do i=1,norb
          !   do j=1,norb
          !     if(dabs(mat2%A(i,j)).lt.1.0d-9) mat2%A(i,j)=0.0d0
          !     if(mat2%A(i,j)-mat2%A(j,i).gt.1.0d-6) then
          !       write(1,*)i,j,"i,j",mat2%A(i,j)-mat2%A(j,i)
          !     end if
          !   end do 
          ! end do

goto 1333 

          T1=0.0d0
          T2=0.0d0
          ioffset=0
          ioffset1=0
          do i=1,orb%nsub
            i1=orb%occ(i)
            koffset=0
            koffset1=0
            do k=1,orb%nsub; do k1=1,orb%occ(i)
              loffset=0
              loffset1=0
              do l=1,orb%nsub; do l1=1,orb%occ(l)
              if(i1.gt.0)then 
              call MXM(i1,U(1+ioffset:i1+ioffset,1+ioffset:i1+ioffset,k1+koffset,l1+loffset),&
                 mat2%p(l1+loffset1,k1+koffset1,1+ioffset1:i1+ioffset1,1+ioffset1:i1+ioffset1),&
                        T1(1+ioffset:i1+ioffset,1+ioffset:i1+ioffset))
                    T2(1+ioffset:i1+ioffset,1+ioffset:i1+ioffset)&
                   =T2(1+ioffset:i1+ioffset,1+ioffset:i1+ioffset)&
                   +T1(1+ioffset:i1+ioffset,1+ioffset:i1+ioffset) 
              end if  
              !write(*,*)k,l,"K,l"  
              end do
              loffset=loffset+orb%total(l)
              loffset1=loffset1+orb%occ(l)
              end do
            end do
            koffset=koffset+orb%total(k)
            koffset1=koffset1+orb%occ(k) 
            end do
            ioffset=ioffset+orb%total(i)
            ioffset1=ioffset1+orb%occ(i)
          end do
          
1333       continue

          if(old_code)then ! === old code ===
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
                do k=1,orb%nsub; do k1=1,orb%total(k)
                  nr=n1+noffset1
                  mr=m1+moffset1
                  jr=j1+joffset1
                  kr=k1+koffset1
                  ni=n1+noffset
                  mi=m1+moffset
                  ji=j1+joffset
                  ki=k1+koffset 
                  mat2%G(ni,mi,ji,ki)&
                    =mat2%G(ni,mi,ji,ki)+2.0d0*T(ni,ki)*mat2%D(mr,jr)
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

            !  ==============

            ! The G=hD+sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)] --  MKL & point-group
            ! only hD part, i.e. h*D_ij
          else ! === new code ===

            ! allocate(GM1(norb,norb,norb,norb)) ! the k,ij,l | The i,j can be further reduced
            ! GM1=0.0d0
            ! ioffset=0
            ! ioffset1=0
            ! do i=1,orb%nsub; do ii=1,orb%occ(i)
            !   i0=ioffset
            !   i1=orb%total(i) 
            !   i2=orb%occ(i) 
            !   ix=ioffset1
    
            !   joffset=0
            !   joffset1=0
            !   do j=1,orb%nsub; do jj=1,orb%occ(j)             
            !     j0=joffset
            !     j1=orb%total(j) 
            !     j2=orb%occ(j) 
            !     jx=joffset1
              
            !     GM1(:,i0+ii,j0+jj,:)=T*mat2%D(ix+ii,jx+jj)*2.0d0
            !   end do
            !   joffset=joffset+orb%total(j)
            !   joffset1=joffset1+orb%occ(j)          
            !   end do 
            ! end do
            ! ioffset=ioffset+orb%total(i)
            ! ioffset1=ioffset1+orb%occ(i)          
            ! end do
            ! mat2%G=GM1
            ! deallocate(GM1)

            ! noffset=0
            ! noffset1=0
            ! do n=1,orb%nsub
            !   moffset=0
            !   moffset1=0 
            !   do m=1,orb%nsub
            !     joffset=0
            !     joffset1=0
            !     do j=1,orb%nsub
            !       koffset=0
            !       koffset1=0
            !       do k=1,orb%nsub
            !         do n1=1,orb%total(n)
            !           do m1=1,orb%occ(m)
            !             do j1=1,orb%occ(j)
            !               do k1=1,orb%total(k)
            !                 nr=n1+noffset1
            !                 mr=m1+moffset1
            !                 jr=j1+joffset1
            !                 kr=k1+koffset1
            !                 ni=n1+noffset
            !                 mi=m1+moffset
            !                 ji=j1+joffset
            !                 ki=k1+koffset 
            !                 mat2%G(ni,mi,ji,ki)&
            !                   =mat2%G(ni,mi,ji,ki)+2.0d0*T(ni,ki)*mat2%D(mr,jr)
            !               end do
            !             end do
            !           end do
            !         end do
            !       koffset=koffset+orb%total(k)
            !       koffset1=koffset1+orb%occ(k)
            !       end do
            !       joffset=joffset+orb%total(j)
            !       joffset1=joffset1+orb%occ(j)
            !     end do
            !     moffset=moffset+orb%total(m)
            !     moffset1=moffset1+orb%occ(m)
            !   end do
            !   noffset=noffset+orb%total(n)
            !   noffset1=noffset1+orb%occ(n) 
            ! end do


            ! noffset=0
            ! noffset1=0
            ! do n=1,orb%nsub; do n1=1,orb%total(n)
            !   moffset=0
            !   moffset1=0 
            !   do m=1,orb%nsub; do m1=1,orb%occ(m)
            !     joffset=0
            !     joffset1=0
            !     do j=1,orb%nsub; do j1=1,orb%occ(j)
            !       koffset=0
            !       koffset1=0 
            !       do k=1,orb%nsub; do k1=1,orb%total(k)
            !         nr=n1+noffset1
            !         mr=m1+moffset1
            !         jr=j1+joffset1
            !         kr=k1+koffset1
            !         ni=n1+noffset
            !         mi=m1+moffset
            !         ji=j1+joffset
            !         ki=k1+koffset 
            !         mat2%G(ni,mi,ji,ki)&
            !           =mat2%G(ni,mi,ji,ki)+2.0d0*T(ni,ki)*mat2%D(mr,jr)
            !         end do
            !       koffset=koffset+orb%total(k)
            !       koffset1=koffset1+orb%occ(k) 
            !       end do
            !     end do
            !     joffset=joffset+orb%total(j)
            !     joffset1=joffset1+orb%occ(j)
            !     end do
            !   end do
            !   moffset=moffset+orb%total(m)
            !   moffset1=moffset1+orb%occ(m)
            !   end do
            ! end do
            ! noffset=noffset+orb%total(n)
            ! noffset1=noffset1+orb%occ(n) 
            ! end do 

          end if ! ===
          if(old_code)then ! === old code ===
          noffset=0
          noffset1=0
          do n=1,orb%nsub; do n1=1,orb%total(n)
            moffset=0
            moffset1=0 
            do m=1,orb%nsub; do m1=1,orb%occ(m)
              kooffset=0
              kooffset1=0
              do ko=1,orb%nsub; do ko1=1,orb%total(ko)
                joffset=0
                joffset1=0
                do j=1,orb%nsub; do j1=1,orb%occ(j)
                  koffset=0
                  koffset1=0
                  do k=1,orb%nsub; do k1=1,orb%occ(k)
                    loffset=0
                    loffset1=0
                    do l=1,orb%nsub; do l1=1,orb%occ(l)

                      nr=n1+noffset1
                      mr=m1+moffset1
                      kor=ko1+kooffset1
                      jr=j1+joffset1
                      kr=k1+koffset1
                      lr=l1+loffset1

                      ni=n1+noffset
                      mi=m1+moffset
                      koi=ko1+kooffset
                      ji=j1+joffset
                      ki=k1+koffset
                      li=l1+loffset

                      mat2%G(ni,mi,ji,koi)=mat2%G(ni,mi,ji,koi)&
                      +2.0d0*mat2%P(mr,kr,lr,jr)*U(ni,koi,ki,li)&
                +2.0d0*2.0d0*mat2%P(mr,jr,lr,kr)*U(ni,ki,li,koi)
                    ! write(201,*)ni,mi,ji,koi,mat2%G(ni,mi,ji,koi)


                    end do
                    loffset=loffset+orb%total(l)
                    loffset1=loffset1+orb%occ(l)
                    end do
                  end do
                  koffset=koffset+orb%total(k)
                  koffset1=koffset1+orb%occ(k)
                  end do
                !write(*,*)ni,mi,ji,koi,mat2%G(ni,mi,ji,koi)
                !stop
                end do
                joffset=joffset+orb%total(j)
                joffset1=joffset1+orb%occ(j)
                end do
              end do
              kooffset=kooffset+orb%total(ko)
              kooffset1=kooffset1+orb%occ(ko)
              end do
            end do
            moffset=moffset+orb%total(m)
            moffset1=moffset1+orb%occ(m)
            end do 
          end do
          noffset=noffset+orb%total(n)
          noffset1=noffset1+orb%occ(n)
          end do

            !      G=hD+"sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"
            !      sum_(kl)[J^(kl)P^(lk)+2K^(kl)Q^(lk)]"  part ---- mkl version  
                  ! stage-1 :  matrix form ready, partly point-group 
                  ! stage-2 :  dgemm and fully point group (later) since it will be too time-consuming
          else ! === new code ===      
  
          ! allocate(GM1(norb,norb,norb,norb))
          ! GM1=0.0d0

          ! ioffset=0
          ! ioffset1=0
          ! do i=1,orb%nsub; do ii=1,orb%occ(i)
          !   i0=ioffset
          !   i1=orb%total(i)
          !   i2=orb%occ(i)
          !   ix=ioffset1

          !   joffset=0
          !   joffset1=0 
          !   do j=1,orb%nsub; do jj=1,orb%occ(j)
          !     j0=joffset
          !     j1=orb%total(j)
          !     j2=orb%occ(j)
          !     jx=joffset1
      
          !     allocate(TM1(norb,norb));TM1=0.0d0

          !     koffset=0
          !     koffset1=0
          !     do k=1,orb%nsub; do kk=1,orb%occ(k)
          !       k0=koffset
          !       k1=orb%total(k)
          !       k2=orb%occ(k)
          !       kx=koffset1

          !       loffset=0
          !       loffset1=0
          !       do l=1,orb%nsub

          !         if(orb%grouptable(i,j).eq.orb%grouptable(k,l))then
          !           do ll=1,orb%occ(l)
          !             l0=loffset
          !             l1=orb%total(l)
          !             l2=orb%occ(l)
          !             lx=loffset1
          !             dtmp=mat2%P(ix+ii,kx+kk,lx+ll,jx+jj)

          !             TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp*2.0d0 
          !           end do
          !         end if

          !         if(orb%grouptable(i,k).eq.orb%grouptable(j,l))then
          !           do ll=1,orb%occ(l)
          !             l0=loffset
          !             l1=orb%total(l)
          !             l2=orb%occ(l)
          !             lx=loffset1
          !             dtmp=mat2%P(ix+ii,jx+jj,lx+ll,kx+kk)
          !             TM1=TM1+U(:,k0+kk,l0+ll,:)*dtmp*4.0d0
          !           end do
          !         end if
          !         loffset=loffset+orb%total(l)
          !         loffset1=loffset1+orb%occ(l)
          !       end do
          !     end do
          !     koffset=koffset+orb%total(k)
          !     koffset1=koffset1+orb%occ(k)
          !     end do
          !     GM1(:,i0+ii,j0+jj,:)=TM1 
          !     deallocate(TM1) 

          !   end do
          !   joffset=joffset+orb%total(j)
          !   joffset1=joffset1+orb%occ(j)
          !   end do 
          ! end do
          ! ioffset=ioffset+orb%total(i)
          ! ioffset1=ioffset1+orb%occ(i)
          ! end do         

          ! mat2%G=mat2%G+GM1
          ! deallocate(GM1)
            ! walltime(20) = wtime()
            
            ! noffset=0
            ! noffset1=0
            ! do n=1,orb%nsub; do n1=1,orb%total(n)
            !   moffset=0
            !   moffset1=0 
            !   do m=1,orb%nsub; do m1=1,orb%occ(m)
            !     kooffset=0
            !     kooffset1=0
            !     do ko=1,orb%nsub; do ko1=1,orb%total(ko)
            !       joffset=0
            !       joffset1=0
            !       do j=1,orb%nsub; do j1=1,orb%occ(j)
            !         koffset=0
            !         koffset1=0
            !         do k=1,orb%nsub; do k1=1,orb%occ(k)
            !           loffset=0
            !           loffset1=0
            !           do l=1,orb%nsub; do l1=1,orb%occ(l)

            !             nr=n1+noffset1
            !             mr=m1+moffset1
            !             kor=ko1+kooffset1
            !             jr=j1+joffset1
            !             kr=k1+koffset1
            !             lr=l1+loffset1

            !             ni=n1+noffset
            !             mi=m1+moffset
            !             koi=ko1+kooffset
            !             ji=j1+joffset
            !             ki=k1+koffset
            !             li=l1+loffset

            !             mat2%G(ni,mi,ji,koi)=mat2%G(ni,mi,ji,koi)&
            !             +2.0d0*mat2%P(mr,kr,lr,jr)*U(ni,koi,ki,li)&
            !       +2.0d0*2.0d0*mat2%P(mr,jr,lr,kr)*U(ni,ki,li,koi)
            !           ! write(201,*)ni,mi,ji,koi,mat2%G(ni,mi,ji,koi)


            !           end do
            !           loffset=loffset+orb%total(l)
            !           loffset1=loffset1+orb%occ(l)
            !           end do
            !         end do
            !         koffset=koffset+orb%total(k)
            !         koffset1=koffset1+orb%occ(k)
            !         end do
            !       !write(*,*)ni,mi,ji,koi,mat2%G(ni,mi,ji,koi)
            !       !stop
            !       end do
            !       joffset=joffset+orb%total(j)
            !       joffset1=joffset1+orb%occ(j)
            !       end do
            !     end do
            !     kooffset=kooffset+orb%total(ko)
            !     kooffset1=kooffset1+orb%occ(ko)
            !     end do
            !   end do
            !   moffset=moffset+orb%total(m)
            !   moffset1=moffset1+orb%occ(m)
            !   end do 
            ! end do
            ! noffset=noffset+orb%total(n)
            ! noffset1=noffset1+orb%occ(n)
            ! end do

            ! noffset=0
            ! noffset1=0
            ! do n=1,orb%nsub
            !   moffset=0
            !   moffset1=0 
            !   do m=1,orb%nsub
            !     kooffset=0
            !     kooffset1=0
            !     do ko=1,orb%nsub
            !       joffset=0
            !       joffset1=0
            !       do j=1,orb%nsub
            !         koffset=0
            !         koffset1=0
            !         do k=1,orb%nsub
            !           loffset=0
            !           loffset1=0
            !           do l=1,orb%nsub
            !             do n1=1,orb%total(n)
            !               do m1=1,orb%occ(m)               
            !                 do ko1=1,orb%total(ko)
            !                   do j1=1,orb%occ(j)
            !                     do k1=1,orb%occ(k)
            !                       do l1=1,orb%occ(l)
            !                         nr=n1+noffset1
            !                         mr=m1+moffset1
            !                         kor=ko1+kooffset1
            !                         jr=j1+joffset1
            !                         kr=k1+koffset1
            !                         lr=l1+loffset1

            !                         ni=n1+noffset
            !                         mi=m1+moffset
            !                         koi=ko1+kooffset
            !                         ji=j1+joffset
            !                         ki=k1+koffset
            !                         li=l1+loffset

            !                         mat2%G(ni,mi,ji,koi)=mat2%G(ni,mi,ji,koi)&
            !                         +2.0d0*mat2%P(mr,kr,lr,jr)*U(ni,koi,ki,li)&
            !                         +2.0d0*2.0d0*mat2%P(mr,jr,lr,kr)*U(ni,ki,li,koi)
            !                       end do
            !                     end do
            !                   end do
            !                 end do
            !               end do
            !             end do
            !             loffset=loffset+orb%total(l)
            !             loffset1=loffset1+orb%occ(l)
            !           end do
            !           koffset=koffset+orb%total(k)
            !           koffset1=koffset1+orb%occ(k)
            !         end do
            !         joffset=joffset+orb%total(j)
            !         joffset1=joffset1+orb%occ(j)
            !       end do
            !       kooffset=kooffset+orb%total(ko)
            !       kooffset1=kooffset1+orb%occ(ko)
            !     end do
            !     moffset=moffset+orb%total(m)
            !     moffset1=moffset1+orb%occ(m)
            !   end do
            !   noffset=noffset+orb%total(n)
            !   noffset1=noffset1+orb%occ(n)
            ! end do
            ! walltime(21) = wtime()
            ! write(*,*)"*******************G_2_time*****",walltime(21)-walltime(20)
          end if ! ===
 
          if(method(icycle).eq."microit")then
            call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag)
          end if

          ! Too lazy, actually only diagional should be enough in this case (also microit)
          if(method(icycle).eq."cp-inte")then
            call Hessian(norb,mat2%A,mat2%G,mat2%H,mat2%Hdiag) 
          end if

         !call Hessian_test(norb,A,mat2%G,mat2%H,mat2%Hdiag)
         !stop

        else  ! This is a very old version!!!!! need to be rewritten 

          !----------closed shell HF-----------
          !          do i=1,orb%nsub
          !            do j=1,orb%closed(i)
          !              RDM_ENERGY=RDM_ENERGY+2*T(ioffset+j,ioffset+j)
          !              koffset=0
          !              do k=1,orb%nsub
          !                do l=1,orb%closed(k)
          !                  RDM_ENERGY=RDM_ENERGY+2*U(ioffset+j,ioffset+j,koffset+l,koffset+l)&
          !                             -U(ioffset+j,koffset+l,ioffset+j,koffset+l)
          !                end do
          !                koffset=koffset+orb%total(k)
          !              end do
          !            end do
          !            ioffset=ioffset+orb%total(i)
          !          end do

          ! ---------------------------------------
          ! With closed shell approximation
          ! ----------  active CI   ---------------

          old_code=.true.   ! temporary active it 

          if(old_code)then

            allocate(Fc(nocc,nocc))  
            Fc=0.0d0
            ! First we need to calculate the closed shell Fock matrix            
            ioffset=0
            do i=1,orb%nsub  
              koffset=0 
              do k=1,orb%nsub
                if(orb%grouptable(i,orb%grouptable(i,orb%grouptable(k,k))).eq.1)then
                  do i1=1,orb%occ(i)
                    ix=i1+ioffset
                    do j1=1,orb%occ(i)
                      jx=j1+ioffset
                      Fc(ix,jx)=T(ix,jx)
                      do k1=1,orb%closed(k)
                        kx=k1+koffset 
                        Fc(ix,jx)=Fc(ix,jx)+2.0d0*U(ix,jx,kx,kx)&
                                           -1.0d0*U(ix,kx,kx,jx)
                      end do
                    end do
                  end do
                end if
                koffset=koffset+orb%total(k)
              end do 
              ioffset=ioffset+orb%total(i)
            end do
              
            !            call print_mat(nocc,nocc,Fc)    
            !            call print_mat(norb,norb,T)    
            !            stop               

            ! The sum of occ <i|h+Fc|i>
            ioffset=0
            do i=1,orb%nsub;  
              do i1=1,orb%closed(i)
                ix=i1+ioffset
                RDM_ENERGY=RDM_ENERGY+T(ix,ix)+Fc(ix,ix)
              end do
              ioffset=ioffset+orb%total(i)
            end do            

            ! + the sum of val <i|Fc*D|i>
            ioffset=0
            do i=1,orb%nsub              
              do i1=1,orb%act(i)
                ix=i1+ioffset+orb%closed(i)
                do j1=1,orb%act(i)
                  jx=j1+ioffset+orb%closed(i)
                  RDM_ENERGY=RDM_ENERGY+Fc(ix,jx)*mat1(i)%D(i1,j1) 
                end do
              end do
              ioffset=ioffset+orb%total(i)
            end do

            write(*,*)RDM_ENERGY
            call flush(6)

            ! + the sum of cal <i|J*P|j>
            ioffset=0
            ioffset1=0
            do i=1,orb%nsub
             joffset=0
             joffset1=0
             do j=1,orb%nsub
              koffset=0
              koffset1=0
              do k=1,orb%nsub
               loffset=0
               loffset1=0
               do l=1,orb%nsub
                 if(orb%grouptable(i,orb%grouptable(j,orb%grouptable(k,l))).eq.1)then

                   do i1=1,orb%act(i)                    
                    ix=i1+ioffset+orb%closed(i) 
                    io=i1+ioffset1
                    do j1=1,orb%act(j)
                     jx=j1+joffset+orb%closed(j) 
                     jo=j1+joffset1
                     do k1=1,orb%act(k)
                      kx=k1+koffset+orb%closed(k) 
                      ko=k1+koffset1
                      do l1=1,orb%act(l)
                       lx=l1+loffset+orb%closed(l) 
                       lo=l1+loffset1
                        RDM_ENERGY=RDM_ENERGY+0.5d0*U(ix,jx,kx,lx)&
                                              *mat2%P(io,ko,lo,jo)
                      end do
                     end do
                    end do
                   end do

                 end if 
                 loffset=loffset+orb%total(l)
                 loffset1=loffset1+orb%act(l)
               end do
               koffset=koffset+orb%total(k)
               koffset1=koffset1+orb%act(k)
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%act(j)
             end do
             ioffset=ioffset+orb%total(i)
             ioffset1=ioffset1+orb%act(i)
            end do

            !            ioffset=0
            !            do i=1,orb%nsub
            !              do j=1,orb%closed(i)
            !                koffset=0
            !                do k=1,orb%nsub
            !                  do l=1,orb%act(k)
            !                    do m=1,orb%act(k)
            !                      RDM_ENERGY=RDM_ENERGY&
            !                                 +(2*U(koffset+l+orb%closed(k)&
            !                                 ,koffset+m+orb%closed(k),ioffset+j,ioffset+j)&
            !                                 -U(koffset+l+orb%closed(k),ioffset+j&
            !                                 ,koffset+m+orb%closed(k),ioffset+j))&
            !                                 *mat1(k)%D(l,m)
            !                    end do
            !                  end do
            !                  koffset=koffset+orb%total(k)
            !                end do
            !              end do
            !              ioffset=ioffset+orb%total(i)
            !            end do
              
            !            ioffset=0
            !            ioffset1=0
            !            do i=1,orb%nsub; do i1=1,orb%act(i)
            !              joffset=0
            !              joffset1=0
            !              do j=1,orb%nsub; do j1=1,orb%act(j)
            !                koffset=0
            !                koffset1=0
            !                do k=1,orb%nsub; do k1=1,orb%act(k)
            !                  loffset=0
            !                  loffset1=0
            !                  do l=1,orb%nsub; do l1=1,orb%act(l)
            !                    RDM_ENERGY=RDM_ENERGY&
            !                               +0.5d0*U(ioffset+i1+orb%closed(i),joffset+j1+orb%closed(j)&
            !                               ,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l))&
            !                               *mat2%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
            !                  end do
            !                  loffset=loffset+orb%total(l)
            !                  loffset1=loffset1+orb%act(l)
            !                  end do
            !                end do
            !                koffset=koffset+orb%total(k)
            !                koffset1=koffset1+orb%act(k)
            !                end do
            !              end do
            !              joffset=joffset+orb%total(j)
            !              joffset1=joffset1+orb%act(j)
            !              end do
            !            end do
            !            ioffset=ioffset+orb%total(i)
            !            ioffset1=ioffset1+orb%act(i)
            !            end do

            ! write(*,2514,advance="no")"RDM_ENERGY",RDM_ENERGY+FRONRE0
            

          else 
            write(*,*)"Now it should be in closed shell part"  
          end if
 
          Stop

            !!----------Gc--------------
            Gc=0.0d0
            ioffset=0
            do i=1,orb%nsub; do i1=1,orb%total(i)
              joffset=0
              do j=1,orb%nsub; do j1=1,orb%total(j)
              if(i==j) then
                Gc(ioffset+i1,joffset+j1)=Gc(ioffset+i1,joffset+j1)&
                      +2*T(ioffset+i1,joffset+j1)
 
                koffset=0
                do k=1,orb%nsub; do k1=1,orb%closed(k)
                  Gc(ioffset+i1,joffset+j1)=Gc(ioffset+i1,joffset+j1)&
                       +4*U(ioffset+i1,joffset+j1,koffset+k1,koffset+k1)&
                       -2*U(ioffset+i1,koffset+k1,joffset+j1,koffset+k1)
                end do
                koffset=koffset+orb%total(k)
                end do
 
                loffset=0
                do l=1,orb%nsub
                  do m=1,orb%act(l); do n=1,orb%act(l)
                    Gc(ioffset+i1,joffset+j1)=Gc(ioffset+i1,joffset+j1)&  
                    +(2*U(ioffset+i1,joffset+j1,loffset+m+orb%closed(l),loffset+n+orb%closed(l))&
                    -U(ioffset+i1,loffset+m+orb%closed(l),joffset+j1,loffset+n+orb%closed(l)))&
                    *mat1(l)%D(m,n)
                  end do
                  end do
                  loffset=loffset+orb%total(l)
                end do
 
              end if
              end do
              joffset=joffset+orb%total(j)
              end do
            end do
            ioffset=ioffset+orb%total(i)
            end do
  
          !!--------------A------------------

            ioffset=0        
            do i=1,orb%nsub; do i1=1,orb%total(i)
              joffset=0
              do j=1,orb%nsub; do j1=1,orb%closed(j)
                mat2%A(ioffset+i1,joffset+j1)=Gc(ioffset+i1,joffset+j1)
              end do
              joffset=joffset+orb%total(j)
              end do
            end do
            ioffset=ioffset+orb%total(i)
            end do
 
            ioffset=0
            do i=1,orb%nsub; do i1=1,orb%total(i)
              joffset=0
              joffset1=0
              do j=1,orb%nsub; do j1=1,orb%act(j)
              if(i==j) then
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%act(k)
                  if(j==k) then
                    mat2%A(ioffset+i1,joffset+j1+orb%closed(j))=mat2%A(ioffset+i1,joffset+j1+orb%closed(j))&
                         +T(ioffset+i1,koffset+k1+orb%closed(k))*mat1(k)%D(j1,k1)
     
                    loffset=0
                    do l=1,orb%nsub; do l1=1,orb%closed(l)
                      mat2%A(ioffset+i1,joffset+j1+orb%closed(j))=mat2%A(ioffset+i1,joffset+j1+orb%closed(j))&
                      +(2*U(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1,loffset+l1)&
                      -U(ioffset+i1,loffset+l1,koffset+k1+orb%closed(k),loffset+l1))&
                      *mat1(k)%D(j1,k1)
                    end do
                    loffset=loffset+orb%total(l)
                    end do
                  end if
                  loffset=0
                  loffset1=0
                  do l=1,orb%nsub; do l1=1,orb%act(l)
                    moffset=0
                    moffset1=0
                    do m=1,orb%nsub; do m1=1,orb%act(m)
                      mat2%A(ioffset+i1,joffset+j1+orb%closed(j))=mat2%A(ioffset+i1,joffset+j1+orb%closed(j))& 
                                +U(ioffset+i1,koffset+k1+orb%closed(k)&
                                ,loffset+l1+orb%closed(l),moffset+m1+orb%closed(m))&
                                *mat2%P(joffset1+j1,loffset1+l1,moffset1+m1,koffset1+k1)
                    end do
                    moffset=moffset+orb%total(m)
                    moffset1=moffset1+orb%act(m)
                    end do
                  end do
                  loffset=loffset+orb%total(l)
                  loffset1=loffset1+orb%act(l)
                  end do    
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%act(k)
                end do
              end if
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%act(j)
              end do
            end do
            ioffset=ioffset+orb%total(i)
            end do
 
            do i=1,norb
              do j=1,norb
                mat2%A(i,j)=2*mat2%A(i,j)
          
              end do
            end do
          !!sstop
 
          !!--------------G--------------
            ioffset=0
            do i=1,orb%nsub; do i1=1,orb%total(i)
              joffset=0
              do j=1,orb%nsub; do j1=1,orb%total(j)
                koffset=0
                koffset1=0
                do k=1,orb%nsub; do k1=1,orb%act(k)
                  loffset=0
                  loffset1=0
                  do l=1,orb%nsub; do l1=1,orb%act(l)
                    if(k==l) then
                      mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                      =mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                      +T(ioffset+i1,joffset+j1)*mat1(k)%D(k1,l1)
  
                      moffset=0
                      do m=1,orb%nsub; do m1=1,orb%closed(m)
                        mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                        =mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                        +(2*U(ioffset+i1,joffset+j1,moffset+m1,moffset+m1)&
                        -U(ioffset+i1,moffset+m1,joffset+m1,moffset+m1))&
                        *mat1(k)%D(k1,l1)
                      end do
                      moffset=moffset+orb%total(m)
                      end do
                    end if
     
                    moffset=0
                    moffset1=0
                    do m=1,orb%nsub; do m1=1,orb%act(m)
                      noffset=0
                      noffset1=0
                      do n=1,orb%nsub; do n1=1,orb%act(n)
                        mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                        =mat2%G(ioffset+i1,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l),joffset+j1)&
                        +U(ioffset+i1,joffset+j1,moffset+m1+orb%closed(m),noffset+n1+orb%closed(n))&
                        *mat2%P(koffset1+k1,moffset1+m1,noffset1+n1,loffset1+l1)&
                        +2*U(ioffset+i1,moffset+m1+orb%closed(m),joffset+j1,noffset+n1+orb%closed(n))&
                        *mat2%P(koffset1+k1,loffset1+l1,noffset1+n1,moffset1+m1)
                      end do
                      noffset=noffset+orb%total(n)
                      noffset1=noffset1+orb%act(n)
                      end do
                    end do
                    moffset=moffset+orb%total(m)
                    moffset1=moffset1+orb%act(m)
                    end do
      
                  end do
                  loffset=loffset+orb%total(l)
                  loffset1=loffset1+orb%act(l)
                  end do
                end do
                koffset=koffset+orb%total(k)
                koffset1=koffset1+orb%act(k)
                end do
              end do
              joffset=joffset+orb%total(j)
              joffset1=joffset1+orb%act(j)
              end do
            end do
            ioffset=ioffset+orb%total(i)
            end do
 
            do i=1,norb
              do j=1,norb
                do k=1,norb
                  do l=1,norb
                    mat2%G(i,j,k,l)=2*mat2%G(i,j,k,l)
          !!                   write(11,*) i,j,k,l,mat2%G(i,j,k,l)
                  end do
                end do
              end do
            end do 
 
          !!-----------H D------------
            H_temp=0.0d0
            Hdiag_temp=0.0d0 
            if(method(icycle).eq."microit")then
              call Hessian(norb,mat2%A,mat2%G,H_temp,Hdiag_temp)
 
              ioffset=0
              do i=1,orb%nsub; do i1=1,orb%act(i)
                joffset=0
                do j=1,orb%nsub; do j1=1,orb%act(j)
                  koffset=0
                  do k=1,orb%nsub; do k1=1,orb%act(k)
                    loffset=0
                    do l=1,orb%nsub; do l1=1,orb%act(l)
                      mat2%H(ioffset+i1+orb%closed(i),joffset+j1+orb%closed(j)&
                      ,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l))=&
                      H_temp(ioffset+i1+orb%closed(i),joffset+j1+orb%closed(j)&
                      ,koffset+k1+orb%closed(k),loffset+l1+orb%closed(l))
                    end do
                    loffset=loffset+orb%total(l)
                    end do
                  end do
                  koffset=koffset+orb%total(k)
                  end do
                  mat2%Hdiag(ioffset+i1+orb%closed(i),joffset+j1+orb%closed(j))=&
                  Hdiag_temp(ioffset+i1+orb%closed(i),joffset+j1+orb%closed(j))
                end do
                joffset=joffset+orb%total(j)
                end do
              end do
              ioffset=ioffset+orb%total(i)
              end do
    
              ioffset=0
              do i=1,orb%nsub; do i1=1,orb%act(i)
                joffset=0
                do j=1,orb%nsub; do j1=1,orb%closed(j)
                  mat2%Hdiag(ioffset+i1+orb%closed(i),joffset+j1)=2*(mat2%G(joffset+j1,ioffset+i1+orb%closed(i)&
                  ,ioffset+i1+orb%closed(i),joffset+j1)+Gc(ioffset+i1+orb%closed(i),ioffset+i1+orb%closed(i))&
                  +6*U(ioffset+i1+orb%closed(i),joffset+j1,ioffset+i1+orb%closed(i),joffset+j1)&
                  -2*U(ioffset+i1+orb%closed(i),ioffset+i1+orb%closed(i),joffset+j1,joffset+j1))&
                  -mat2%A(ioffset+i1+orb%closed(i),ioffset+i1+orb%closed(i))-mat2%A(joffset+j1,joffset+j1)
    
                  koffset=0
                  do k=1,orb%nsub; do k1=1,orb%act(k)
                    if(i==k) then
                      mat2%Hdiag(ioffset+i1+orb%closed(i),joffset+j1)=mat2%Hdiag(ioffset+i1+orb%closed(i),joffset+j1)&
                      -4*(3*U(ioffset+i1+orb%closed(i),joffset+j1,koffset+k1+orb%closed(k),joffset+j1)&
                      -U(ioffset+i1+orb%closed(i),koffset+k1+orb%closed(k),joffset+j1,joffset+j1))&
                      *mat1(i)%D(i1,k1)
                    end if
                  end do
                  koffset=koffset+orb%total(k)
                  end do
                end do
                joffset=joffset+orb%total(j)
                end do
              end do
              ioffset=ioffset+orb%total(i)
              end do
            end if 

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

      end subroutine PreOneSCF

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


