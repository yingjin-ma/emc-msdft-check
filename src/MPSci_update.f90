      Subroutine mpsci_update(icycle)

        use global_control
        use coupling_terms
        use matrix


        call coupling_prepare()

        !call MPSci_gradient()
        if(.true.)then 
          ! Current read in the full Hamiltonian  
          call coupling_CICI() 
        else
          write(6,*)"It should only read in diagional elements"
        end if
        !stop
        call MPSci_perturb()

        call coupling_finish()  

      end Subroutine mpsci_update

! 1st order perturbation update
      Subroutine MPSci_perturb()

        use global_control
        use coupling_terms
        use matrix

        double precision HamiCK(nC,nC)
        double precision eigCK(nC)

! Used in the 2 particle grad part
        double precision tmp_KL, E_reformed 
        double precision E2,T1,T2,V1,V2
        double precision Heig(nC) ! diagional elements of H
        double precision rdm1check(nact,nact)

        double precision updated_Davidson_c(nC_all)

        E2=0.0d0
        E_reformed=0.0d0

        write(2,*)"before MPSci perturbation"
        call flush(2) 

        do iC=1,nC
          Ci(iC)%gi=0.0d0

! This is actually too fussy!! 
!   RDM1 -> RDM1%sym
          j0=0; k0=0
          rdm1check=0.0d0
          do i=1,orb%nsub
            mat1(i)%D=0.0d0
            !write(6,*)"i th iter", i
            call flush(2)
            do j=1,orb%act(i)
              do k=1,orb%act(i)
!               write(6,*)j,k
!               call flush(6)
                mat1(i)%D(j,k)=Ci(iC)%d(j0+j,k0+k)               
              end do
            end do
!           rdm1check(1+j0:orb%act(i)+j0,1+j0:orb%act(i)+j0)=mat1(i)%D 
            j0=j0+orb%act(i)
            k0=k0+orb%act(i)
          end do       
!          write(2,*)"RDM-deri for", iC,"-th elements"
!          call print_mat(nact,nact,rdm1check,2)

! MPSci gradient from one particle part
          ioffset=0
          do i=1,orb%nsub
            do j=1,orb%occ(i)
              koffset=0
              do k=1,orb%nsub
                if(i.eq.k)then
                  do l=1,orb%occ(k)
                    Ci(iC)%gi=Ci(iC)%gi+T2nd(ioffset+j,koffset+l)*mat1(i)%D(j,l)
                  end do
                end if
                koffset=koffset+orb%act(k)
              end do
            end do
            ioffset=ioffset+orb%act(i)
          end do 

! MPSci gradient from two particle part
          ioffset1=0
!         RDM_IJKL=0.0d0
          do i=1,orb%nsub; do i1=1,orb%occ(i)
            joffset1=0
            do j=1,orb%nsub; do j1=1,orb%occ(j)
              tmp_KL=0.0d0
              koffset1=0
              do k=1,orb%nsub; do k1=1,orb%occ(k)
                loffset1=0
                do l=1,orb%nsub; do l1=1,orb%occ(l)
                  tmp_KL=tmp_KL+&
                    0.5d0*U2nd(ioffset1+i1,joffset1+j1,koffset1+k1,loffset1+l1)&
              *Ci(iC)%P(ioffset1+i1,koffset1+k1,loffset1+l1,joffset1+j1)
                end do
                loffset1=loffset1+orb%occ(l)
                end do
              end do
              koffset1=koffset1+orb%occ(k)
              end do
              Ci(iC)%gi=Ci(iC)%gi+tmp_KL
            end do
            joffset1=joffset1+orb%occ(j)
            end do
          end do
          ioffset1=ioffset1+orb%occ(i)
          end do

          !Ci(iC)%gi = Ci(iC)%gi - FRONRE0
          !Ci(iC)%gi = Ci(iC)%gi - FRONRE0
!          write(6,*)iC,"Ci(iC)%gi", Ci(iC)%gi
!          call flush(6)

        end do

        V1=0 
        do iC=1,nC
          E_reformed=E_reformed+Ci(iC)%gi*Ci(iC)%dv
          V1=V1+Ci(iC)%dv**2
        end do

        E2=E_reformed/V1
        write(2,*)"E2 is",E2+FRONRE0
        call flush(2)
 
! Read in the diagional elements of Hami
        do iC=1,nC 
          ! add soon 
          Heig(iC)=Hami(iC,iC) 
          !write(2,*)"Diagional matrix ",ic,Heig(iC)
        end do
!        write(2,*)"The Hamiltonian matrix (only doagional is used)"
!        call print_maT(nC,nC,Hami,2)
!        HamiCK=Hami
!        call eigtql2(nC,HamiCK,eigCK)
!        write(2,*)"The eigen-value for check "
!        eigCK=eigCK+FRONRE0
!        call print_mat(1,nC,eigCK,2)  

! Perturbation the MPSci vector
        do iC=1,nC
          if(iC.eq.1)then
            Ci(iC)%dvp=Ci(iC)%dv
          else 
            V1 = Ci(iC)%gi - E2*Ci(iC)%dv
            V2 = Heig(i) - E2
            Ci(iC)%dvp = Ci(iC)%dv - V1/V2
            write(2,*)"grad, Ci(iC)%dv, Ci(iC)%dvp", V1, Ci(iC)%dv, Ci(iC)%dvp
          end if
        end do

        call MODIFIED_GRAM_SCHMIDT(nC,1,Ci(:)%dvp) 

!        write(*,*)Ci(:)%dvp
 
        updated_Davidson_c=0.0d0
        do iC=1,nC 
          updated_Davidson_c(Ci(iC)%num)=Ci(iC)%dvp
        end do

        !write(*,*)"updated_Davidson_c"
        !write(*,*) updated_Davidson_c
        open(unit=100,file="updated_davidson.vec")
          do iC=1,nC_all       
            write(100,*)1.0d0*updated_Davidson_c(iC)
            write(2,*)ic,1.0d0*updated_Davidson_c(iC)
          end do 
        close(100)

      end Subroutine MPSci_perturb


