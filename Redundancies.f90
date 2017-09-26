      Subroutine detect_Redundancies()

        use global_control
        use matrix

        call detect_degeneracy()

      End Subroutine detect_Redundancies

      ! detect degenracy basing on 1-integrals
      Subroutine detect_degeneracy()

        use global_control
        use matrix

        double precision dv
        integer,allocatable::pair(:,:)

        allocate(redu%irreps(orb%nsub))
        allocate(redu%reflct(orb%nsub))
        allocate(redu%offset(orb%nsub))
        allocate(redu%orbirp(norb))
        allocate(redu%sign_mat(norb,norb))
        ij=0
        do i=1,orb%nsub
          redu%irreps(i)=i
          do j=1,orb%total(i)
            ij=ij+1 
            if(j.le.orb%occ(i))then
              redu%orbirp(ij)=-i
            else
              redu%orbirp(ij)= i
            end if
            if(j.le.orb%closed(i))then
              redu%orbirp(ij)=-i-10
            end if
          end do
        end do
        redu%sign_mat=0 

        write(1,*)"In detect_degeneracy"
! Detect all suspicious degeneracies
        is0=0
        do isym=1,orb%nsub
          ks0=0
          do ksym=1,orb%nsub
            ! Same total orbitals but different irreps
            if(orb%total(isym).eq.orb%total(ksym).and.isym.lt.ksym)then

              do i=1,orb%total(isym)
                ii=is0+i
                kk=ks0+i
                if(T(ii,ii).eq.T(kk,kk))then
!                  write(*,*)"degeneracies",ii,kk,T(ii,kk)  
                  redu%irreps(ksym)=redu%irreps(isym)
                end if
              end do

              if(dabs(T(ii,ii)-T(kk,kk)).lt.thrs%s)then
                do i=1,orb%total(isym)
                  do j=1,orb%total(isym)
                    ii=is0+i;  jj=is0+j
                    kk=ks0+i;  ll=ks0+j
!                    write(*,*)ii,jj,T(ii,jj),kk,ll,T(kk,ll)
                    if(dabs(T(ii,jj)+T(kk,ll)).lt.thrs%s)then
                      dv=-1.0d0
                    else
                      dv= 1.0d0
                    end if
                    redu%sign_mat(kk,ll)=dv
                  end do           
                end do
              end if 

            end if
            ks0=ks0+orb%total(ksym)            
          end do 
          is0=is0+orb%total(isym)
        end do
!        call print_mat(norb,norb,T,1)
!        write(*,*)"T // Sign"
        do i=1,norb
          do j=1,norb
            write(1,"(f3.0)",advance='no')redu%sign_mat(i,j)
          end do
          write(1,*)
        end do
!        call print_mat(norb,norb,redu%sign_mat,1)
!        write(1,*)redu%irreps

        ndegenrate_pair=0
        do i=1,orb%nsub
          do j=1,orb%nsub
            if(redu%irreps(i).eq.redu%irreps(j).and.i.lt.j)then
              ndegenrate_pair=ndegenrate_pair+1
            end if
          end do
        end do
        np = ndegenrate_pair ! save the degenrate pair

        allocate(pair(2,np))
        ip=0
        do i=1,orb%nsub
          do j=1,orb%nsub
            if(redu%irreps(i).eq.redu%irreps(j).and.i.lt.j)then
              ip=ip+1
              pair(1,ip)=i
              pair(2,ip)=j
            end if
          end do
        end do

        ioffset=0
        do i=1,orb%nsub
          redu%reflct(i)=i
          redu%offset(i)=ioffset
          ioffset=ioffset+orb%occ(i)           
        end do
        do i=1,np
          redu%reflct(pair(1,i))=pair(2,i) 
          redu%reflct(pair(2,i))=pair(1,i) 
        end do

        write(1,*)"redu%irreps",redu%irreps
        write(1,*)"redu%reflct",redu%reflct       
        write(1,*)"redu%offset",redu%offset 
        write(1,*)"redu%orbirp",redu%orbirp 

        deallocate(pair)
!        stop 

      end subroutine


      Subroutine symmetrize()
! Symmetrize the one- and two-RDMs

        use global_control
        use matrix

        double precision,allocatable::TM1(:,:),TM2(:,:),TM3(:,:)
!        double precision,allocatable::GM1(:,:,:,:)

        type::average
          double precision,allocatable::D(:,:)
        end type average
        type(average),allocatable::TM(:),TN(:)

        double precision dv1,dv2,dv3

!        double precision DM1(nact**2,nact**2) 
!        allocate(GM1(nact,nact,nact,nact)) 

!        call print_mat(nact,nact,mat2%D,1)

! For 1-RDMs
        allocate(TM(orb%nsub)) 
 
        is0=0
        do isym=1,orb%nsub
          allocate(TM(isym)%D(orb%act(isym),orb%act(isym)))
          TM(isym)%D=0.0d0
          TM(isym)%D=&
          mat2%D(is0+1:is0+orb%act(isym),is0+1:is0+orb%act(isym)) 
          js0=0
          do jsym=1,isym-1
            if(redu%irreps(isym).eq.redu%irreps(jsym))then
! return as order and averaged
              if(.true.)then
                allocate(TM1(orb%act(isym),orb%act(isym)))
                TM1=0.0d0
                TM1=(dabs(TM(isym)%D)+dabs(TM(jsym)%D))/2.0d0
                do i=1,orb%act(isym)
                  do j=1,orb%act(isym)
                    TM(isym)%D(i,j) = DSIGN(TM1(i,j),TM(isym)%D(i,j))
                    TM(jsym)%D(i,j) = DSIGN(TM1(i,j),TM(jsym)%D(i,j)) 
                  end do
                end do
                mat2%D(is0+1:is0+orb%act(isym),is0+1:is0+orb%act(isym))&
                =TM(isym)%D
                mat2%D(js0+1:js0+orb%act(jsym),js0+1:js0+orb%act(jsym))&
                =TM(jsym)%D
                deallocate(TM1)
! return as reflect
              else
                allocate(TM1(orb%act(isym),orb%act(isym)))
                allocate(TM2(orb%act(isym),orb%act(isym)))
                TM1=0.0d0
                TM2=0.0d0
                TM1=dabs(TM(isym)%D)
                TM2=dabs(TM(jsym)%D)
                do i=1,orb%act(isym)
                  do j=1,orb%act(isym)
                    TM(isym)%D(i,j) = DSIGN(TM2(i,j),TM(isym)%D(i,j))
                    TM(jsym)%D(i,j) = DSIGN(TM1(i,j),TM(jsym)%D(i,j))
                  end do
                end do
                mat2%D(is0+1:is0+orb%act(isym),is0+1:is0+orb%act(isym))&
                =TM(isym)%D
                mat2%D(js0+1:js0+orb%act(jsym),js0+1:js0+orb%act(jsym))&
                =TM(jsym)%D
                deallocate(TM1)
                deallocate(TM2)
              end if 
            end if
            js0=js0+orb%act(jsym)
          end do          
          is0=is0+orb%act(isym)
        end do

        deallocate(TM)

!        call print_mat(nact,nact,mat2%D,1)
!        call print_gat(nact,nact,nact,nact,mat2%P,11)

! For the 2-RDMs
        do is=1,orb%nsub
          do js=1,orb%nsub
            do ks=1,orb%nsub
              do ls=1,orb%nsub

                if(orb%grouptable(is,orb%grouptable&
                                 (js,orb%grouptable(ks,ls))).eq.1)then

                  do i=1,orb%act(is) 
                    do j=1,orb%act(js) 
                      do k=1,orb%act(ks) 
                        do l=1,orb%act(ls) 
                          ! the starting index
                          ii=i+redu%offset(is) 
                          jj=j+redu%offset(js) 
                          kk=k+redu%offset(ks) 
                          ll=l+redu%offset(ls) 
                          ! the reflected index 
                          ix=i+redu%offset(redu%reflct(is)) 
                          jx=j+redu%offset(redu%reflct(js)) 
                          kx=k+redu%offset(redu%reflct(ks)) 
                          lx=l+redu%offset(redu%reflct(ls)) 

                          dv1=mat2%P(ii,kk,ll,jj)
                          dv2=mat2%P(ix,kx,lx,jx)
   
                          dv3=(dabs(dv1)+dabs(dv2))/2.0d0
! return as order and averaged
                          if(.true.)then
                            mat2%P(ii,kk,ll,jj)=DSIGN(dv3,dv1)
                            mat2%P(ix,kx,lx,jx)=DSIGN(dv3,dv2)
! return as reflect
                          else
                            mat2%P(ii,kk,ll,jj)=DSIGN(dabs(dv2),dv1)
                            mat2%P(ix,kx,lx,jx)=DSIGN(dabs(dv1),dv2)
                          end if 
 
                        end do
                      end do 
                    end do
                  end do
                 
                end if 

              end do
            end do
          end do
        end do 

        write(1,*)""
        write(1,*)""
        write(1,*)""
!        call print_gat(nact,nact,nact,nact,mat2%P,12)

!        stop 

      end Subroutine symmetrize




