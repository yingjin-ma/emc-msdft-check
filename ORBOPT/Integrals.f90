      Subroutine integrals_active(fcidump1,if_2nd_trans)

! Write integrals of the active space into FILE (e.g. FCIDUMP)

        use global_control
        use matrix 

        character*18::fcidump1  
        
        character*10 A1
        character*10 A2
        character*10 A3

! Only available when active the 2-index transformation of Integrals        
        logical::if_2nd_trans

        integer,allocatable::base(:),base1(:),base2(:)
        double precision d1,d2

        double precision,allocatable::PartT(:,:)
        !fcidump1=""

        A1=''
        A2=''
        A3='' 

!        write(6,*)"integrals_active",nact  
!        call flush(6) 

!        call print_gat(nact,nact,nact,nact,U2nd,6)
!        call flush(6) 

        allocate(base(orb%nsub))
        allocate(base1(orb%nsub))
        allocate(base2(orb%nsub))
        base=0
        base1=0
        base2=0

        ioffset=0
        ioffset1=0
        ioffset2=0
        do i=1,orb%nsub
          base(i)=ioffset
          base1(i)=ioffset1
          base2(i)=ioffset2
          ioffset=ioffset+orb%total(i)
          ioffset1=ioffset1+orb%occ(i)!-orb%closed(i)
          ioffset2=ioffset2+orb%occ(i)-orb%closed(i)
        end do
 
        !write(*,*)orb%closed
!        write(*,*)fcidump0,fcidump1
!        write(*,*)"ifcore, ifcore2", ifcore, ifcore2

        if(ifcore2.eqv..false.)then
          open(unit=117,file=fcidump1)
            write(117,"(A12,I3,A11)")" &FCI NORB= ",nact,",NELEC=  2," 
            write(117,"(A9)",advance='no')"  ORBSYM=" 
            do i=1,orb%nsub
              do j=1,orb%occ(i)
                write(117,"(I2,A1)",advance='no')i,","
              end do
            end do
            write(117,*)
            write(117,"(A8)")"  ISYM=1" 
            write(117,"(A5)")" &END"
            if(orb%sym.eq."d2h")then
             open(unit=113,file="INTEGRAL_ORDER_8_IRREPS.grp") 
            else if(orb%sym.eq."c2v")then ! also C2h, D2
             open(unit=113,file="INTEGRAL_ORDER_4_IRREPS.grp")
            else if(orb%sym.eq."Cs")then  ! also C2, Ci
             open(unit=113,file="INTEGRAL_ORDER_2_IRREPS.grp")
            else if(orb%sym.eq."c1")then
             open(unit=113,file="INTEGRAL_ORDER_1_IRREPS.grp")
            end if
!            write(*,*)"Start FCIDUMP operations"   
              do 
                read(113,*,iostat=ierr)i0,j0,k0,l0
                if(ierr.ne.0)exit

                if(i0.eq.j0.and.k0.eq.l0)then
                  ij=0
                  do i=1,orb%occ(i0)
                    do j=1,i
                      ij=ij+1
                      kl=0
                      do k=1,orb%occ(k0)
                        do l=1,k
                          kl=kl+1
                          i2=i+base(i0)
                          j2=j+base(j0)
                          k2=k+base(k0)
                          l2=l+base(l0)
                          i1=i+base1(i0)
                          j1=j+base1(j0)
                          k1=k+base1(k0)
                          l1=l+base1(l0) 
                          if(if_2nd_trans)then
                            d2=U2nd(i1,j1,k1,l1)
                          else
                            d2=U(i2,j2,k2,l2)
                          end if
                          if(i0.eq.k0)then
                            if(ij.ge.kl)then
                              write(117,112)d2,i1,j1,k1,l1
                            end if
                          else
                            write(117,112)d2,i1,j1,k1,l1
                          end if
                        end do
                      end do
                    end do
                  end do
                else 
                  ij=0
                  do i=1,orb%occ(i0)
                    do j=1,orb%occ(j0)
                      ij=ij+1
                      kl=0
                      do k=1,orb%occ(k0)
                        do l=1,orb%occ(l0)
                          kl=kl+1
                          i2=i+base(i0)
                          j2=j+base(j0)
                          k2=k+base(k0)
                          l2=l+base(l0)
                          i1=i+base1(i0)
                          j1=j+base1(j0)
                          k1=k+base1(k0)
                          l1=l+base1(l0)
                          if(if_2nd_trans)then
                            d2=U2nd(i1,j1,k1,l1)
                          else
                            d2=U(i2,j2,k2,l2)
                          end if
                          if(i0.eq.k0)then
                            if(ij.ge.kl)then
                              write(117,112)d2,i1,j1,k1,l1
                            end if
                          else
                            write(117,112)d2,i1,j1,k1,l1
                          end if
                        end do
                      end do
                    end do
                  end do
                end if
              end do 
            close(113)

            ioffset=0
            ioffset1=0
            do i=1,orb%nsub
              do i1=1,orb%occ(i)
                do i2=1,i1
                  if(if_2nd_trans)then
                    write(117,112)T2nd(ioffset1+i1,ioffset1+i2),ioffset1+i1,ioffset1+i2,0,0
                  else
                    write(117,112)   T(ioffset+i1,ioffset+i2),ioffset1+i1,ioffset1+i2,0,0
                  end if
                end do
              end do
              ioffset1=ioffset1+orb%occ(i)
              ioffset=ioffset+orb%total(i)
            end do
            write(117,112)FRONRE0,0,0,0,0

            FRONREA=FRONRE0
!            write(*,*)"fcidump_active exactly active ,FRONREA",FRONREA 

          close(117)

!          write(*,*)"nocore"
!          call flush(6) 

        else

!         write(6,*)" core/closed shell "
!         call flush(6) 
!          write(6,*)"orb%frozen",orb%frozen
!          write(6,*)"orb%closed",orb%closed
!          write(6,*)"orb%act",orb%act
!          write(6,*)"orb%occ",orb%occ

          open(unit=117,file=fcidump1)
            write(117,"(A12,I3,A11)")" &FCI NORB= ",nact,&
                                     ",NELEC=  2," 
            write(117,"(A9)",advance='no')"  ORBSYM=" 
            do i=1,orb%nsub
              do j=1,orb%occ(i)-orb%closed(i)
                write(117,"(I2,A1)",advance='no')i,","
              end do
            end do
            write(117,*)
            write(117,"(A8)")"  ISYM=1" 
            write(117,"(A5)")" &END"
            if(orb%sym.eq."d2h")then
             open(unit=113,file="INTEGRAL_ORDER_8_IRREPS.grp") 
            else if(orb%sym.eq."c2v")then
             open(unit=113,file="INTEGRAL_ORDER_4_IRREPS.grp")
            else if(orb%sym.eq."Cs")then
             open(unit=113,file="INTEGRAL_ORDER_2_IRREPS.grp")
            else if(orb%sym.eq."c1")then
             open(unit=113,file="INTEGRAL_ORDER_1_IRREPS.grp")
            end if   
              do 
                read(113,*,iostat=ierr)i0,j0,k0,l0
                if(ierr.ne.0)exit

                if(i0.eq.j0.and.k0.eq.l0)then
                  ij=0
                  do i=1+orb%closed(i0),orb%occ(i0)
                    do j=1+orb%closed(i0),i
                      ij=ij+1
                      kl=0
                      do k=1+orb%closed(k0),orb%occ(k0)
                        do l=1+orb%closed(k0),k
                          kl=kl+1
                          i2=i+base(i0)
                          j2=j+base(j0)
                          k2=k+base(k0)
                          l2=l+base(l0)
                          i1=i+base1(i0)!-orb%closed(i0)
                          j1=j+base1(j0)!-orb%closed(j0)
                          k1=k+base1(k0)!-orb%closed(k0)
                          l1=l+base1(l0)!-orb%closed(l0)
                          i3=i+base2(i0)-orb%closed(i0)
                          j3=j+base2(j0)-orb%closed(j0)
                          k3=k+base2(k0)-orb%closed(k0)
                          l3=l+base2(l0)-orb%closed(l0)
                          if(if_2nd_trans)then
                            d2=U2nd(i1,j1,k1,l1) ! old
                            !d2=U2nd(i2,j2,k2,l2) 
                          else
                            d2=   U(i2,j2,k2,l2)
                          end if
                          if(i0.eq.k0)then
                            if(ij.ge.kl)then
                              write(117,112)d2,i3,j3,k3,l3
                              !write(6,*)d2,i2,j2,k2,l2
                              !call flush(6) 
                            end if
                          else
                            write(117,112)d2,i3,j3,k3,l3
                          end if
                        end do
                      end do
                    end do
                  end do
                else 
                  ij=0
                  do i=1+orb%closed(i0),orb%occ(i0)
                    do j=1+orb%closed(j0),orb%occ(j0)
                      ij=ij+1
                      kl=0
                      do k=1+orb%closed(k0),orb%occ(k0)
                        do l=1+orb%closed(l0),orb%occ(l0)
                          kl=kl+1
                          i2=i+base(i0)
                          j2=j+base(j0)
                          k2=k+base(k0)
                          l2=l+base(l0)
                          i1=i+base1(i0)!-orb%closed(i0)
                          j1=j+base1(j0)!-orb%closed(j0)
                          k1=k+base1(k0)!-orb%closed(k0)
                          l1=l+base1(l0)!-orb%closed(l0)
                          i3=i+base2(i0)-orb%closed(i0)
                          j3=j+base2(j0)-orb%closed(j0)
                          k3=k+base2(k0)-orb%closed(k0)
                          l3=l+base2(l0)-orb%closed(l0)
                          if(if_2nd_trans)then 
                            d2=U2nd(i1,j1,k1,l1)              
                          else
                            d2=   U(i2,j2,k2,l2)
                          end if
                          if(i0.eq.k0)then
                            if(ij.ge.kl)then
                              write(117,112)d2,i3,j3,k3,l3
                            end if
                          else
                            write(117,112)d2,i3,j3,k3,l3
                          end if
                        end do
                      end do
                    end do
                  end do
                end if
              end do 
            close(113)

            allocate(PartT(nact,nact))
            PartT=0.0d0

            
        ioffset=0
        ioffset1=0
        do i=1,orb%nsub
          do j=1+orb%closed(i),orb%occ(i)
            koffset=0
            koffset1=0
            do k=1,orb%nsub
              do l=1+orb%closed(k),orb%occ(k)

                ! h 
                !ij1=j+ioffset1!-orb%closed(i)
                !kl1=l+koffset1!-orb%closed(k)
                ij1=j+base2(i)-orb%closed(i)
                kl1=l+base2(k)-orb%closed(k)
                if(if_2nd_trans)then 
                  PartT(ij1,kl1)=T2nd(j+ioffset1,l+koffset1)  
                else 
                  PartT(ij1,kl1)=   T(j+ioffset,l+koffset)
                end if 

                ! 2J-K
                moffset=0
                moffset1=0
                !write(*,*)PartT(ij1,kl1),ij1-orb%closed(i),kl1-orb%closed(k)
                do m=1,orb%nsub
                  do n=1,orb%closed(m) 
                    if(if_2nd_trans)then 
                      PartT(ij1,kl1)=PartT(ij1,kl1)&
                          +2.0d0*U2nd(j+ioffset1,l+koffset1,n+moffset1,n+moffset1)&
                          -1.0d0*U2nd(j+ioffset1,n+moffset1,n+moffset1,l+koffset1)
                    else
                      PartT(ij1,kl1)=PartT(ij1,kl1)&
                             +2.0d0*U(j+ioffset,l+koffset,n+moffset,n+moffset)&
                             -1.0d0*U(j+ioffset,n+moffset,n+moffset,l+koffset)
                    end if 
                  end do
                  moffset=moffset+orb%total(m)
                  moffset1=moffset1+orb%occ(m)
                end do  
              end do
              koffset=koffset+orb%total(k) 
              !koffset1=koffset1+orb%occ(k)-orb%closed(k) 
              koffset1=koffset1+orb%occ(k)
            end do
          end do 
          ioffset=ioffset+orb%total(i)
          !ioffset1=ioffset1+orb%occ(i)-orb%closed(i)
          ioffset1=ioffset1+orb%occ(i)
        end do
      
            !do i=1,nact
            !  do j=1,nact
            !    write(*,*)i,j,PartT(i,j)
            !  end do
            !end do 

            ioffset=0
            ioffset1=0
            do i=1,orb%nsub
              do i1=1+orb%closed(i),orb%occ(i)
                do i2=1+orb%closed(i),i1
                  ii1=base2(i)+i1-orb%closed(i)
                  ii2=base2(i)+i2-orb%closed(i)
                  write(117,112)PartT(ii1,ii2),&
                  base2(i)+i1-orb%closed(i),base2(i)+i2-orb%closed(i),0,0
                end do
              end do
              !ioffset1=ioffset1+orb%occ(i)-orb%closed(i)
              ioffset1=ioffset1+orb%occ(i)
              ioffset=ioffset+orb%total(i)
            end do

            FRONRE=FRONRE0
            ioffset=0
            ioffset1=0
            do i=1,orb%nsub
              do j=1,orb%closed(i)
                if(if_2nd_trans)then 
                  FRONRE=FRONRE+2.0d0*T2nd(j+ioffset1,j+ioffset1)
                else
                  FRONRE=FRONRE+2.0d0*T(j+ioffset,j+ioffset)
                end if 
                koffset=0
                koffset1=0
                do k=1,orb%nsub
                  do l=1,orb%closed(k)
                    if(if_2nd_trans)then 
             FRONRE=FRONRE+2.0d0*U2nd(j+ioffset1,j+ioffset1,l+koffset1,l+koffset1)&
                          -1.0d0*U2nd(j+ioffset1,l+koffset1,l+koffset1,j+ioffset1)
                    else
             FRONRE=FRONRE+2.0d0*U(j+ioffset,j+ioffset,l+koffset,l+koffset)&
                          -1.0d0*U(j+ioffset,l+koffset,l+koffset,j+ioffset)
                    end if
                  end do
                  koffset=koffset+orb%total(k)
                  koffset1=koffset1+orb%occ(k)
                end do 
              end do
              ioffset=ioffset+orb%total(i)
              ioffset1=ioffset1+orb%occ(i)
            end do  

            write(117,112)FRONRE,0,0,0,0
            FRONREA=FRONRE
!            write(*,*)"fcidump_active else ,FRONREA",FRONREA 

          close(117)
          !stop

        end if 

112     format(G21.12,I4,I4,I4,I4)

        call flush(117) 
 
!        write(*,*)"fcidump_active before finish,FRONREA",FRONREA 
!        call flush(6)
!        stop

      end Subroutine integrals_active

