      Subroutine transform_INT(if_trans)

        use global_control 
        use matrix

        logical::if_trans
        
        character*10 A2,A3,A4 
        integer,allocatable::base(:),IL(:) 
        double precision d1,d2
        
   
        allocate(base(orb%nsub))

        ioffset=0
        do i=1,orb%nsub
          base(i)=ioffset
          ioffset=ioffset+orb%total(i)
        end do         

        ! This is the old whole orbital space transform
        !        call transform_RAM(norb,mat2%U,T,U)
              
        ! The is the subspace integrals transform
        if(if_trans)then
          call transform_subspaces&
             (orb%nsub,orb%occ,orb%total,norb,MAT2%U,T,U,orb%grouptable)
        end if

        !        call print_mat(norb,norb,MAT2%U)    
        !        stop 
                !write(*,*)T
        open(unit=111,file="FCIDUMP_NEW")
        open(unit=121,file=FCIDUMP)
          read(121,*)A2,A3,A4
          if (A3(1:5).ne."NORB=") then
            A4=""
            A4=trim(A3(6:))
            read(A4,*) norb
          else
            read(A4,*) norb
          end if
          write(111,"(1X,A4,1X,A5,I4)")A2,A3,norb 
          allocate(IL(norb))                         !define  
          !write(*,*)norb
          read(121,"(A9)",advance='no')A2
          write(111,"(A9)",advance='no')A2
          read(121,*)(IL(i),i=1,norb)
          do i=1,norb
            write(111,"(I2,A1)",advance='no')IL(i),","
          end do
          write(111,*)
          read(121,"(A10)")A3
          write(111,"(A10)")A3
          write(111,*)"&END"
        close(121)
 
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
              do i=1,orb%total(i0)
                do j=1,i
                  ij=ij+1
                  kl=0
                  do k=1,orb%total(k0)
                    do l=1,k
                      kl=kl+1
                      i2=i+base(i0)
                      j2=j+base(j0)
                      k2=k+base(k0)
                      l2=l+base(l0)
                      d2=U(i2,j2,k2,l2)
                      if(i0.eq.k0)then
                        if(ij.ge.kl)then
                          write(111,112)d2,i2,j2,k2,l2
                        end if
                      else
                        write(111,112)d2,i2,j2,k2,l2
                      end if
                    end do
                  end do
                end do
              end do
            else 
              ij=0
              do i=1,orb%total(i0)
                do j=1,orb%total(j0)
                  ij=ij+1
                  kl=0
                  do k=1,orb%total(k0)
                    do l=1,orb%total(l0)
                      kl=kl+1
                      i2=i+base(i0)
                      j2=j+base(j0)
                      k2=k+base(k0)
                      l2=l+base(l0)
                      d2=U(i2,j2,k2,l2)
                      if(i0.eq.k0)then
                        if(ij.ge.kl)then
                          write(111,112)d2,i2,j2,k2,l2
                        end if
                      else
                        write(111,112)d2,i2,j2,k2,l2
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do 
        close(113)

        ioffset=0
        do i=1,orb%nsub
          do i1=1,orb%total(i)
            do i2=1,i1
              write(111,112)T(ioffset+i1,ioffset+i2),ioffset+i1,ioffset+i2,0,0
            end do
          end do
          ioffset=ioffset+orb%total(i)
        end do

        write(111,112)FRONRE0,0,0,0,0

        close(111)
 
        112     format(G21.12,I4,I4,I4,I4)

        deallocate(base)

        call integrals_active("FCIDUMP_NEW_ACTIVE",.false.)

      end subroutine transform_INT


