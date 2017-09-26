      Subroutine hfenergy

        use global_control 

!        write(*,*)orb%ele 

        double precision ONE_ELE_ENERGY
        double precision TWO_ELE_ENERGY
 
        integer,allocatable::occ(:)

        allocate(occ(orb%nsub)); occ=0
        occ=(orb%ele)/2
    
        ONE_ELE_ENERGY=0.0d0
        TWO_ELE_ENERGY=0.0d0 

        ioffset=0
        do i=1,orb%nsub
          do j=1,occ(i)
            ONE_ELE_ENERGY=ONE_ELE_ENERGY+2.0d0*T(ioffset+j,ioffset+j)
          end do  
          ioffset=ioffset+orb%total(i)
        end do 

        ioffset=0
        do i=1,orb%nsub
          do j=1,occ(i)
            koffset=0
            do k=1,orb%nsub
              do l=1,occ(k)
                TWO_ELE_ENERGY=TWO_ELE_ENERGY &
                        +2.0d0*U(j+ioffset,j+ioffset,l+koffset,l+koffset)&
                        -1.0d0*U(j+ioffset,l+koffset,l+koffset,j+ioffset)
!                write(*,*)j+ioffset,l+koffset
!                write(*,*)U(j+ioffset,j+ioffset,l+koffset,l+koffset),&
!                          U(j+ioffset,l+koffset,l+koffset,j+ioffset)
              end do
              koffset=koffset+orb%total(k) 
            end do
          end do 
          ioffset=ioffset+orb%total(i)
        end do
       
!        if(ifcore2)then
        HF_ENERGY=FRONRE0+ONE_ELE_ENERGY+TWO_ELE_ENERGY   
!        else 
!          HF_ENERGY=FRONRE+ONE_ELE_ENERGY+TWO_ELE_ENERGY    
!        end if       
        write(*,*)""
        write(*,*)"RHF STATE PART:"
        write(*,*)" One-electron   energy :",ONE_ELE_ENERGY
        write(*,*)" Two-electron   energy :",TWO_ELE_ENERGY
!        if(ifcore2)then
        write(*,*)" Core & Nuclear energy :",FRONRE0
!        else
!          write(*,*)" Core & Nuclear energy :",FRONRE
!        end if
        write(*,*)" [- RHF state energy ] :",HF_ENERGY
!        Stop

      end subroutine hfenergy

